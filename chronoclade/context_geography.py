"""Offline country composition of frozen, annotated context catalogues.

No dates, tree tips, or live services participate in the public denominator.
"""
from __future__ import annotations

import base64
import csv
import hashlib
import html
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

DEPTHS = (5, 6, 7)
CAUTION = (
    'Geographical composition of available sequenced records, affected by surveillance, '
    'submission coverage and selection. These figures do not estimate country prevalence, '
    'incidence, migration direction or transmission. Country metadata is not country of acquisition.'
)
FIELDS = ('source', 'snapshot', 'scope', 'filters', 'scheme', 'scheme_version', 'prefix_depth',
          'prefix_key', 'group_label', 'assignment_category', 'assignment_status', 'cohort',
          'country', 'count', 'denominator', 'known_country_n', 'unknown_n', 'percentage',
          'counting_unit', 'raw_record_count')


def country_colour(country: str) -> str:
    """Stable across runs, cohorts and depths; Unknown and Other remain distinct."""
    if country == 'Unknown':
        return '#9b9b9b'
    if country == 'Other':
        return '#303030'
    import colorsys
    digest = hashlib.sha256(country.encode()).digest()
    hue = int.from_bytes(digest[:4], 'big') / 2**32
    rgb = colorsys.hsv_to_rgb(hue, .48 + digest[4] / 255 * .22, .65 + digest[5] / 255 * .2)
    return '#' + ''.join(f'{round(channel * 255):02x}' for channel in rgb)


def _source_id(row: Mapping[str, Any]) -> str:
    return str(row.get('source_genome_id') or row.get('genome_id') or row.get('id') or '')


def _aliases(row: Mapping[str, Any]) -> set[str]:
    aliases = set()
    for field, kind in (('biosample', 'biosample'), ('biosample_accession', 'biosample'),
                        ('assembly_accession', 'assembly'), ('assembly', 'assembly')):
        value = row.get(field)
        if value:
            aliases.add(f'{kind}:{str(value).strip().upper()}')
    for field, kind in (('assembly_accessions', 'assembly'), ('run_accessions', 'run')):
        for value in row.get(field, []) or []:
            aliases.add(f'{kind}:{str(value).strip().upper()}')
    for value in row.get('aliases', []) or []:
        aliases.add('accession:' + str(value).strip().upper())
    if row.get('sample_unit_id'):
        aliases.add('unit:' + str(row['sample_unit_id']))
    if _source_id(row):
        aliases.add('source:' + str(row.get('source', 'pathogenwatch')) + ':' + _source_id(row))
    return aliases


def _eligible(row: Mapping[str, Any]) -> bool:
    if 'qc_pass' in row:
        return row['qc_pass'] is True
    if 'qc_passed' in row:
        return row['qc_passed'] is True
    return str(row.get('qc_status', '')).lower() in {'pass', 'passed'}


def _deduplicate(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Connected accession aliases handle assembly chains without counting them twice."""
    parent = list(range(len(rows)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen = {}
    for i, row in enumerate(rows):
        for alias in _aliases(row):
            if alias in seen:
                parent[find(i)] = find(seen[alias])
            seen[alias] = i
    clusters = defaultdict(list)
    for i, row in enumerate(rows):
        clusters[find(i)].append(row)
    units, audit = [], []
    for members in clusters.values():
        members.sort(key=lambda r: (str(r.get('source', 'pathogenwatch')), _source_id(r)))
        unit = dict(members[0])
        aliases = sorted(set().union(*(_aliases(r) for r in members)))
        biosamples = [a for a in aliases if a.startswith('biosample:')]
        unit['_unit'] = (biosamples or aliases or ['unlinked:' + str(len(units))])[0]
        unit['_aliases'] = aliases
        unit['_source_ids'] = sorted({identifier for r in members for identifier in
            (r.get('source_genome_ids') or [_source_id(r)])})
        unit['_raw_n'] = sum(int(r.get('raw_genome_count', 1)) for r in members)
        countries = {str(r.get('country') or 'Unknown') for r in members}
        if len(countries) > 1:
            unit['country'] = 'Unknown'
        # Conflicting assignments must not quietly borrow the first record's code.
        for depth in DEPTHS:
            keys = {r.get(f'cglin_group_{depth}') for r in members
                    if r.get(f'cglin_group_{depth}')}
            if len(keys) > 1:
                unit[f'cglin_group_{depth}'] = ''
                unit[f'cglin_status_{depth}'] = 'conflict'
        for member in members:
            audit.append({'sample_unit': unit['_unit'], 'source_genome_id': _source_id(member),
                          'representative_source_id': _source_id(unit),
                          'sample_id': member.get('sample_id', ''),
                          'raw_country': member.get('country_raw', member.get('country', '')),
                          'normalised_country': member.get('country') or 'Unknown',
                          'country_conflict': len(countries) > 1,
                          'biosample_conflict': len(biosamples) > 1,
                          'identity_resolved': bool(biosamples) or bool(unit.get('identity_resolved')),
                          'raw_records_in_unit': unit['_raw_n']})
        units.append(unit)
    return sorted(units, key=lambda r: r['_unit']), audit


def geography_tables(catalogue_rows: Iterable[Mapping[str, Any]], *,
                     selected_source_ids: Iterable[str] = (),
                     focal_rows: Iterable[Mapping[str, Any]] = (),
                     depths: Iterable[int] = DEPTHS, scope: Mapping[str, Any] | None = None) -> dict:
    """Return auditable full catalogue, selected and separate focal composition tables."""
    scope = dict(scope or {})
    rows = [dict(r) for r in catalogue_rows]
    eligible = [r for r in rows if _eligible(r)]
    units, audit = _deduplicate(eligible)
    selected = set(map(str, selected_source_ids))
    selected_units = [r for r in units if selected.intersection(r['_source_ids'])]
    focal_units, focal_audit = _deduplicate([dict(r) for r in focal_rows])
    public_aliases = set().union(*(set(r['_aliases']) for r in units)) if units else set()
    overlap = []
    for unit in focal_units:
        matches = sorted(set(unit['_aliases']) & public_aliases)
        if matches:
            overlap.append({'focal_sample_unit': unit['_unit'], 'matched_aliases': matches})
    tables = []
    summaries = []
    for depth in depths:
        if depth not in DEPTHS:
            raise ValueError('Supported cgLIN prefix depths are 5, 6 and 7')
        for cohort, cohort_units in (('public_catalogue', units), ('selected_context', selected_units),
                                     ('focal_survey', focal_units)):
            groups = defaultdict(list)
            for unit in cohort_units:
                key = unit.get(f'cglin_group_{depth}') or ''
                status = unit.get(f'cglin_status_{depth}', unit.get('cglin_status', 'missing'))
                valid = bool(key) and status in {'resolved', 'provisional', 'complete', 'partial'}
                groups[(key if valid else '', 'lineage' if valid else 'assignment_coverage',
                        '' if valid else str(status), str(unit.get('cglin_scheme', 'unknown')),
                        str(unit.get('cglin_scheme_version', 'unknown')),
                        str(unit.get('source', 'pathogenwatch')))].append(unit)
            assigned = 0
            for (key, category, missing_status, scheme, version, source), members in sorted(groups.items()):
                if category == 'lineage':
                    assigned += len(members)
                counts = Counter(str(r.get('country') or 'Unknown') for r in members)
                raw_counts = Counter()
                for member in members:
                    raw_counts[str(member.get('country') or 'Unknown')] += member['_raw_n']
                representative = members[0]
                label = key
                if key:
                    try:
                        parts = json.loads(key)
                        label = f'{parts[0]} / {parts[1]}: ' + '.'.join(map(str, parts[2]))
                    except (ValueError, TypeError, IndexError):
                        label = key
                else:
                    label = f'Unresolved: {missing_status} ({scheme} / {version}; {source})'
                for country, count in sorted(counts.items()):
                    tables.append(dict(zip(FIELDS, (
                        source,
                        scope.get('snapshot', representative.get('snapshot', 'unspecified')),
                        scope.get('description', 'Frozen public same-ST catalogue'),
                        scope.get('filters', 'Explicit QC pass; deduplicated; includes undated records'),
                        scheme, version, depth, key, label,
                        category, '|'.join(sorted({str(r.get(f'cglin_status_{depth}',
                            r.get('cglin_status', 'missing'))) for r in members})),
                        cohort, country, count, len(members), len(members) - counts['Unknown'],
                        counts['Unknown'], count / len(members) * 100,
                        'deduplicated biological sample where accessions support identity; otherwise record unit',
                        raw_counts[country]))))
            summaries.append({'depth': depth, 'cohort': cohort, 'total_units': len(cohort_units),
                              'assigned_units': assigned, 'unresolved_units': len(cohort_units) - assigned,
                              'lineage_groups': sum(k[1] == 'lineage' for k in groups),
                              'singleton_groups': sum(k[1] == 'lineage' and len(v) == 1
                                                      for k, v in groups.items())})
    return {'rows': tables, 'summaries': summaries, 'duplicate_audit': audit,
            'focal_duplicate_audit': focal_audit, 'focal_overlap': overlap,
            'raw_catalogue_records': sum(int(r.get('raw_genome_count', 1)) for r in rows),
            'qc_eligible_records': sum(int(r.get('raw_genome_count', 1)) for r in eligible),
            'sample_units': len(units), 'qc_excluded_records': sum(int(r.get('raw_genome_count', 1)) for r in rows if not _eligible(r)),
            'identity_unresolved_units': sum(not (r.get('identity_resolved') or
                any(a.startswith('biosample:') for a in r['_aliases'])) for r in units), 'scope': scope}


def _write_table(path: Path, rows: list[dict], fields: Iterable[str], delimiter=','):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def _plot(rows: list[dict], path: Path, title: str, percentage: bool, caption: str):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    countries = sorted({r['country'] for r in rows})
    if len(countries) > 13:
        raise ValueError('Plot country grouping was not applied')
    groups = list(dict.fromkeys((r['prefix_key'], r['group_label'], r['denominator'],
                                r['known_country_n'], r['unknown_n']) for r in rows))
    fig, ax = plt.subplots(figsize=(14, max(5, len(groups) * .48 + 3)))
    lookup = {(r['prefix_key'], r['group_label'], r['country']): r for r in rows}
    left = [0.] * len(groups)
    for country in countries:
        values = [lookup.get((key, label, country), {}).get('percentage' if percentage else 'count', 0)
                  for key, label, *_ in groups]
        ax.barh(range(len(groups)), values, left=left, label=country, color=country_colour(country),
                edgecolor='white', linewidth=.4)
        left = [a + b for a, b in zip(left, values)]
    ax.set_yticks(range(len(groups)),
                  [f'{label}  N={n}; known={known}; Unknown={unknown}'
                   + (' (singleton)' if n == 1 else ' (small)' if n < 5 else '')
                   for _, label, n, known, unknown in groups], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Within-group percentage: n/N, including Unknown (%)' if percentage
                  else 'Deduplicated sample units (n); Unknown included in N')
    if percentage:
        ax.set_xlim(0, 100)
    ax.set_title(title, fontsize=13, pad=15)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.16), ncol=min(5, len(countries)),
              frameon=False, fontsize=9)
    fig.text(.02, .015, caption, fontsize=9, va='bottom', wrap=True)
    fig.subplots_adjust(left=.35, right=.98, bottom=.28, top=.88)
    fig.savefig(path.with_suffix('.svg'), bbox_inches='tight')
    fig.savefig(path.with_suffix('.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)


def generate_context_geography(catalogue_rows: Iterable[Mapping[str, Any]], output_dir: str | Path,
                               *, selected_source_ids: Iterable[str] = (),
                               focal_rows: Iterable[Mapping[str, Any]] = (),
                               depths: Iterable[int] = DEPTHS,
                               scope: Mapping[str, Any] | None = None) -> dict:
    """Write CSV/TSV, JSON audit, figures, and self-contained embeddable HTML fragment."""
    depths = tuple(depths)
    result = geography_tables(catalogue_rows, selected_source_ids=selected_source_ids,
                              focal_rows=focal_rows, depths=depths, scope=scope)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_table(output_dir / 'country_composition.csv', result['rows'], FIELDS)
    _write_table(output_dir / 'country_composition.tsv', result['rows'], FIELDS, '\t')
    matrix = {}
    countries = sorted({r['country'] for r in result['rows']})
    matrix_fields = [f for f in FIELDS if f not in {'country', 'count', 'percentage', 'raw_record_count'}]
    for item in result['rows']:
        key = tuple(item[f] for f in matrix_fields)
        if key not in matrix:
            matrix[key] = {f: item[f] for f in matrix_fields}
            matrix[key].update({country: 0 for country in countries})
        matrix[key][item['country']] = item['count']
    _write_table(output_dir / 'country_composition_matrix.csv', list(matrix.values()), matrix_fields + countries)
    (output_dir / 'country_composition_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    outputs = [str(output_dir / name) for name in ('country_composition.csv',
               'country_composition.tsv', 'country_composition_matrix.csv', 'country_composition_audit.json')]
    fragments = ['<section id="context-geography"><h2>Country composition within cgLIN groups</h2>',
                 '<p>' + html.escape(CAUTION) + '</p>',
                 '<p>Public denominators include QC-passing undated records before selection. '
                 'Focal denominators are separate; public overlap is reported in the audit.</p>']
    if not result['rows']:
        fragments.append('<p>No eligible annotated context or focal records are available.</p>')
    fragments.append(f"<p>Focal sample units overlapping public catalogue: {len(result['focal_overlap'])}; "
                     "their survey counts remain separate from public proportions.</p>")
    for name in ('country_composition.csv', 'country_composition.tsv', 'country_composition_matrix.csv', 'country_composition_audit.json'):
        encoded = base64.b64encode((output_dir / name).read_bytes()).decode()
        fragments.append(f'<p><a download="{name}" href="data:application/octet-stream;base64,{encoded}">Download {name}</a></p>')
    country_totals = Counter()
    for item in result['rows']:
        if item['cohort'] == 'public_catalogue' and item['prefix_depth'] == min(depths):
            country_totals[item['country']] += item['count']
    if not country_totals:
        for item in result['rows']:
            if item['prefix_depth'] == min(depths):
                country_totals[item['country']] += item['count']
    visible_countries = {c for c, _ in sorted(country_totals.items(), key=lambda pair: (-pair[1], pair[0]))
                         if c != 'Unknown'}
    visible_countries = set(sorted(visible_countries, key=lambda c: (-country_totals[c], c))[:11])
    visible_countries.add('Unknown')
    for summary in result['summaries']:
        cohort, depth = summary['cohort'], summary['depth']
        fragments.append('<h3>' + html.escape(f'{cohort}: full-prefix depth {depth}') + '</h3>')
        fragments.append(f"<p>Total N={summary['total_units']}; resolved N={summary['assigned_units']}; "
                         f"unresolved N={summary['unresolved_units']}. Prefix depths are not SNP cutoffs.</p>")
        for category in ('lineage', 'assignment_coverage'):
            panel = [r for r in result['rows'] if r['cohort'] == cohort
                     and r['prefix_depth'] == depth and r['assignment_category'] == category]
            if not panel:
                if category == 'lineage':
                    fragments.append('<p>No resolved comparable cgLIN groups at this depth.</p>')
                continue
            # Fixed page size retains singletons without making labels unreadable.
            labels = list(dict.fromkeys((r['prefix_key'], r['group_label']) for r in panel))
            for page in range(0, len(labels), 24):
                page_labels = set(labels[page:page + 24])
                page_rows = [r for r in panel if (r['prefix_key'], r['group_label']) in page_labels]
                for mode in ('counts', 'percentages'):
                    stem = f'{cohort}_depth_{depth}_{category}_{mode}_{page // 24 + 1}'
                    visual = {}
                    for item in page_rows:
                        item = dict(item)
                        if item['country'] not in visible_countries:
                            item['country'] = 'Other'
                        merge_key = (item['prefix_key'], item['group_label'], item['country'])
                        if merge_key in visual:
                            for field in ('count', 'percentage', 'raw_record_count'):
                                visual[merge_key][field] += item[field]
                        else:
                            visual[merge_key] = item
                    path = output_dir / stem
                    title = f'{cohort} — depth {depth} — {category.replace("_", " ")} — {mode}'
                    caption = ('Scope: ' + str(result['scope'].get('description', 'Frozen public same-ST catalogue'))
                               + '\nSnapshot: ' + str(result['scope'].get('snapshot', 'unspecified'))
                               + '. Source: ' + ', '.join(sorted({str(r['source']) for r in page_rows}))
                               + '\nFilters: ' + str(result['scope'].get('filters', 'Explicit QC pass; deduplicated; undated records included'))
                               + '\nDisplay: top 11 named countries by public-catalogue count + Unknown; remaining countries = Other. Full named-country table downloadable.\n' + CAUTION)
                    _plot(list(visual.values()), path, title, mode == 'percentages', caption)
                    outputs.extend([str(path.with_suffix('.svg')), str(path.with_suffix('.png'))])
                    svg = path.with_suffix('.svg').read_text()
                    fragments.append('<svg' + svg.split('<svg', 1)[1])
    fragments.append('</section>')
    result['outputs'] = outputs
    result['report_html'] = '\n'.join(fragments)
    for name in ('country_composition.html', 'fragment.html'):
        (output_dir / name).write_text(result['report_html'])
        result['outputs'].append(str(output_dir / name))
    document = ('<!doctype html><html lang="en"><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width,initial-scale=1">'
                '<title>cgLIN country composition</title><style>'
                'body{font:16px system-ui,sans-serif;margin:2rem;max-width:1500px}'
                'svg{width:100%;height:auto;display:block;margin:1rem 0}'
                '</style><body>' + result['report_html'] + '</body></html>')
    (output_dir / 'index.html').write_text(document)
    result['outputs'].append(str(output_dir / 'index.html'))
    return result
