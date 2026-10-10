from io import StringIO

from Bio import Phylo

from chronoclade.profile_network import build_profile_network
from chronoclade.location_network.static import draw_location_network as draw_profile_network


def tree():
    return Phylo.read(StringIO('((q:1,c1:1):1,(c2:1,c3:4):1);'), 'newick')


def records():
    return [dict(sample_id=name, origin='query' if name == 'q' else 'context',
                 country=country, cglin_raw=code, cglin_scheme='scgMLST629_S',
                 cglin_scheme_version='v1')
            for name, country, code in [('q', 'Greece', '0,0,0,0,1,0,0,0,0,0'),
                                        ('c1', 'Italy', '0,0,0,0,1,0,0,0,0,1'),
                                        ('c2', 'Italy', '0,0,0,0,1,0,0,0,0,2'),
                                        ('c3', 'France', '0,0,0,0,2,0,0,0,0,0')]]


def test_undirected_pairs_count_once_per_root():
    data = build_profile_network(tree(), records())
    pairs = [(e['source'], e['target']) for e in data['edges']]
    assert len(pairs) == len(set(pairs))
    assert all(a < b for a, b in pairs)
    assert all(0 < e['root_fraction'] <= 1 for e in data['edges'])
    assert all(e['roots_tested'] == len(data['tested_roots']) for e in data['edges'])
    assert data['tested_roots'][0] == 'q'
    assert data['tested_roots'][1] == 'c3'  # Farthest, not first alphabetical ID.
    assert 'not probability or confidence' in data['interpretation']


def test_unknown_country_not_inferred_node_and_known_isolate_preserved(tmp_path):
    rows = records()
    for row in rows:
        row['country'] = 'Unknown'
    rows[0]['country'] = 'Greece'
    data = build_profile_network(tree(), rows)
    assert data['unknown_country_count'] == 3
    assert data['edges'] == []
    assert data['nodes'] == [dict(country='Greece', count=1, input_count=1,
                                   context_count=0, is_input_country=True)]
    path = tmp_path / 'network.svg'
    draw_profile_network(data, path)
    assert 'Greece' in path.read_text()


def test_subgroup_version_compatibility_and_nearest_ties():
    rows = records()
    rows[2]['cglin_scheme_version'] = 'v2'
    nearest = [dict(query_id='q', context_id=cid, status='matched',
                    allele_differences=1, shared_called_loci=629) for cid in ('c1', 'c2')]
    nearest += [nearest[0].copy()]  # Accidental duplicate must not inflate counts.
    data = build_profile_network(tree(), rows, nearest)
    assert data['nearest_edges'][0]['context_count'] == 2
    assert data['nearest_edges'][0]['query_count'] == 1
    view = data['views'][0]
    assert view['sample_count'] == 2
    assert view['nearest_edges'][0]['context_ids'] == ['c1', 'c2']
    assert 'complete selected CG pool' in view['nearest_search_scope']
    assert view['nearest_edges'][0]['shared_loci_min'] == 629
    assert view['scope']['prefix'] == '0,0,0,0,1'


def test_same_country_nearest_is_not_inferred_edge():
    rows = records()
    rows[1]['country'] = 'Greece'
    data = build_profile_network(tree(), rows, [dict(query_id='q', context_id='c1',
            status='matched', allele_differences=0, shared_called_loci=620)])
    assert data['nearest_edges'][0]['source'] == data['nearest_edges'][0]['target'] == 'Greece'
    assert all(e['source'] != e['target'] for e in data['edges'])


def test_all_missing_countries_has_no_root_denominator():
    rows = records()
    for row in rows:
        row['country'] = None
    data = build_profile_network(tree(), rows)
    assert data['root_count'] == 0
    assert data['nodes'] == data['edges'] == []


def test_subgroup_keeps_nearest_ties_outside_subgroup_tree():
    nearest = [dict(query_id='q', context_id=cid, status='matched',
                    allele_differences=2, shared_called_loci=629) for cid in ('c1', 'c3')]
    data = build_profile_network(tree(), records(), nearest)
    view = data['views'][0]
    assert view['sample_count'] == 3
    assert sum(row['comparison_count'] for row in view['nearest_edges']) == 2
    assert {row['target'] for row in view['nearest_edges']} == {'Italy', 'France'}


def test_unknown_nearest_countries_retained_as_observation_not_state():
    rows = records()
    rows[0]['country'] = None
    rows[1]['country'] = None
    data = build_profile_network(tree(), rows, [dict(query_id='q', context_id='c1',
            status='matched', allele_differences=0, shared_called_loci=629)])
    assert data['nearest_edges'][0]['source'] == 'Unknown'
    assert data['nearest_edges'][0]['target'] == 'Unknown'
    assert all(n['country'] != 'Unknown' for n in data['nodes'])


def test_nearest_country_outside_comparability_cohort_kept_from_evidence():
    rows = records()
    data = build_profile_network(tree(), rows, [dict(query_id='q', context_id='outside',
            country='Germany', status='matched', allele_differences=1, shared_called_loci=629)])
    assert data['nearest_edges'][0]['target'] == 'Germany'
    assert data['views'][0]['nearest_edges'][0]['context_ids'] == ['outside']


def _brute_force(t, countries, directed=False):
    """Independent enumeration oracle for tiny trees, including wildcard tips."""
    from collections import Counter
    from itertools import product

    nodes = list(t.find_clades(order='preorder'))
    states = sorted(set(countries.values()) - {None})
    free = [n for n in nodes if not n.is_terminal() or countries.get(n.name) is None]
    optimum, histories, assignments = float('inf'), [], []
    for values in product(states, repeat=len(free)):
        assigned = dict(zip(free, values))
        assigned.update({n: countries[n.name] for n in nodes
                         if n.is_terminal() and countries.get(n.name) is not None})
        counts = Counter((assigned[n], assigned[c]) if directed else
                         tuple(sorted((assigned[n], assigned[c])))
                         for n in nodes for c in n.clades if assigned[n] != assigned[c])
        score = sum(counts.values())
        if score < optimum:
            optimum, histories, assignments = score, [], []
        if score == optimum:
            histories.append(counts)
            assignments.append(assigned)
    pairs = set().union(*(set(h) for h in histories))
    return optimum, {pair: (min(h[pair] for h in histories), max(h[pair] for h in histories))
                     for pair in pairs}, {i: {a[n] for a in assignments}
                                         for i, n in enumerate(nodes)}


def test_exact_ranges_match_all_optimal_tiny_histories():
    from chronoclade.location_network.reconstruction import _reconstruct

    cases = [('((a,b),(c,d));', dict(a='A', b='B', c='C', d='A')),
             ('(a,b,c,d);', dict(a='A', b='B', c='C', d=None)),
             ('((a,b),c);', dict(a='A', b='A', c='B')),
             ('((a,b),c);', dict(a='A', b=None, c='A'))]
    for newick, countries in cases:
        t = Phylo.read(StringIO(newick), 'newick')
        locations = sorted(set(countries.values()) - {None})
        counts, ranges, audit = _reconstruct(t, countries, locations)
        optimum, expected, allowed = _brute_force(t, countries)
        assert ranges == expected
        assert sum(counts.values()) == optimum == audit['optimum_changes']
        assert sum(b['changed'] for b in audit['branches']) == optimum
        assert all(set(n['allowed_states']) == allowed[n['id']] for n in audit['nodes'])
        by_id = {n['id']: n for n in audit['nodes']}
        for branch in audit['branches']:
            assert branch['source'] == by_id[branch['parent_id']]['inferred_state']
            assert branch['target'] == by_id[branch['child_id']]['inferred_state']
            if branch['unknown_tip']:
                assert not branch['changed']
                assert by_id[branch['child_id']]['state'] is None
        assert (counts, ranges, audit) == _reconstruct(t, countries, locations)


def test_tied_histories_keep_possible_only_pairs_and_nonjoint_ranges():
    from chronoclade.location_network.reconstruction import _reconstruct

    t = Phylo.read(StringIO('(a,b,c);'), 'newick')
    counts, ranges, audit = _reconstruct(t, dict(a='A', b='B', c='C'), ['A', 'B', 'C'])
    assert ranges == {('A', 'B'): (0, 1), ('A', 'C'): (0, 1), ('B', 'C'): (0, 1)}
    assert sum(counts.values()) == 2
    assert len(counts) == 2  # The union of three pairs is not one coherent history.
    assert set(audit['nodes'][0]['allowed_states']) == {'A', 'B', 'C'}


def test_weighted_network_counts_are_coherent_and_palette_stable(tmp_path):
    import json
    from chronoclade.location_network.colours import country_palette
    from chronoclade.location_network.static import _weighted_layout

    data = build_profile_network(tree(), records())
    assert sum(e['representative_count'] for e in data['edges']) == data['optimum_changes']
    assert all(e['min_changes'] <= e['representative_count'] <= e['max_changes']
               for e in data['edges'])
    assert all(view['country_colors'][name] == colour for view in data['views']
               for name, colour in data['country_colors'].items()
               if name in view['country_colors'])
    assert country_palette(['Greece'])['Greece'] == data['country_colors']['Greece']
    json.dumps(data)  # No NumPy scalar values or tree objects escape the audit.
    positions = _weighted_layout(data['nodes'], data['edges'])
    assert positions == _weighted_layout(data['nodes'], data['edges'])
    changed_inputs = [dict(n, is_input_country=not n['is_input_country']) for n in data['nodes']]
    assert positions == _weighted_layout(changed_inputs, data['edges'])
    draw_profile_network(data, tmp_path / 'full.svg')
    assert all(n['country'] in (tmp_path / 'full.svg').read_text() for n in data['nodes'])
    draw_profile_network(data, tmp_path / 'possible.svg', include_possible=True)


def test_repeated_changes_produce_weight_two_and_single_tip_zero():
    from chronoclade.location_network.reconstruction import _reconstruct

    t = Phylo.read(StringIO('((a,b),(c,d));'), 'newick')
    counts, ranges, audit = _reconstruct(t, dict(a='A', b='B', c='A', d='B'), ['A', 'B'])
    assert counts == {('A', 'B'): 2}
    assert ranges == {('A', 'B'): (2, 2)}
    assert audit['optimum_changes'] == 2
    one = Phylo.read(StringIO('a;'), 'newick')
    counts, ranges, audit = _reconstruct(one, {'a': 'A'}, ['A'])
    assert counts == ranges == {}
    assert audit['optimum_changes'] == 0


def test_ordered_pair_ranges_match_bruteforce_and_same_coherent_audit():
    from collections import Counter
    from chronoclade.location_network.reconstruction import _reconstruct

    for newick, countries in [
        ('((a,b),(c,d));', dict(a='A', b='B', c='C', d='A')),
        ('(a,b,c,d);', dict(a='A', b='B', c='C', d=None)),
        ('((a,b),(c,d));', dict(a='A', b='B', c='A', d='B')),
        ('((a,b),c);', dict(a='A', b=None, c='A')),
    ]:
        t = Phylo.read(StringIO(newick), 'newick')
        locations = sorted(set(countries.values()) - {None})
        counts, ranges, audit = _reconstruct(t, countries, locations, directed=True)
        optimum, expected, allowed = _brute_force(t, countries, directed=True)
        assert ranges == expected
        assert sum(counts.values()) == optimum == audit['optimum_changes']
        undirected_counts, _, undirected_audit = _reconstruct(t, countries, locations)
        assert audit == undirected_audit
        collapsed = Counter()
        for (a, b), count in counts.items():
            collapsed[tuple(sorted((a, b)))] += count
        assert collapsed == undirected_counts
        assert counts == Counter((b['source'], b['target']) for b in audit['branches'] if b['changed'])
        assert all(set(n['allowed_states']) == allowed[n['id']] for n in audit['nodes'])
        assert all('Unknown' not in pair and None not in pair for pair in ranges)
        assert all(not b['changed'] for b in audit['branches'] if b['unknown_tip'])
        assert (counts, ranges, audit) == _reconstruct(t, countries, locations, directed=True)


def test_directed_count_ranges_depend_on_displayed_root():
    import copy
    from chronoclade.location_network.reconstruction import _reconstruct

    t = Phylo.read(StringIO('((a,b),(c,d));'), 'newick')
    countries = dict(a='A', b='A', c='B', d='B')
    _, original, _ = _reconstruct(t, countries, ['A', 'B'], directed=True)
    assert original == {('A', 'B'): (0, 1), ('B', 'A'): (0, 1)}
    a_root, b_root = copy.deepcopy(t), copy.deepcopy(t)
    a_root.root_with_outgroup('a')
    b_root.root_with_outgroup('c')
    a_counts, a_ranges, a_audit = _reconstruct(a_root, countries, ['A', 'B'], directed=True)
    b_counts, b_ranges, b_audit = _reconstruct(b_root, countries, ['A', 'B'], directed=True)
    assert a_counts == {('A', 'B'): 1} and a_ranges == {('A', 'B'): (1, 1)}
    assert b_counts == {('B', 'A'): 1} and b_ranges == {('B', 'A'): (1, 1)}
    assert a_audit['optimum_changes'] == b_audit['optimum_changes'] == 1


def test_directed_network_metrics_and_renderer(tmp_path):
    from chronoclade.location_network.reconstruction import _network_metrics

    data = build_profile_network(tree(), records())
    assert sum(e['representative_count'] for e in data['directed_edges']) == data['optimum_changes']
    assert sum(m['in_changes'] for m in data['network_metrics']) == data['optimum_changes']
    assert sum(m['out_changes'] for m in data['network_metrics']) == data['optimum_changes']
    metrics = _network_metrics(['A', 'B', 'C', 'isolated'], [
        dict(source='A', target='B', representative_count=2),
        dict(source='B', target='A', representative_count=1),
        dict(source='A', target='C', representative_count=0),
    ])
    assert metrics[0] == dict(country='A', in_degree=1, out_degree=1,
                              degree=2, betweenness=0.0, closeness=1.0,
                              in_changes=1, out_changes=2, source_hub_ratio=0.5)
    assert metrics[-1]['source_hub_ratio'] is None
    assert metrics[-1]['closeness'] is None
    draw_profile_network(data, tmp_path / 'directed.svg')
    assert 'Directed country-state changes' in (tmp_path / 'directed.svg').read_text()
    assert 'root-dependent' in (tmp_path / 'directed.svg').read_text()
    draw_profile_network(data, tmp_path / 'directed_possible.svg', include_possible=True)


def test_network_metrics_match_directed_unweighted_strainhub_definitions():
    from chronoclade.location_network.reconstruction import _network_metrics

    metrics = _network_metrics(['A', 'B', 'C', 'D', 'isolated'], [
        dict(source='C', target='A', representative_count=1),
        dict(source='A', target='B', representative_count=5),
        dict(source='B', target='C', representative_count=2),
        dict(source='B', target='D', representative_count=3),
    ])
    by_country = {metric['country']: metric for metric in metrics}

    # Degree and SHR count links; in/out changes retain representative weights.
    assert (by_country['A']['in_degree'], by_country['A']['out_degree'],
            by_country['A']['degree']) == (1, 1, 2)
    assert by_country['A']['source_hub_ratio'] == 0.5
    assert (by_country['A']['in_changes'], by_country['A']['out_changes']) == (1, 5)
    assert by_country['A']['source_hub_ratio'] != 5 / 6

    # Directed shortest paths pass through A and B; mode="all" closeness uses
    # undirected distances and remains unnormalised.
    assert [by_country[c]['betweenness'] for c in ['A', 'B', 'C', 'D']] == [2.0, 3.0, 1.0, 0.0]
    assert [by_country[c]['closeness'] for c in ['A', 'B', 'C', 'D']] == [0.25, 1 / 3, 0.25, 0.2]
    assert by_country['isolated']['betweenness'] == 0.0
    assert by_country['isolated']['closeness'] is None
    assert by_country['isolated']['source_hub_ratio'] is None
