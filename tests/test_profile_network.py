from io import StringIO

from Bio import Phylo

from chronoclade.profile_network import build_profile_network, draw_profile_network


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
