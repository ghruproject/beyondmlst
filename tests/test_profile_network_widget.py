"""Presentation contracts for offline network widgets, independent of inference."""
import copy
import json
import re
from pathlib import Path

from chronoclade.location_network.widget import interactive_network_html, widget_assets


def row():
    return dict(
        nodes=[dict(country=country, count=3, is_input_country=country == 'Greece')
               for country in ('Greece', 'Italy', 'France', 'Germany')],
        country_colors={'Greece': '#123456', 'Italy': '#aabbcc'},
        network_metrics=[dict(country=country, in_degree=i, out_degree=1, degree=i+1,
                              betweenness=0.5, closeness=None, source_hub_ratio=1/(i+1))
                         for i, country in enumerate(('Greece', 'Italy', 'France', 'Germany'))],
        directed_edges=[
            dict(source='Greece', target='Italy', representative_count=2,
                 min_changes=1, max_changes=3, count_ambiguous=True),
            dict(source='Italy', target='France', representative_count=1,
                 min_changes=1, max_changes=1, count_ambiguous=False),
            dict(source='France', target='Greece', representative_count=0,
                 min_changes=0, max_changes=2, count_ambiguous=True),
        ],
        reconstruction={'nodes': ['coherent history remains untouched']},
    )


def payload(markup):
    return json.loads(re.search(r'<script type="application/json" data-network-data>(.*?)</script>',
                               markup, flags=re.S).group(1))


def test_all_countries_and_representative_metrics_are_preserved_without_mutation():
    data = row()
    original = copy.deepcopy(data)
    markup = interactive_network_html(data, 'cohort / LIN')
    result = payload(markup)
    assert data == original
    assert len(result['nodes']) == 4  # Include isolated Germany.
    assert not result['includePossible']
    assert result['nodes'][0]['color'] == '#123456'
    assert result['nodes'][0]['input']
    assert result['nodes'][2]['metrics']['degree'] == 3
    assert result['nodes'][2]['metrics']['closeness'] is None
    assert result['edges'][0]['count'] == 2
    assert result['edges'][0]['minimum'] == 1
    assert result['edges'][0]['maximum'] == 3
    assert result['edges'][0]['ambiguous']
    assert not result['edges'][1]['ambiguous']
    # Alternatives are available for the toggle, but excluded initially by the viewer.
    assert result['edges'][2]['count'] == 0
    assert 'reconstruction' not in result
    assert 'data-network-country' in markup
    assert 'role="status" aria-live="polite"' in markup
    assert 'role="img"' in markup
    assert markup.index('value="in_degree"') < markup.index('value="out_degree"')
    assert all(f'value="{metric}"' in markup for metric in
               ('degree', 'betweenness', 'closeness', 'source_hub_ratio'))


def test_possible_edges_and_optional_input_focus_keep_original_counts():
    markup = interactive_network_html(row(), 'focus', focus_inputs=True, include_possible=True)
    result = payload(markup)
    assert result['includePossible']
    assert {node['id'] for node in result['nodes']} == {'Greece', 'Italy', 'France'}
    assert [(edge['source'], edge['target'], edge['count']) for edge in result['edges']] == [
        ('Greece', 'Italy', 2), ('France', 'Greece', 0)]
    assert 'data-network-possible checked' in markup


def test_hostile_country_labels_and_identifiers_cannot_close_scripts_or_inject_html():
    hostile = '</script><img src=x onerror=alert(1)>&\u2028'
    data = dict(nodes=[dict(country=hostile, count=1)],
                network_metrics=[dict(country=hostile, in_degree=float('nan'))])
    markup = interactive_network_html(data, hostile)
    assert hostile not in markup
    assert '<img' not in markup
    assert markup.count('</script>') == 1
    assert payload(markup)['nodes'][0]['country'] == hostile
    assert payload(markup)['nodes'][0]['metrics']['in_degree'] is None
    assert '\\u003c/script' in markup


def test_repeated_identifiers_have_unique_label_targets():
    first, second = [interactive_network_html(row(), 'same subgroup') for _ in range(2)]
    assert set(re.findall(r'id="([^"]+)"', first)).isdisjoint(
        re.findall(r'id="([^"]+)"', second))


def test_empty_and_legacy_edge_only_views_are_valid():
    assert payload(interactive_network_html({}, 'empty'))['nodes'] == []
    legacy = dict(edges=[dict(source='A', target='B', representative_count=1)])
    assert {node['id'] for node in payload(interactive_network_html(legacy, 'legacy'))['nodes']} == {'A', 'B'}


def test_offline_assets_expose_repulsion_navigation_and_cleanup_api():
    assets = widget_assets()
    assert 'solver:\'repulsion\'' in assets
    assert 'dragNodes:true,dragView:true,zoomView:true' in assets
    assert 'window.ChronoCladeNetworks={mount,destroy}' in assets
    assert 'instance.network.destroy()' in assets
    assert 'abort.abort()' in assets
    assert 'Math.max(edge.minimum,edge.maximum)>0' in assets
    assert 'element.textContent=text' in assets  # Tooltips are safe DOM text, never HTML metadata.
    assert 'new Option(node.country,node.id)' in assets
    assert 'new vis.Network' in assets
    assert assets.count('data-cc-network-library') == 1
    assert '<script src=' not in assets  # No runtime CDN or other network dependency.
    vendor = Path(__file__).resolve().parents[1] / 'chronoclade' / 'vendor'
    assert 'Permission is hereby granted' in (vendor / 'vis-network.LICENSE-MIT').read_text()
    assert 'Apache License' in (vendor / 'vis-network.LICENSE-APACHE-2.0').read_text()
