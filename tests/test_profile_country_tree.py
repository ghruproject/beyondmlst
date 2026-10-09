from io import StringIO

from Bio import Phylo

from chronoclade.profile_country_tree import draw_country_tree
from chronoclade.profile_network import build_profile_network


def test_country_tree_uses_full_reconstruction_when_display_tips_hidden(tmp_path):
    tree = Phylo.read(StringIO("((q:0.1,c:0.1):0.2,(x:0.1,u:-0.01):0.2);"), "newick")
    rows = [dict(sample_id=ident, country=country, origin=origin, collection_date=date)
            for ident, country, origin, date in [
                ("q", "Greece", "query", "2019-01-01"),
                ("c", "Germany", "context", "2018"),
                ("x", "Italy", "context", "2017"),
                ("u", None, "context", "2351")]]
    network = build_profile_network(tree, rows)
    before = [dict(node) for node in network["reconstruction"]["nodes"]]
    path = tmp_path / "country_tree.svg"
    draw_country_tree(network, rows, path, display_ids=["q", "u"], nearest_ids=["c"])
    svg = path.read_text()
    assert "Showing 2 of 4 profiles" in svg
    assert "q | Greece | 2019-01-01" in svg
    assert "u | Unknown | 2351 (excluded from date analysis)" in svg
    assert "c | Germany" not in svg
    assert network["reconstruction"]["nodes"] == before
    assert network["country_colors"]["Greece"] in svg


def test_country_tree_without_reconstruction_has_no_misleading_image(tmp_path):
    path = tmp_path / "country_tree.svg"
    assert draw_country_tree({}, [], path) is None
    assert not path.exists()


def test_empty_display_selection_does_not_show_the_complete_pool(tmp_path):
    tree = Phylo.read(StringIO("(q:0.1,c:0.1);"), "newick")
    rows = [dict(sample_id="q", origin="query", country="Greece"),
            dict(sample_id="c", origin="context", country="Italy")]
    network = build_profile_network(tree, rows)
    path = tmp_path / "empty.svg"
    assert draw_country_tree(network, rows, path, display_ids=[]) is None
    assert not path.exists()
