"""Shared foundations must remain usable independently of analysis stages."""

import subprocess
import sys

import pytest

from chronoclade.errors import WorkflowError
from chronoclade.selections import selected_context_records


def test_shared_style_and_selection_imports_do_not_load_scientific_stages():
    subprocess.run(
        [sys.executable, "-c", """
import sys
from chronoclade.report_components.styles import report_styles
from chronoclade.selections import selected_context_records
from chronoclade.datasets import PreparedDataset, load_dataset
assert "@font-face" in report_styles()
for name in (
    "chronoclade.report", "chronoclade.temporal_report", "chronoclade.profile_analysis",
    "chronoclade.staged_workflow", "chronoclade.lineage", "chronoclade.pathogenwatch",
    "chronoclade.context_refinement", "torch", "esm",
):
    assert name not in sys.modules, name
"""],
        check=True,
    )


def test_saved_selection_preserves_identity_and_does_not_mutate_pool():
    queries = [{"sample_id": "query"}]
    context = [{"sample_id": "a", "accession": "SAMN1"}, {"sample_id": "b"}]
    selection = dict(
        available_context_ids=["a", "b"], query_ids=["query"],
        selected_context_ids=["a"], selected_sample_ids=["query", "a"],
        decisions=[dict(sample_id="a", reason="user_requested")],
    )
    assert selected_context_records(selection, queries, context) == [
        dict(sample_id="a", accession="SAMN1", selection_reason="user_requested")
    ]
    assert "selection_reason" not in context[0]
    selection["decisions"][0].pop("reason")
    with pytest.raises(WorkflowError, match="decisions"):
        selected_context_records(selection, queries, context)
