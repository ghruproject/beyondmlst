from chronoclade.sample_labels import label_analysis, sample_label, sample_labels


def test_accession_precedence_focal_identity_and_fallback():
    row = {"sample_id": "PW_1181HF1hzDGp5yG8xz4i5h", "origin": "context",
           "aliases": ["SAMN46159676", "SRR32641190", "GCA_123.1"]}
    assert sample_label(row) == "SRR32641190"
    assert sample_label({**row, "aliases": ["SAMN46159676"]}) == "SAMN46159676"
    assert sample_label({**row, "aliases": []}) == row["sample_id"]
    assert sample_label({**row, "sample_id": "Hospital-A", "origin": "local"}) == "Hospital-A"


def test_repeated_accessions_do_not_merge_genomes():
    labels = sample_labels([
        {"sample_id": "PW_a", "run_accessions": ["SRR123"]},
        {"sample_id": "PW_b", "run_accessions": ["SRR123"]},
    ])
    assert labels == {"PW_a": "SRR123 (PW_a)", "PW_b": "SRR123 (PW_b)"}
    analysis = {"nearest_neighbours": [{"query_id": "PW_a", "context_id": "PW_b"}],
                "paths": {"summary": "PW_a"}}
    rendered = label_analysis(analysis, labels)
    assert rendered["nearest_neighbours"][0]["context_id"] == "SRR123 (PW_b)"
    assert rendered["paths"] == analysis["paths"]
    assert analysis["nearest_neighbours"][0]["context_id"] == "PW_b"
