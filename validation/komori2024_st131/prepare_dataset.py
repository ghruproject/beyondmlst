#!/usr/bin/env python3
"""Prepare the accession-defined Komori et al. ST131-H30 validation set.

The authors' Figure S5 workbook is the source of truth. In particular, this
script never searches NCBI by strain name: strain names are not identifiers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import http.cookiejar
import re
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from xml.etree import ElementTree


ARTICLE_DOI = "10.1128/aac.00817-24"
SUPPLEMENT_URL = (
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/11373201/bin/aac.00817-24-s0005.xlsx"
)
SUPPLEMENT_SHA256 = "20735bdb7b03db302c3642d6a12cb31408a9bd60bcb3ff05c75e68d99241ebe7"
EXPECTED_ALL_STRAINS = 145
EXPECTED_C0_C1_STRAINS = 96
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
USER_AGENT = "ChronoClade-validation/0.1 (https://github.com/ghruproject/chronoclade)"
XML_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _request_bytes(
    url: str,
    *,
    data: bytes | None = None,
    headers: dict[str, str] | None = None,
    timeout: int = 240,
) -> bytes:
    request_headers = {"User-Agent": USER_AGENT, **(headers or {})}
    for attempt in range(8):
        try:
            request = urllib.request.Request(url, data=data, headers=request_headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (OSError, TimeoutError):
            if attempt == 7:
                raise
            time.sleep(min(60, 2 ** (attempt + 1)))
    raise RuntimeError("unreachable")


def download_supplement(destination: Path) -> None:
    """Download the workbook through PMC's lightweight proof-of-work page."""
    cookie_jar = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cookie_jar))
    request = urllib.request.Request(SUPPLEMENT_URL, headers={"User-Agent": USER_AGENT})
    with opener.open(request, timeout=60) as response:
        first_response = response.read()

    if first_response.startswith(b"PK"):
        workbook = first_response
    else:
        page = first_response.decode("utf-8", errors="replace")
        challenge_match = re.search(r'const POW_CHALLENGE = "([^"]+)"', page)
        difficulty_match = re.search(r'const POW_DIFFICULTY = "(\d+)"', page)
        cookie_match = re.search(r'const POW_COOKIE_NAME = "([^"]+)"', page)
        if not (challenge_match and difficulty_match and cookie_match):
            raise RuntimeError("PMC did not return the workbook or its download challenge")
        challenge = challenge_match.group(1)
        prefix = "0" * int(difficulty_match.group(1))
        nonce = 0
        while not hashlib.sha256(f"{challenge}{nonce}".encode()).hexdigest().startswith(prefix):
            nonce += 1
        request = urllib.request.Request(
            SUPPLEMENT_URL,
            headers={
                "User-Agent": USER_AGENT,
                "Cookie": f"{cookie_match.group(1)}={challenge},{nonce}",
            },
        )
        with opener.open(request, timeout=60) as response:
            workbook = response.read()

    digest = hashlib.sha256(workbook).hexdigest()
    if digest != SUPPLEMENT_SHA256:
        raise RuntimeError(
            "Figure S5 workbook checksum changed; inspect the new source before continuing "
            f"(expected {SUPPLEMENT_SHA256}, received {digest})"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(workbook)


def _column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if not letters:
        raise ValueError(f"Invalid spreadsheet cell reference: {reference}")
    number = 0
    for character in letters.group():
        number = number * 26 + ord(character) - 64
    return number - 1


def read_figure_s5(workbook: Path) -> list[dict[str, str]]:
    """Read Figure S5 without adding an Excel library to the runtime."""
    with zipfile.ZipFile(workbook) as archive:
        shared = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = [
            "".join(node.text or "" for node in item.findall(".//m:t", XML_NS))
            for item in shared.findall("m:si", XML_NS)
        ]
        sheet = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    raw_rows: list[dict[int, str]] = []
    for row in sheet.findall(".//m:sheetData/m:row", XML_NS):
        values: dict[int, str] = {}
        for cell in row.findall("m:c", XML_NS):
            value_node = cell.find("m:v", XML_NS)
            value = "" if value_node is None else value_node.text or ""
            if cell.get("t") == "s" and value:
                value = strings[int(value)]
            values[_column_index(cell.get("r", ""))] = value.strip()
        raw_rows.append(values)

    header_row = next(
        row
        for row in raw_rows
        if "Strain name" in row.values() and "GenBank accession no." in row.values()
    )
    header = {column: value for column, value in header_row.items() if value}
    strain_column = next(column for column, name in header.items() if name == "Strain name")
    rows = [
        {name: raw.get(column, "") for column, name in header.items()}
        for raw in raw_rows[raw_rows.index(header_row) + 1 :]
        if raw.get(strain_column, "")
    ]
    if len(rows) != EXPECTED_ALL_STRAINS:
        raise RuntimeError(f"Expected {EXPECTED_ALL_STRAINS} Figure S5 strains, found {len(rows)}")
    subset = [row for row in rows if row["Clade"] in {"C0", "C1"}]
    if len(subset) != EXPECTED_C0_C1_STRAINS:
        raise RuntimeError(f"Expected {EXPECTED_C0_C1_STRAINS} C0/C1 strains, found {len(subset)}")
    accessions = [row["GenBank accession no."].split(".", 1)[0] for row in subset]
    if len(accessions) != len(set(accessions)):
        raise RuntimeError("Figure S5 C0/C1 accession identifiers are not unique")
    return subset


def _post_efetch(parameters: dict[str, str], *, timeout: int = 240) -> bytes:
    data = urllib.parse.urlencode({**parameters, "tool": "chronoclade"}).encode()
    return _request_bytes(NCBI_EFETCH, data=data, timeout=timeout)


def fetch_ncbi_records(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Fetch the exact nucleotide records named by the paper."""
    records: dict[str, dict[str, str]] = {}
    accessions = [row["GenBank accession no."].split(".", 1)[0] for row in rows]
    for start in range(0, len(accessions), 20):
        batch = accessions[start : start + 20]
        root = ElementTree.fromstring(
            _post_efetch(
                {
                    "db": "nuccore",
                    "id": ",".join(batch),
                    "rettype": "gb",
                    "retmode": "xml",
                }
            )
        )
        for sequence in root.findall("GBSeq"):
            primary = sequence.findtext("GBSeq_primary-accession", "")
            qualifiers = {
                item.findtext("GBQualifier_name", ""): item.findtext("GBQualifier_value", "")
                for item in sequence.findall(
                    "./GBSeq_feature-table/GBFeature[GBFeature_key='source']/"
                    "GBFeature_quals/GBQualifier"
                )
            }
            xrefs = {
                item.findtext("GBXref_dbname", ""): item.findtext("GBXref_id", "")
                for item in sequence.findall("./GBSeq_xrefs/GBXref")
            }
            records[primary] = {
                "accession_version": sequence.findtext("GBSeq_accession-version", ""),
                "organism": sequence.findtext("GBSeq_organism", ""),
                "definition": sequence.findtext("GBSeq_definition", ""),
                "sequence_length": sequence.findtext("GBSeq_length", ""),
                "strain": qualifiers.get("strain", qualifiers.get("isolate", "")),
                "collection_date": qualifiers.get("collection_date", ""),
                "country": qualifiers.get("geo_loc_name", qualifiers.get("country", "")),
                "biosample": xrefs.get("BioSample", ""),
                "assembly": xrefs.get("Assembly", ""),
            }
        time.sleep(0.4)
    missing = sorted(set(accessions) - set(records))
    if missing:
        raise RuntimeError("NCBI returned no nucleotide record for: " + ", ".join(missing))
    return records


def _normalise(value: str) -> str:
    return "".join(character.lower() for character in value if character.isalnum())


def _same_country(paper: str, ncbi: str) -> bool:
    aliases = {"uk": "unitedkingdom", "usa": "unitedstates"}
    paper_value = aliases.get(_normalise(paper), _normalise(paper))
    ncbi_value = _normalise(ncbi.split(":", 1)[0])
    ncbi_value = aliases.get(ncbi_value, ncbi_value)
    return bool(ncbi_value) and paper_value == ncbi_value


def build_provenance(
    rows: list[dict[str, str]], records: dict[str, dict[str, str]]
) -> list[dict[str, str]]:
    provenance: list[dict[str, str]] = []
    for row in rows:
        accession = row["GenBank accession no."].split(".", 1)[0]
        record = records[accession]
        paper_year = row["Year"]
        ncbi_year_match = re.search(r"(?:19|20)\d{2}", record["collection_date"])
        ncbi_year = ncbi_year_match.group() if ncbi_year_match else ""
        strain_matches = _normalise(row["Strain name"]) == _normalise(record["strain"])
        provenance.append(
            {
                "strain": row["Strain name"],
                "paper_accession": accession,
                "paper_year": paper_year,
                "paper_country": row["Country"],
                "paper_source": row["Source"],
                "paper_clade": row["Clade"],
                "paper_subclade": row["C1-subclade"],
                **{f"ncbi_{key}": value for key, value in record.items()},
                "species_check": "pass"
                if record["organism"].startswith("Escherichia coli")
                else "fail",
                "strain_name_check": "pass" if strain_matches else "review",
                "year_check": "pass" if ncbi_year == paper_year else "review",
                "country_check": "pass"
                if _same_country(row["Country"], record["country"])
                else "review",
                "source_doi": ARTICLE_DOI,
            }
        )
    species_failures = [row["strain"] for row in provenance if row["species_check"] == "fail"]
    if species_failures:
        raise RuntimeError(
            "The paper accession set contains non-E. coli records: " + ", ".join(species_failures)
        )
    return provenance


def _safe_name(value: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in value).strip("_")


def write_tables(rows: list[dict[str, str]], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    with (output / "accession_provenance.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with (output / "metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "sample_id",
            "assembly",
            "collection_date",
            "location",
            "species",
            "lineage",
            "origin",
            "is_reference",
            "patient_id",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            sample_id = _safe_name(row["strain"])
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "assembly": f"assemblies/{sample_id}.fasta",
                    "collection_date": row["paper_year"],
                    "location": row["paper_country"],
                    "species": "Escherichia_coli",
                    "lineage": "ST131-H30-C0-C1",
                    "origin": "retrospective",
                    "is_reference": str(row["strain"] == "CD306").lower(),
                    "patient_id": "",
                }
            )


def _fasta_records(payload: bytes) -> dict[str, bytes]:
    records: dict[str, bytes] = {}
    for section in payload.split(b">")[1:]:
        header, _, sequence = section.partition(b"\n")
        accession = header.split(maxsplit=1)[0].decode().split(".", 1)[0]
        records[accession] = b">" + header + b"\n" + sequence
    return records


def download_chromosomes(rows: list[dict[str, str]], output: Path) -> None:
    destination = output / "assemblies"
    destination.mkdir(parents=True, exist_ok=True)
    pending = [
        row for row in rows if not (destination / f"{_safe_name(row['strain'])}.fasta").is_file()
    ]
    for start in range(0, len(pending), 8):
        batch = pending[start : start + 8]
        payload = _post_efetch(
            {
                "db": "nuccore",
                "id": ",".join(row["paper_accession"] for row in batch),
                "rettype": "fasta",
                "retmode": "text",
            },
            timeout=360,
        )
        fasta = _fasta_records(payload)
        for row in batch:
            accession = row["paper_accession"]
            if accession not in fasta:
                raise RuntimeError(f"NCBI FASTA response omitted {accession}")
            path = destination / f"{_safe_name(row['strain'])}.fasta"
            path.write_bytes(fasta[accession])
        time.sleep(0.4)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "verified")
    parser.add_argument("--supplement", type=Path)
    parser.add_argument("--download", action="store_true", help="Download the 96 exact chromosomes")
    args = parser.parse_args()

    supplement = args.supplement or args.output / "source" / "aac.00817-24-s0005.xlsx"
    if not supplement.is_file():
        download_supplement(supplement)
    digest = hashlib.sha256(supplement.read_bytes()).hexdigest()
    if digest != SUPPLEMENT_SHA256:
        raise RuntimeError(f"Unexpected Figure S5 workbook checksum: {digest}")

    figure_rows = read_figure_s5(supplement)
    records = fetch_ncbi_records(figure_rows)
    provenance = build_provenance(figure_rows, records)
    write_tables(provenance, args.output)
    if args.download:
        download_chromosomes(provenance, args.output)

    reviews = sum(
        any(
            row[field] == "review" for field in ("strain_name_check", "year_check", "country_check")
        )
        for row in provenance
    )
    print(
        f"Prepared {len(provenance)} accession-defined C0/C1 chromosomes; "
        f"{reviews} records have metadata differences to review"
    )


if __name__ == "__main__":
    main()
