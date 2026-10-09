"""Public-only Pathogenwatch catalogues with auditable, offline replayable metadata.

Search/details never send credentials: authenticated visibility is inappropriate for
public context. The API key is available separately for authenticated downloads.
"""

from __future__ import annotations

import calendar
import hashlib
import json
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


class PathogenwatchError(RuntimeError):
    """Upstream failure or unsafe/incomplete catalogue."""


def content_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    ).hexdigest()


# English CLDR territory labels, frozen for deterministic offline country aliases.
_COUNTRIES = {
    "AC": "Ascension Island",
    "AD": "Andorra",
    "AE": "United Arab Emirates",
    "AF": "Afghanistan",
    "AG": "Antigua & Barbuda",
    "AI": "Anguilla",
    "AL": "Albania",
    "AM": "Armenia",
    "AO": "Angola",
    "AQ": "Antarctica",
    "AR": "Argentina",
    "AS": "American Samoa",
    "AT": "Austria",
    "AU": "Australia",
    "AW": "Aruba",
    "AX": "Åland Islands",
    "AZ": "Azerbaijan",
    "BA": "Bosnia & Herzegovina",
    "BB": "Barbados",
    "BD": "Bangladesh",
    "BE": "Belgium",
    "BF": "Burkina Faso",
    "BG": "Bulgaria",
    "BH": "Bahrain",
    "BI": "Burundi",
    "BJ": "Benin",
    "BL": "St. Barthélemy",
    "BM": "Bermuda",
    "BN": "Brunei",
    "BO": "Bolivia",
    "BQ": "Caribbean Netherlands",
    "BR": "Brazil",
    "BS": "Bahamas",
    "BT": "Bhutan",
    "BV": "Bouvet Island",
    "BW": "Botswana",
    "BY": "Belarus",
    "BZ": "Belize",
    "CA": "Canada",
    "CC": "Cocos (Keeling) Islands",
    "CD": "Congo - Kinshasa",
    "CF": "Central African Republic",
    "CG": "Congo - Brazzaville",
    "CH": "Switzerland",
    "CI": "Côte d’Ivoire",
    "CK": "Cook Islands",
    "CL": "Chile",
    "CM": "Cameroon",
    "CN": "China",
    "CO": "Colombia",
    "CP": "Clipperton Island",
    "CQ": "Sark",
    "CR": "Costa Rica",
    "CU": "Cuba",
    "CV": "Cape Verde",
    "CW": "Curaçao",
    "CX": "Christmas Island",
    "CY": "Cyprus",
    "CZ": "Czechia",
    "DE": "Germany",
    "DG": "Diego Garcia",
    "DJ": "Djibouti",
    "DK": "Denmark",
    "DM": "Dominica",
    "DO": "Dominican Republic",
    "DZ": "Algeria",
    "EA": "Ceuta & Melilla",
    "EC": "Ecuador",
    "EE": "Estonia",
    "EG": "Egypt",
    "EH": "Western Sahara",
    "ER": "Eritrea",
    "ES": "Spain",
    "ET": "Ethiopia",
    "EU": "European Union",
    "EZ": "Eurozone",
    "FI": "Finland",
    "FJ": "Fiji",
    "FK": "Falkland Islands",
    "FM": "Micronesia",
    "FO": "Faroe Islands",
    "FR": "France",
    "GA": "Gabon",
    "GB": "United Kingdom",
    "GD": "Grenada",
    "GE": "Georgia",
    "GF": "French Guiana",
    "GG": "Guernsey",
    "GH": "Ghana",
    "GI": "Gibraltar",
    "GL": "Greenland",
    "GM": "Gambia",
    "GN": "Guinea",
    "GP": "Guadeloupe",
    "GQ": "Equatorial Guinea",
    "GR": "Greece",
    "GS": "South Georgia & South Sandwich Islands",
    "GT": "Guatemala",
    "GU": "Guam",
    "GW": "Guinea-Bissau",
    "GY": "Guyana",
    "HK": "Hong Kong SAR China",
    "HM": "Heard & McDonald Islands",
    "HN": "Honduras",
    "HR": "Croatia",
    "HT": "Haiti",
    "HU": "Hungary",
    "IC": "Canary Islands",
    "ID": "Indonesia",
    "IE": "Ireland",
    "IL": "Israel",
    "IM": "Isle of Man",
    "IN": "India",
    "IO": "British Indian Ocean Territory",
    "IQ": "Iraq",
    "IR": "Iran",
    "IS": "Iceland",
    "IT": "Italy",
    "JE": "Jersey",
    "JM": "Jamaica",
    "JO": "Jordan",
    "JP": "Japan",
    "KE": "Kenya",
    "KG": "Kyrgyzstan",
    "KH": "Cambodia",
    "KI": "Kiribati",
    "KM": "Comoros",
    "KN": "St. Kitts & Nevis",
    "KP": "North Korea",
    "KR": "South Korea",
    "KW": "Kuwait",
    "KY": "Cayman Islands",
    "KZ": "Kazakhstan",
    "LA": "Laos",
    "LB": "Lebanon",
    "LC": "St. Lucia",
    "LI": "Liechtenstein",
    "LK": "Sri Lanka",
    "LR": "Liberia",
    "LS": "Lesotho",
    "LT": "Lithuania",
    "LU": "Luxembourg",
    "LV": "Latvia",
    "LY": "Libya",
    "MA": "Morocco",
    "MC": "Monaco",
    "MD": "Moldova",
    "ME": "Montenegro",
    "MF": "St. Martin",
    "MG": "Madagascar",
    "MH": "Marshall Islands",
    "MK": "North Macedonia",
    "ML": "Mali",
    "MM": "Myanmar (Burma)",
    "MN": "Mongolia",
    "MO": "Macao SAR China",
    "MP": "Northern Mariana Islands",
    "MQ": "Martinique",
    "MR": "Mauritania",
    "MS": "Montserrat",
    "MT": "Malta",
    "MU": "Mauritius",
    "MV": "Maldives",
    "MW": "Malawi",
    "MX": "Mexico",
    "MY": "Malaysia",
    "MZ": "Mozambique",
    "NA": "Namibia",
    "NC": "New Caledonia",
    "NE": "Niger",
    "NF": "Norfolk Island",
    "NG": "Nigeria",
    "NI": "Nicaragua",
    "NL": "Netherlands",
    "NO": "Norway",
    "NP": "Nepal",
    "NR": "Nauru",
    "NU": "Niue",
    "NZ": "New Zealand",
    "OM": "Oman",
    "PA": "Panama",
    "PE": "Peru",
    "PF": "French Polynesia",
    "PG": "Papua New Guinea",
    "PH": "Philippines",
    "PK": "Pakistan",
    "PL": "Poland",
    "PM": "St. Pierre & Miquelon",
    "PN": "Pitcairn Islands",
    "PR": "Puerto Rico",
    "PS": "Palestinian Territories",
    "PT": "Portugal",
    "PW": "Palau",
    "PY": "Paraguay",
    "QA": "Qatar",
    "QO": "Outlying Oceania",
    "RE": "Réunion",
    "RO": "Romania",
    "RS": "Serbia",
    "RU": "Russia",
    "RW": "Rwanda",
    "SA": "Saudi Arabia",
    "SB": "Solomon Islands",
    "SC": "Seychelles",
    "SD": "Sudan",
    "SE": "Sweden",
    "SG": "Singapore",
    "SH": "St. Helena",
    "SI": "Slovenia",
    "SJ": "Svalbard & Jan Mayen",
    "SK": "Slovakia",
    "SL": "Sierra Leone",
    "SM": "San Marino",
    "SN": "Senegal",
    "SO": "Somalia",
    "SR": "Suriname",
    "SS": "South Sudan",
    "ST": "São Tomé & Príncipe",
    "SV": "El Salvador",
    "SX": "Sint Maarten",
    "SY": "Syria",
    "SZ": "Eswatini",
    "TA": "Tristan da Cunha",
    "TC": "Turks & Caicos Islands",
    "TD": "Chad",
    "TF": "French Southern Territories",
    "TG": "Togo",
    "TH": "Thailand",
    "TJ": "Tajikistan",
    "TK": "Tokelau",
    "TL": "Timor-Leste",
    "TM": "Turkmenistan",
    "TN": "Tunisia",
    "TO": "Tonga",
    "TR": "Türkiye",
    "TT": "Trinidad & Tobago",
    "TV": "Tuvalu",
    "TW": "Taiwan",
    "TZ": "Tanzania",
    "UA": "Ukraine",
    "UG": "Uganda",
    "UM": "U.S. Outlying Islands",
    "UN": "United Nations",
    "US": "United States",
    "UY": "Uruguay",
    "UZ": "Uzbekistan",
    "VA": "Vatican City",
    "VC": "St. Vincent & Grenadines",
    "VE": "Venezuela",
    "VG": "British Virgin Islands",
    "VI": "U.S. Virgin Islands",
    "VN": "Vietnam",
    "VU": "Vanuatu",
    "WF": "Wallis & Futuna",
    "WS": "Samoa",
    "XA": "Pseudo-Accents",
    "XB": "Pseudo-Bidi",
    "XK": "Kosovo",
    "YE": "Yemen",
    "YT": "Mayotte",
    "ZA": "South Africa",
    "ZM": "Zambia",
    "ZW": "Zimbabwe",
    "ZZ": "Unknown Region",
}
_COUNTRY_ALIASES = {
    "uk": "United Kingdom",
    "u.k.": "United Kingdom",
    "great britain": "United Kingdom",
    "united states of america": "United States",
    "usa": "United States",
    "u.s.a.": "United States",
    "south korea": "South Korea",
    "republic of korea": "South Korea",
    "russian federation": "Russia",
    "viet nam": "Vietnam",
    "czech republic": "Czechia",
    "turkey": "Türkiye",
    "unknown region": "Unknown",
}
_MISSING = {"", "unknown", "none", "null", "na", "n/a", "not provided", "missing", "not collected"}
_ACCESSION = re.compile(
    r"\b(?:SAM[NED]\d+|[SED]RR\d+|GC[AF]_\d+(?:\.\d+)?|PRJ[END][AB]\d+|[SED]RP\d+)\b", re.I
)


def normalize_country(value: Any) -> str:
    if isinstance(value, dict):
        value = value.get("country") or value.get("name") or value.get("code") or ""
    text = str(value or "").strip()
    if text.casefold() in _MISSING:
        return "Unknown"
    # ENA country fields may append a place after a colon.
    text = text.split(":", 1)[0].strip()
    decorated_code = re.match(r"^([A-Z]{2})\s+-\s+", text)
    if decorated_code and decorated_code[1] in _COUNTRIES:
        return _COUNTRIES[decorated_code[1]]
    canonical_names = {name.casefold(): name for name in _COUNTRIES.values()}
    return _COUNTRIES.get(
        text.upper(),
        _COUNTRY_ALIASES.get(text.casefold(), canonical_names.get(text.casefold(), text)),
    )


def _metadata(detail: dict) -> dict[str, Any]:
    metadata = detail.get("metadata") or {}
    if isinstance(metadata, list):
        if any(not isinstance(item, (list, tuple)) or len(item) != 2 for item in metadata):
            raise PathogenwatchError("Unsupported genome metadata pairs")
        return {str(k): v for k, v in metadata if isinstance(k, str)}
    if not isinstance(metadata, dict):
        raise PathogenwatchError("Unsupported genome metadata shape")
    return metadata


def _date_bounds(value: Any) -> tuple[str, str, str]:
    text = str(value or "").strip()
    if not text or text.casefold() in _MISSING:
        return "", "", "missing"
    try:
        if re.fullmatch(r"\d{4}", text):
            date(int(text), 1, 1)
            return f"{text}-01-01", f"{text}-12-31", "year"
        if re.fullmatch(r"\d{4}-\d{2}", text):
            year, month = map(int, text.split("-"))
            date(year, month, 1)
            return f"{text}-01", f"{text}-{calendar.monthrange(year, month)[1]:02}", "month"
        day = date.fromisoformat(text[:10])
        return day.isoformat(), day.isoformat(), "day"
    except ValueError:
        return "", "", "invalid"


def normalize_dates(detail: dict, metadata: dict, search: dict) -> dict:
    raw = next(
        (
            metadata[k]
            for k in ("Collection date", "collection_date", "Date", "date")
            if metadata.get(k)
        ),
        "",
    )
    if isinstance(raw, (list, tuple)) and len(raw) == 2:
        parts = list(raw)
    elif isinstance(raw, str) and "/" in raw:
        parts = raw.split("/", 1)
    else:
        parts = []
    if parts:
        start, _, first = _date_bounds(parts[0])
        _, end, last = _date_bounds(parts[1])
        precision = "interval" if start and end and start <= end else "invalid"
    elif raw:
        start, end, precision = _date_bounds(raw)
    else:
        start, _, first = _date_bounds(detail.get("startDate") or search.get("metadataDate"))
        _, end, last = _date_bounds(
            detail.get("endDate") or detail.get("startDate") or search.get("metadataDate")
        )
        if start and end:
            # Bounds supplied by the service may encode year/month precision.
            year = start[:4]
            month = start[:7]
            if start == f"{year}-01-01" and end == f"{year}-12-31":
                raw, precision = year, "year"
            elif (
                start[:7] == end[:7]
                and start.endswith("-01")
                and end.endswith(f"-{calendar.monthrange(int(year), int(start[5:7]))[1]:02}")
            ):
                raw, precision = month, "month"
            else:
                raw, precision = (
                    start if start == end else f"{start}/{end}",
                    "day" if start == end else "interval",
                )
        else:
            precision = "missing" if not start and not end else "invalid"
    if precision == "invalid":
        start, end = "", ""
    return {
        "collection_date": "/".join(str(p) for p in parts) if parts else str(raw or ""),
        "date_start": start,
        "date_end": end,
        "date_precision": precision,
        "date_raw": raw,
        "dated_cohort_eligible": precision in {"day", "month", "year", "interval"},
    }


def normalize_record(
    search: dict, detail: dict, *, organism_id: str, st: str, mlst_scheme: str = "mlst"
) -> dict:
    source_id = str(search.get("uuid") or "")
    if not source_id or detail.get("uuid") != source_id:
        raise PathogenwatchError("Genome details do not match requested source UUID")
    if search.get("projectAccess") != "PUBLIC" or detail.get("projectAccess", "PUBLIC") != "PUBLIC":
        raise PathogenwatchError("Public catalogue contains non-public genome")
    if str(detail.get("organismId")) != str(organism_id) or str(detail.get(mlst_scheme)) != str(st):
        raise PathogenwatchError("Genome details do not match organism/ST query")
    metadata = _metadata(detail)
    countries = []
    for key in ("Country", "country", "country of origin", "geo_loc_name"):
        if metadata.get(key):
            countries.append({"field": f"metadata.{key}", "value": metadata[key]})
    for field, value in (
        ("details.location", detail.get("location")),
        ("search.location", search.get("location")),
    ):
        if value:
            countries.append({"field": field, "value": value})
    known = [
        (item, normalize_country(item["value"]))
        for item in countries
        if normalize_country(item["value"]) != "Unknown"
    ]
    country = known[0][1] if known else "Unknown"
    accession_text = " ".join(
        str(v)
        for v in [
            detail.get("name", ""),
            detail.get("sampleAccession", ""),
            detail.get("runAccession", ""),
            detail.get("assemblyAccession", ""),
            detail.get("studyAccession", ""),
            *metadata.values(),
        ]
    )
    accessions = sorted(set(_ACCESSION.findall(accession_text.upper())))
    biosamples = [v for v in accessions if v.startswith("SAM")]
    qc = detail.get("qc")
    if qc not in (True, False, None) or isinstance(qc, str):
        raise PathogenwatchError("Unsupported Pathogenwatch QC value")
    row = {
        "source": "pathogenwatch",
        "source_genome_id": source_id,
        "numeric_source_id": detail.get("id"),
        "sample_id": str(detail.get("name") or source_id),
        "species": detail.get("species") or search.get("organism") or "",
        "organism_id": str(organism_id),
        "mlst_scheme": mlst_scheme,
        "mlst_st": str(st),
        "lineage": f"ST{st}",
        "biosample": biosamples[0] if len(biosamples) == 1 else "",
        "biosample_accessions": biosamples,
        "run_accessions": [v for v in accessions if re.match(r"[SED]RR", v)],
        "assembly_accessions": [v for v in accessions if v.startswith(("GCA_", "GCF_"))],
        "study_accessions": [v for v in accessions if v.startswith(("PRJ", "ERP", "SRP", "DRP"))],
        "aliases": [
            source_id,
            *[a for a in accessions if not a.startswith(("PRJ", "ERP", "SRP", "DRP"))],
        ],
        "country": country,
        "country_raw": countries,
        "country_provenance": known[0][0]["field"] if known else "missing",
        "country_conflict": len({v for _, v in known}) > 1,
        "host": str(metadata.get("Host") or metadata.get("host") or ""),
        "isolation_source": str(
            metadata.get("Isolation source") or metadata.get("isolation_source") or ""
        ),
        "qc_pass": qc,
        "source_qc": qc,
        "source_length": detail.get("length"),
        "source_n50": detail.get("n50"),
        "source_contigs": detail.get("contigs"),
        "source_checksum": detail.get("checksum"),
        "source_project_id": search.get("projectId"),
        "public_access": "PUBLIC",
        "raw_details_sha256": content_hash(detail),
    }
    row.update(normalize_dates(detail, metadata, search))
    return row


def _identity_aliases(row: dict) -> set[str]:
    # Study IDs describe a cohort and never identify one biological sample.
    aliases = [
        *row.get("biosample_accessions", []),
        *row.get("run_accessions", []),
        *row.get("assembly_accessions", []),
    ]
    if row.get("biosample"):
        aliases.append(row["biosample"])
    return {re.sub(r"\.[0-9]+$", "", str(v).upper()) for v in aliases if v}


def deduplicate_catalogue(
    rows: list[dict], focal_aliases: Iterable[str] = ()
) -> tuple[list[dict], dict]:
    """Deterministic accession connected-components, preferring QC-pass representatives."""
    ordered = sorted(rows, key=lambda r: r["source_genome_id"])
    parent = list(range(len(ordered)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    seen = {}
    for i, row in enumerate(ordered):
        for alias in _identity_aliases(row) | {f"source:{row['source_genome_id']}"}:
            if alias in seen:
                parent[find(i)] = find(seen[alias])
            seen[alias] = i
    groups: dict[int, list[dict]] = {}
    for i, row in enumerate(ordered):
        groups.setdefault(find(i), []).append(row)
    focal = {re.sub(r"\.[0-9]+$", "", str(v).upper()) for v in focal_aliases if v}
    kept, decisions = [], []
    for members in groups.values():
        members.sort(key=lambda r: (r.get("qc_pass") is not True, r["source_genome_id"]))
        representative = dict(members[0])
        aliases = set().union(*(_identity_aliases(r) for r in members))
        ids = sorted({r["source_genome_id"] for r in members})
        samples = sorted({a for a in aliases if a.startswith("SAM")})
        conflicting = len(samples) > 1
        unit = (
            samples[0]
            if len(samples) == 1
            else min(aliases)
            if aliases
            else f"pathogenwatch:{ids[0]}"
        )
        overlap = sorted((aliases | {v.upper() for v in ids}) & focal)
        conflicts = {
            field: sorted({str(r.get(field)) for r in members})
            for field in ("country", "collection_date", "biosample")
            if len({str(r.get(field)) for r in members}) > 1
        }
        representative.update(
            sample_unit_id=unit,
            source_genome_ids=ids,
            raw_genome_count=len(members),
            identity_resolved=bool(aliases) and not conflicting,
            identity_conflict=conflicting,
        )
        decision = {
            "sample_unit_id": unit,
            "source_genome_ids": ids,
            "representative": representative["source_genome_id"],
            "raw_genome_count": len(members),
            "reason": "focal_overlap"
            if overlap
            else "conflicting_biosamples"
            if conflicting
            else "duplicate"
            if len(members) > 1
            else "retained",
            "focal_matches": overlap,
            "conflicts": conflicts,
        }
        decisions.append(decision)
        if not overlap:
            kept.append(representative)
    kept.sort(key=lambda r: r["source_genome_id"])
    return kept, {
        "raw_records": len(rows),
        "sample_units": len(groups),
        "retained_units": len(kept),
        "focal_excluded_units": len(groups) - len(kept),
        "decisions": decisions,
        "counting_unit": "accession-deduplicated sample units; unresolved source records retained",
    }


def load_api_key(config_path: Path | str | None = None) -> str | None:
    """Read environment or a user-only JSON config without logging its contents."""
    value = os.environ.get("PATHOGENWATCH_API_KEY")
    if value:
        return value.strip()
    path = (
        Path(config_path) if config_path else Path.home() / ".config/chronoclade/pathogenwatch.json"
    )
    if not path.exists():
        return None
    if path.stat().st_mode & 0o077:
        raise PathogenwatchError(
            "Pathogenwatch credential file must have user-only permissions (0600)"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        value = payload.get("api_key")
    except (OSError, ValueError, AttributeError):
        raise PathogenwatchError("Cannot read Pathogenwatch credential configuration") from None
    if not isinstance(value, str) or not value.strip():
        raise PathogenwatchError("Pathogenwatch credential configuration has no api_key")
    return value.strip()


class PathogenwatchClient:
    def __init__(
        self,
        base_url: str = "https://pathogen.watch",
        *,
        api_key: str | None = None,
        timeout: float = 30,
        retries: int = 3,
        workers: int = 6,
        transport: Callable | None = None,
        sleep: Callable = time.sleep,
    ):
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key if api_key is not None else load_api_key()
        self.timeout, self.retries, self.workers = timeout, retries, workers
        self.transport, self.sleep = transport, sleep

    def request_json(
        self,
        method: str,
        path: str,
        *,
        body: dict | None = None,
        params: dict | None = None,
        authenticated: bool = False,
    ):
        url = self.base_url + path + ("?" + urlencode(params) if params else "")
        headers = {"Accept": "application/json", "Content-Type": "application/json"}
        if authenticated and self.api_key:
            headers["X-API-Key"] = self.api_key
        for attempt in range(self.retries + 1):
            try:
                if self.transport:
                    return self.transport(
                        method, path, body=body, params=params or {}, headers=headers
                    )
                request = Request(
                    url,
                    data=json.dumps(body).encode() if body is not None else None,
                    headers=headers,
                    method=method,
                )
                with urlopen(request, timeout=self.timeout) as response:
                    return json.load(response)
            except HTTPError as error:
                if error.code not in (429, 500, 502, 503, 504) or attempt == self.retries:
                    raise PathogenwatchError(f"Pathogenwatch {path} HTTP {error.code}") from None
            except (URLError, TimeoutError, OSError):
                if attempt == self.retries:
                    raise PathogenwatchError(
                        f"Pathogenwatch {path} transport failed after retries"
                    ) from None
            except (ValueError, json.JSONDecodeError):
                raise PathogenwatchError(f"Pathogenwatch {path} returned invalid JSON") from None
            self.sleep(min(2**attempt, 30))

    def supported_organisms(self) -> list[dict]:
        response = self.request_json("GET", "/api/organisms/supported")
        if not isinstance(response, list):
            raise PathogenwatchError("Unsupported organisms response shape")
        return response

    def search(
        self, organism_id: str, st: str, mlst_scheme: str = "mlst", *, page_size: int = 100
    ) -> dict:
        if mlst_scheme not in {"mlst", "mlst2"}:
            raise PathogenwatchError(
                "MLST scheme must be the explicit Pathogenwatch field mlst or mlst2"
            )
        if not organism_id or not str(st).strip() or page_size < 1:
            raise PathogenwatchError("Organism, ST and positive page size are required")
        query = {"organismId": str(organism_id), mlst_scheme: [str(st)], "qc": [True, False]}
        pages, genomes, ids, cursors = [], [], set(), set()
        cursor, total = "", None
        for _ in range(100000):
            params = {"limit": page_size, "sort": "id"}
            if cursor:
                params["after"] = cursor
            page = self.request_json("POST", "/api/search/genomes", body=query, params=params)
            if (
                not isinstance(page, dict)
                or not isinstance(page.get("genomes"), list)
                or not isinstance(page.get("meta"), dict)
            ):
                raise PathogenwatchError("Unsupported search response shape")
            meta = page["meta"]
            count = meta.get("count")
            if type(count) is not int or count < 0 or (total is not None and total != count):
                raise PathogenwatchError("Search count missing or changed during pagination")
            total = count
            pages.append(page)
            for row in page["genomes"]:
                if not row.get("uuid") or row.get("projectAccess") != "PUBLIC":
                    raise PathogenwatchError("Search returned missing UUID or non-public record")
                if row["uuid"] in ids:
                    raise PathogenwatchError(
                        "Repeated genome across search pages; incomplete snapshot"
                    )
                ids.add(row["uuid"])
                genomes.append(row)
            if len(genomes) == total:
                break
            if len(genomes) > total or not page["genomes"]:
                raise PathogenwatchError("Search pagination did not reconcile with total count")
            cursor = meta.get("endCursor")
            if not cursor or cursor in cursors:
                raise PathogenwatchError("Search pagination has missing or repeated cursor")
            cursors.add(cursor)
        else:
            raise PathogenwatchError("Search pagination safety limit exceeded")
        return {
            "query": query,
            "raw_search_pages": pages,
            "genomes": genomes,
            "expected_count": total,
            "complete": True,
        }

    def freeze_catalogue(
        self,
        path: Path | str,
        organism_id: str = "573",
        st: str = "147",
        mlst_scheme: str = "mlst",
        *,
        page_size: int = 100,
    ) -> dict:
        started = time.monotonic()
        organisms = self.supported_organisms()
        organism = next(
            (r for r in organisms if str(r.get("organismId")) == str(organism_id)), None
        )
        if organism is None or "MLST" not in (organism.get("typing") or []):
            raise PathogenwatchError("Unsupported Pathogenwatch organism/MLST capability")
        result = self.search(organism_id, st, mlst_scheme, page_size=page_size)

        def fetch(search):
            detail = self.request_json("GET", "/api/genomes/details", params={"id": search["uuid"]})
            if not isinstance(detail, dict):
                raise PathogenwatchError("Unsupported details response shape")
            row = normalize_record(
                search, detail, organism_id=organism_id, st=st, mlst_scheme=mlst_scheme
            )
            return search["uuid"], detail, row

        with ThreadPoolExecutor(max_workers=max(1, min(self.workers, 16))) as pool:
            fetched = list(pool.map(fetch, result["genomes"]))
        rows = sorted([row for _, _, row in fetched], key=lambda r: r["source_genome_id"])
        details = {source_id: detail for source_id, detail, _ in fetched}
        retrieved_at = datetime.now(timezone.utc).isoformat()
        unique, audit = deduplicate_catalogue(rows)
        provenance = {
            "source": "pathogenwatch",
            "base_url": self.base_url,
            "query": result["query"],
            "query_sha256": content_hash(result["query"]),
            "retrieved_at": retrieved_at,
            "species": organism["fullName"],
            "organism_capabilities": organism,
            "service_version": "unavailable",
            "database_version": "unavailable",
            "public_only": True,
            "complete": True,
            "expected_count": result["expected_count"],
            "raw_record_count": len(rows),
            "sample_unit_count": len(unique),
            "search_page_count": len(result["raw_search_pages"]),
            "details_count": len(details),
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "raw_search_sha256": content_hash(result["raw_search_pages"]),
            "raw_details_sha256": content_hash(details),
            "rows_sha256": content_hash(rows),
            "metadata_completeness": catalogue_summary(rows),
        }
        envelope = {
            "schema_version": 1,
            "provenance": provenance,
            "rows": rows,
            "raw_search_pages": result["raw_search_pages"],
            "raw_details": details,
            "raw_supported_organisms": organisms,
            "deduplication_audit": audit,
        }
        envelope["snapshot_sha256"] = content_hash(envelope)
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        temporary.write_text(json.dumps(envelope, sort_keys=True, indent=2), encoding="utf-8")
        temporary.replace(destination)
        return envelope


def catalogue_summary(rows: list[dict]) -> dict:
    return {
        "records": len(rows),
        "qc_pass": sum(r.get("qc_pass") is True for r in rows),
        "qc_fail": sum(r.get("qc_pass") is False for r in rows),
        "qc_unknown": sum(r.get("qc_pass") is None for r in rows),
        "known_country": sum(r.get("country", "Unknown") != "Unknown" for r in rows),
        "country_conflicts": sum(bool(r.get("country_conflict")) for r in rows),
        "dated": sum(bool(r.get("dated_cohort_eligible")) for r in rows),
        "date_precision": {
            p: sum(r.get("date_precision") == p for r in rows)
            for p in ("day", "month", "year", "interval", "missing", "invalid")
        },
        "host": sum(bool(r.get("host")) for r in rows),
        "isolation_source": sum(bool(r.get("isolation_source")) for r in rows),
        "biosample": sum(bool(r.get("biosample")) for r in rows),
    }


def load_catalogue(path: Path | str) -> dict:
    """Read a complete snapshot, verify hashes, and re-normalise the raw records."""
    try:
        envelope = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise PathogenwatchError("Cannot read Pathogenwatch frozen catalogue") from error
    if not isinstance(envelope, dict) or envelope.get("schema_version") != 1:
        raise PathogenwatchError("Unsupported frozen catalogue schema")
    snapshot_hash = envelope.get("snapshot_sha256")
    body = {key: value for key, value in envelope.items() if key != "snapshot_sha256"}
    if snapshot_hash != content_hash(body):
        raise PathogenwatchError("Frozen catalogue content hash mismatch")
    provenance = envelope["provenance"]
    rows, details = envelope["rows"], envelope["raw_details"]
    if provenance.get("complete") is not True or provenance.get("public_only") is not True:
        raise PathogenwatchError("Catalogue is not explicitly complete and public-only")
    if not (
        len(rows)
        == len(details)
        == provenance.get("expected_count")
        == provenance.get("raw_record_count")
    ):
        raise PathogenwatchError("Catalogue records do not reconcile")
    if (
        content_hash(rows) != provenance.get("rows_sha256")
        or content_hash(details) != provenance.get("raw_details_sha256")
        or content_hash(envelope["raw_search_pages"]) != provenance.get("raw_search_sha256")
    ):
        raise PathogenwatchError("Frozen catalogue section hash mismatch")
    query = provenance["query"]
    scheme = "mlst" if "mlst" in query else "mlst2"
    search_rows = [r for page in envelope["raw_search_pages"] for r in page["genomes"]]
    replay = sorted(
        [
            normalize_record(
                r,
                details[r["uuid"]],
                organism_id=query["organismId"],
                st=query[scheme][0],
                mlst_scheme=scheme,
            )
            for r in search_rows
        ],
        key=lambda r: r["source_genome_id"],
    )
    if replay != rows:
        raise PathogenwatchError("Frozen normalised rows differ from raw replay")
    return envelope
