"""Collection-date intervals for analysis and shared tree presentations."""

import calendar
import datetime as dt
import re
from datetime import date
from typing import Any


def date_interval(row, *, allow_future=False):
    value = str(row.get("collection_date") or row.get("collection_year") or "").strip()
    try:
        parts = value.split("-")
        year = int(parts[0])
        if len(parts) == 1:
            start, end, precision = dt.date(year, 1, 1), dt.date(year, 12, 31), "year"
        elif len(parts) == 2:
            month = int(parts[1])
            start = dt.date(year, month, 1)
            end = dt.date(year, month, calendar.monthrange(year, month)[1])
            precision = "month"
        elif len(parts) == 3:
            start = end = dt.date.fromisoformat(value)
            precision = "day"
        else:
            return None

        if not allow_future and start > dt.date.today():
            return None

        def decimal(day):
            return day.year + (day - dt.date(day.year, 1, 1)).days / (
                366 if calendar.isleap(day.year) else 365
            )

        return {
            "start": start.isoformat(),
            "end": end.isoformat(),
            "precision": precision,
            "year_min": decimal(start),
            "year_max": decimal(end),
            "year": year,
        }
    except (ValueError, TypeError, OverflowError):
        return None


_MISSING = {"", "unknown", "none", "null", "na", "n/a", "not provided", "missing", "not collected"}


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
