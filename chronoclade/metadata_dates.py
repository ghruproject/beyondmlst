"""Collection-date intervals for analysis and shared tree presentations."""

import calendar
from collections.abc import Mapping
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


def _decimal_year(day: date) -> float:
    return day.year + (day - date(day.year, 1, 1)).days / (
        366 if calendar.isleap(day.year) else 365
    )


def _bounds(value: str) -> tuple[date, date, str]:
    if re.fullmatch(r"\d{4}", value):
        year = int(value)
        return date(year, 1, 1), date(year, 12, 31), "year"
    if re.fullmatch(r"\d{4}-\d{2}", value):
        year, month = map(int, value.split("-"))
        return date(year, month, 1), date(year, month, calendar.monthrange(year, month)[1]), "month"
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        day = date.fromisoformat(value)
        return day, day, "day"
    raise ValueError("unsupported collection date")


def sample_date_interval(sample: Mapping) -> dict:
    """Read canonical bounds first, preserving month/year/interval uncertainty.

    Intervals entirely after today are excluded. A current month/year interval
    retains its full bounds rather than silently acquiring day precision.
    """
    result = dict.fromkeys(("start", "end", "precision", "year_min", "year_max",
                          "date_lower", "date_upper", "midpoint_year", "year"))
    precision = sample.get("date_precision")
    raw = str(sample.get("collection_date") or sample.get("collection_year") or "").strip()
    if precision in {"missing", "invalid"}:
        return {**result, "precision": precision, "status": precision,
                "reason": f"{precision}_collection_date"}
    try:
        if sample.get("date_start") or sample.get("date_end"):
            start, end = (date.fromisoformat(sample[key]) for key in ("date_start", "date_end"))
            if precision is None:
                precision = "day" if start == end else "interval"
            if precision not in {"day", "month", "year", "interval"}:
                raise ValueError("invalid date precision")
            if precision == "day" and start != end:
                raise ValueError("day precision has differing bounds")
            if precision == "month" and (
                start.day != 1 or (start.year, start.month) != (end.year, end.month)
                or end.day != calendar.monthrange(end.year, end.month)[1]
            ):
                raise ValueError("month precision does not cover a full month")
            if precision == "year" and (
                start != date(start.year, 1, 1) or end != date(start.year, 12, 31)
            ):
                raise ValueError("year precision does not cover a full year")
        elif raw.casefold() in _MISSING:
            return {**result, "precision": "missing", "status": "missing",
                    "reason": "missing_collection_date"}
        elif "/" in raw:
            first, last = raw.split("/")
            start, _, _ = _bounds(first.strip())
            _, end, _ = _bounds(last.strip())
            precision = "interval"
        else:
            start, end, precision = _bounds(raw)
        if start > end:
            raise ValueError("date interval is reversed")
    except (ValueError, TypeError, KeyError, OverflowError):
        return {**result, "precision": "invalid", "status": "invalid",
                "reason": "invalid_collection_date"}
    lower, upper = _decimal_year(start), _decimal_year(end)
    future = start > date.today()
    return {
        "status": "future" if future else "valid",
        "reason": "future_collection_date" if future else None,
        "start": start.isoformat(), "end": end.isoformat(), "precision": precision,
        "year_min": lower, "year_max": upper, "date_lower": lower, "date_upper": upper,
        "midpoint_year": (lower + upper) / 2, "year": start.year,
    }
