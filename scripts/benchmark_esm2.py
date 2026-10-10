#!/usr/bin/env python3
"""Measure local ESM-2 mean-pooled embedding runtime without network in timed phases.

Use an isolated environment with torch, fair-esm, numpy, and psutil. Download
checkpoints separately with --download; each model/device run should be a fresh
process so its sampled memory is independent. Timing includes token conversion,
device transfer, inference, mean pooling (residue tokens only), and CPU copy.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import statistics
import subprocess
import threading
import time
import urllib.request
from pathlib import Path

MODELS = {"8M": "esm2_t6_8M_UR50D", "35M": "esm2_t12_35M_UR50D"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path):
    records, name, parts = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name, parts = line[1:].split()[0], []
        elif line.strip():
            parts.append(line.strip().upper())
    if name is not None:
        records.append((name, "".join(parts)))
    if not records or any(not seq for _, seq in records):
        raise ValueError("Nonempty protein FASTA required")
    unique, seen = [], set()
    for name, seq in records:
        if seq not in seen:
            seen.add(seq)
            unique.append((name, seq))
    return records, unique


def download(model, cache):
    name = MODELS[model]
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / (name + ".pt")
    url = f"https://dl.fbaipublicfiles.com/fair-esm/models/{name}.pt"
    cached = path.exists()
    start = time.perf_counter()
    if not cached:
        temporary = path.with_suffix(".part")
        with urllib.request.urlopen(url, timeout=90) as source, temporary.open("wb") as target:
            while block := source.read(1024 * 1024):
                target.write(block)
        temporary.replace(path)
    return {
        "model": model,
        "url": url,
        "path": str(path.resolve()),
        "cache_hit": cached,
        "download_or_cache_check_seconds": time.perf_counter() - start,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def benchmark(args):
    import esm
    import numpy as np
    import psutil
    import torch

    torch.set_num_threads(args.cpu_threads)
    torch.set_num_interop_threads(1)
    if args.device == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MPS is unavailable on this host")
    records, unique = read_fasta(args.fasta)
    lengths = [len(s) for _, s in unique]
    ordered = sorted(unique, key=lambda item: (len(item[1]), item[0]))
    batches, current, max_length = [], [], 0
    for record in ordered:
        size = len(record[1]) + 2
        if size > args.token_budget:
            raise ValueError("A sequence exceeds the padded token budget")
        if current and max(max_length, size) * (len(current) + 1) > args.token_budget:
            batches.append(current)
            current, max_length = [], 0
        current.append(record)
        max_length = max(max_length, size)
    if current:
        batches.append(current)

    process = psutil.Process()
    phase = ["load"]
    peaks = {"load": {}, "warmup": {}, "timed": {}, "correctness": {}}
    done = threading.Event()

    def sample():
        values = {"rss_bytes": process.memory_info().rss}
        if args.device == "mps":
            values["mps_allocated_bytes"] = torch.mps.current_allocated_memory()
            values["mps_driver_allocated_bytes"] = torch.mps.driver_allocated_memory()
        target = peaks[phase[0]]
        for key, value in values.items():
            target[key] = max(target.get(key, 0), value)

    def monitor():
        while not done.wait(0.02):
            sample()

    def sync():
        if args.device == "mps":
            torch.mps.synchronize()

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    start = time.perf_counter()
    model_data = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model, alphabet = esm.pretrained.load_model_and_alphabet_core(
        MODELS[args.model], model_data, None
    )
    del model_data
    model.eval().to(args.device)
    sync()
    load_seconds = time.perf_counter() - start
    sample()
    converter = alphabet.get_batch_converter()
    layer = model.num_layers

    def extract(batch_list):
        embeddings = {}
        with torch.inference_mode():
            for batch in batch_list:
                _, _, tokens = converter(batch)
                result = model(tokens.to(args.device), repr_layers=[layer], return_contacts=False)
                representation = result["representations"][layer]
                vectors = torch.stack(
                    [
                        representation[i, 1 : len(seq) + 1].mean(0)
                        for i, (_, seq) in enumerate(batch)
                    ]
                )
                vectors = vectors.cpu().numpy()
                for (name, _), vector in zip(batch, vectors):
                    embeddings[name] = vector
        return (
            np.stack([embeddings[name] for name, _ in unique])
            if batch_list == batches
            else embeddings
        )

    phase[0] = "warmup"
    start = time.perf_counter()
    for _ in range(args.warmups):
        extract(batches)
        sync()
    warmup_seconds = time.perf_counter() - start
    phase[0] = "timed"
    durations, vectors, repeat_difference = [], None, 0.0
    for repeat in range(args.repeats):
        sync()
        start = time.perf_counter()
        value = extract(batches)
        sync()
        duration = time.perf_counter() - start
        durations.append(duration)
        if vectors is not None:
            repeat_difference = max(repeat_difference, float(np.max(np.abs(value - vectors))))
        vectors = value
        sample()
        print(
            json.dumps(
                {
                    "event": "repeat",
                    "model": args.model,
                    "device": args.device,
                    "repeat": repeat + 1,
                    "seconds": duration,
                }
            ),
            flush=True,
        )

    phase[0] = "correctness"
    first_name, first = unique[0]
    changed = first[:10] + ("A" if first[10] != "A" else "G") + first[11:]
    check = extract(
        [[(first_name, first), ("exact_duplicate", first), ("synthetic_variant", changed)]]
    )
    sync()
    duplicate_difference = float(np.max(np.abs(check[first_name] - check["exact_duplicate"])))
    variant_difference = float(np.linalg.norm(check[first_name] - check["synthetic_variant"]))
    sample()
    done.set()
    thread.join()
    np.savez_compressed(
        args.output.with_suffix(".npz"),
        ids=np.array([name for name, _ in unique]),
        embeddings=vectors,
    )
    info = {
        "model": args.model,
        "model_name": MODELS[args.model],
        "device": args.device,
        "checkpoint": str(args.checkpoint.resolve()),
        "checkpoint_sha256": sha256(args.checkpoint),
        "checkpoint_bytes": args.checkpoint.stat().st_size,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cpu": (
                subprocess.check_output(
                    ["sysctl", "-n", "machdep.cpu.brand_string"], text=True
                ).strip()
                if platform.system() == "Darwin"
                else platform.processor()
            ),
            "architecture": platform.machine(),
            "logical_cpu_count": psutil.cpu_count(),
            "physical_memory_bytes": psutil.virtual_memory().total,
            "versions": {
                p: importlib.metadata.version(p) for p in ["torch", "fair-esm", "numpy", "psutil"]
            },
            "torch_cpu_threads": torch.get_num_threads(),
            "interop_threads": torch.get_num_interop_threads(),
            "mps_available": torch.backends.mps.is_available(),
            "mps_fallback_env": os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK"),
        },
        "input": {
            "path": str(args.fasta.resolve()),
            "sha256": sha256(args.fasta),
            "records": len(records),
            "unique_sequences": len(unique),
            "residues": sum(lengths),
            "length_min": min(lengths),
            "length_median": statistics.median(lengths),
            "length_mean": statistics.mean(lengths),
            "length_max": max(lengths),
            "lengths": lengths,
            "token_budget": args.token_budget,
            "batches": len(batches),
            "padded_tokens": sum((max(len(s) for _, s in b) + 2) * len(b) for b in batches),
        },
        "timing": {
            "model_load_seconds": load_seconds,
            "warmup_full_cohort_count": args.warmups,
            "warmup_seconds": warmup_seconds,
            "repeat_seconds": durations,
            "median_seconds": statistics.median(durations),
            "min_seconds": min(durations),
            "proteins_per_second_median": len(unique) / statistics.median(durations),
            "proteins_per_second_min": min(len(unique) / x for x in durations),
            "residues_per_second_median": sum(lengths) / statistics.median(durations),
        },
        "sampled_peak_memory": peaks,
        "memory_sample_interval_seconds": 0.02,
        "correctness": {
            "shape": list(vectors.shape),
            "all_finite": bool(np.isfinite(vectors).all()),
            "repeat_max_abs_difference": repeat_difference,
            "exact_duplicate_max_abs_difference": duplicate_difference,
            "synthetic_single_substitution_l2_difference": variant_difference,
            "variant_is_runtime_sanity_check_only": True,
        },
        "scope": "Float32 standalone fair-esm, last-layer residue-only mean pooling; no biological ranking-quality or cluster benchmark.",
    }
    args.output.write_text(json.dumps(info, indent=2))
    print(
        json.dumps(
            {
                "event": "complete",
                "output": str(args.output),
                "timing": info["timing"],
                "correctness": info["correctness"],
            }
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--device", choices=["cpu", "mps"], default="cpu")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--token-budget", type=int, default=4096)
    parser.add_argument("--cpu-threads", type=int, default=4)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        if args.download:
            if args.cache_dir is None:
                parser.error("--download requires --cache-dir")
            result = download(args.model, args.cache_dir)
            args.output.write_text(json.dumps(result, indent=2))
            print(json.dumps(result), flush=True)
        else:
            if args.fasta is None or args.checkpoint is None:
                parser.error("benchmark requires --fasta and --checkpoint")
            benchmark(args)
    except Exception as exc:
        args.output.write_text(
            json.dumps(
                {
                    "model": args.model,
                    "device": args.device,
                    "failure": f"{type(exc).__name__}: {exc}",
                },
                indent=2,
            )
        )
        raise


if __name__ == "__main__":
    main()
