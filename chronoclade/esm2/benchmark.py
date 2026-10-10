"""Runtime comparison, separate from extraction and biological validation."""

import argparse
import gc
from pathlib import Path
import platform
import resource
import statistics
import sys
import time

from .engine import run_embeddings
from .proteins import EmbeddingError, file_sha256
from .runtime import ESMBackend, model_spec
from .storage import write_json


def benchmark_embeddings(fasta, out, *, models=("8M", "35M"), devices=("cpu",),
                         checkpoints=None, allow_download=False, token_budget=4096,
                         repeats=3):
    """Compare models on exactly the same proteins, keeping cache reuse separate.

    The first timed load includes optional dependency imports in this process;
    later loads have warm imports and potentially warm operating-system file pages.
    No biological adequacy or genome-neighbour agreement is measured here.
    """
    if repeats < 1:
        raise EmbeddingError("repeats must be at least 1")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    report = {"schema": "chronoclade.esm2.runtime-benchmark", "schema_version": 1,
              "status": "running", "input_sha256": file_sha256(Path(fasta)),
              "hardware": {"platform": platform.platform(), "machine": platform.machine(),
                           "processor": platform.processor()},
              "parameters": {"models": list(models), "devices": list(devices),
                             "token_budget": token_budget, "repeats": repeats},
              "limitations": ["Runtime evidence only; no biological validation",
                              "Checkpoint loads after first run reuse imported libraries",
                              "Operating-system file cache is not cleared between runs",
                              "Process peak RSS is cumulative within the benchmark process"],
              "runs": rows}
    for model in models:
        spec = model_spec(model)
        checkpoint = (checkpoints or {}).get(model, (checkpoints or {}).get(spec["name"]))
        for device in devices:
            namespace = out / f"{model}-{device}"
            try:
                backend = ESMBackend(model=model, checkpoint=checkpoint,
                                     allow_download=allow_download, device=device)
                # An unmeasured inference pass warms kernel/operator initialisation.
                run_embeddings(fasta, namespace / "warmup", model=model, device=device,
                               token_budget=token_budget, _backend=backend)
                runs = []
                for repeat in range(repeats):
                    result = run_embeddings(fasta, namespace / f"fresh-{repeat + 1}", model=model,
                                            device=device, token_budget=token_budget,
                                            _backend=backend)
                    runs.append(result.manifest["runtime"])
                cache_dir = namespace / "vector-cache"
                populated = run_embeddings(fasta, namespace / "cache-population", model=model,
                                           device=device, token_budget=token_budget,
                                           cache_dir=cache_dir, _backend=backend)
                backend.synchronize()
                started = time.perf_counter()
                reused = run_embeddings(fasta, namespace / "cached-reuse", model=model,
                                        device=device, token_budget=token_budget,
                                        cache_dir=cache_dir, _backend=backend)
                backend.synchronize()
                cache_seconds = time.perf_counter() - started
                durations = [run["inference_seconds"] for run in runs]
                residues = runs[0]["residues_embedded"]
                peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                row = {"model": spec["name"], "requested_device": device, "status": "complete",
                       "provenance": result.manifest["provenance"],
                       "checkpoint_load_seconds": backend.load_seconds,
                       "sequence_lengths": result.manifest["sequence_lengths"],
                       "counts": result.manifest["counts"], "fresh_inference_runs": runs,
                       "median_inference_seconds": statistics.median(durations),
                       "median_residues_per_second": residues / statistics.median(durations),
                       "process_peak_rss_bytes": peak_rss if sys.platform == "darwin"
                       else peak_rss * 1024,
                       "accelerator_memory": backend.accelerator_memory(),
                       "cached_reuse_wall_seconds": cache_seconds,
                       "cached_reuse": reused.manifest["runtime"],
                       "cache_population": populated.manifest["runtime"]}
                del backend
                gc.collect()
            except EmbeddingError as exc:
                row = {"model": spec["name"], "requested_device": device, "status": "failed",
                       "reason": str(exc)}
            rows.append(row)
            write_json(out / "benchmark.json", report)
    report["status"] = "complete" if all(r["status"] == "complete" for r in rows) else "partial"
    write_json(out / "benchmark.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--checkpoint-8m", type=Path)
    parser.add_argument("--checkpoint-35m", type=Path)
    parser.add_argument("--allow-download", action="store_true")
    parser.add_argument("--devices", nargs="+", default=["cpu"], choices=["cpu", "mps", "cuda"])
    parser.add_argument("--token-budget", type=int, default=4096)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    report = benchmark_embeddings(args.fasta, args.out, devices=args.devices,
                                  checkpoints={"8M": args.checkpoint_8m,
                                               "35M": args.checkpoint_35m},
                                  allow_download=args.allow_download,
                                  token_budget=args.token_budget, repeats=args.repeats)
    print(args.out / "benchmark.json")
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
