#!/usr/bin/env python3

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


REQUIRED_FIELDS = {
    "phase",
    "population_size",
    "arrangement_frequency",
    "context_size",
    "lineage_count",
    "eligible_pair_count",
    "pair_rate_used",
    "total_c_used",
}


def close_enough(observed: float, expected: float) -> bool:
    return math.isclose(observed, expected, rel_tol=1e-5, abs_tol=1e-12)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate the coalescence rates recorded by the current SMC implementation."
    )
    parser.add_argument(
        "path",
        type=Path,
        help="Run directory or smc_coalescence_diagnostics.csv path.",
    )
    args = parser.parse_args()

    csv_path = args.path
    if csv_path.is_dir():
        csv_path = csv_path / "smc_coalescence_diagnostics.csv"

    failures = Counter()
    pair_counts = Counter()
    skipped_generalized_rows = 0
    rows = 0

    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_FIELDS.difference(reader.fieldnames or [])
        if missing:
            parser.error(f"missing fields in {csv_path}: {', '.join(sorted(missing))}")

        for row in reader:
            rows += 1
            population_size = float(row["population_size"])
            arrangement_frequency = float(row["arrangement_frequency"])
            context_size = float(row["context_size"])
            phase = row["phase"]
            lineage_count = int(row["lineage_count"])
            eligible_pairs = int(row["eligible_pair_count"])
            pair_rate = float(row["pair_rate_used"])
            total_c = float(row["total_c_used"])

            pair_counts[eligible_pairs] += 1
            if population_size < 0.0:
                failures["negative population size"] += 1
            if not 0.0 <= arrangement_frequency <= 1.0:
                failures["arrangement frequency outside [0,1]"] += 1
            if lineage_count < 0:
                failures["negative lineage count"] += 1
            if eligible_pairs < 0:
                failures["negative eligible-pair count"] += 1

            expected_context_size = population_size * arrangement_frequency
            if not close_enough(context_size, expected_context_size):
                failures["context_size != population_size * arrangement_frequency"] += 1

            expected_pair_rate = 1.0 / context_size if context_size > 0.0 else 0.0
            if not close_enough(pair_rate, expected_pair_rate):
                failures["pair_rate_used != 1 / context_size"] += 1

            # A generalized above-root row could contain compatible pairs
            # from several contexts with different rates. The current
            # single-cut case has two lineages, so this identity still holds.
            if phase == "below_root" or lineage_count == 2:
                expected_total_c = eligible_pairs * pair_rate
                if not close_enough(total_c, expected_total_c):
                    failures["total_c_used != eligible_pair_count * pair_rate_used"] += 1
            else:
                skipped_generalized_rows += 1

    print(f"Validated {rows:,} rows from {csv_path}")
    print(
        "Eligible-pair counts: "
        + ", ".join(f"{count}={pair_counts[count]:,}" for count in sorted(pair_counts))
    )
    if skipped_generalized_rows:
        print(
            f"Skipped aggregate-rate identity for {skipped_generalized_rows:,} "
            "above-root rows with more than two lineages"
        )
    if failures:
        for message, count in failures.items():
            print(f"FAIL: {message}: {count:,} rows")
        return 1

    print("PASS: all coalescence-rate identities hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
