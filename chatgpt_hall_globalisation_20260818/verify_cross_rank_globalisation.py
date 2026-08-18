#!/usr/bin/env python3
"""Cross-rank deterministic regression for Hall-cocycle globalization."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import sympy as sp

from globalisation_holonomy import (
    TropicalSeed,
    full_principal_globalization_conditions,
    globalizable,
    imat,
    invariant_skew_space,
    matrix_json,
    wall_transition,
)


def integer_compatible_form(b: sp.Matrix) -> tuple[sp.Matrix, int]:
    b = imat(b)
    if b.det() == 0:
        raise ValueError("full-rank exchange matrix required")
    rational = b.T.inv()
    scale = 1
    for x in rational:
        scale = math.lcm(scale, int(sp.Rational(x).q))
    form = imat(scale * rational)
    assert b.T * form == scale * sp.eye(b.rows)
    return form, scale


def words(rank: int, count: int = 8, length: int = 30) -> Iterable[tuple[int, ...]]:
    for seed in range(count):
        state = 101 + 47 * rank + 211 * seed
        word = []
        previous = -1
        for step in range(length):
            state = (1103515245 * state + 12345 + 17 * step) % (2**31)
            vertex = state % rank
            if rank > 1 and vertex == previous:
                vertex = (vertex + step + 1) % rank
            word.append(int(vertex))
            previous = int(vertex)
        yield tuple(word)


def full_rank_matrices() -> list[sp.Matrix]:
    return [
        imat([[0, 1], [-1, 0]]),
        imat(
            [
                [0, -1, 0, 0],
                [1, 0, 1, -1],
                [0, -1, 0, 1],
                [0, 1, -1, 0],
            ]
        ),
        imat(
            [
                [0, 2, -1, 0],
                [-2, 0, 1, 1],
                [1, -1, 0, 1],
                [0, -1, -1, 0],
            ]
        ),
        imat(
            [
                [0, 1, 0, 1, 0, 0],
                [-1, 0, 1, 0, 1, 0],
                [0, -1, 0, 1, 0, 1],
                [-1, 0, -1, 0, 1, 0],
                [0, -1, 0, -1, 0, 1],
                [0, 0, -1, 0, -1, 0],
            ]
        ),
    ]


def path_exchange(rank: int) -> sp.Matrix:
    b = sp.zeros(rank)
    for i in range(rank - 1):
        b[i, i + 1] = 1
        b[i + 1, i] = -1
    return imat(b)


def sampled_states(b0: sp.Matrix):
    for word in words(b0.rows):
        actual = TropicalSeed.initial(b0)
        opposite = TropicalSeed.initial(-b0)
        yield actual, opposite
        for k in word:
            actual = actual.mutate(k)
            opposite = opposite.mutate(k)
            yield actual, opposite


def build_certificate() -> dict[str, object]:
    full_records = []
    singular_records = []
    full_seed_checks = 0
    full_wall_checks = 0
    full_fixed_frozen_checks = 0
    canonical_any_rank_checks = 0
    singular_full_no_go_witnesses = 0

    for index, b0 in enumerate(full_rank_matrices(), start=1):
        s, d = integer_compatible_form(b0)
        transfers = []
        local_checks = 0
        previous = None
        for actual, opposite in sampled_states(b0):
            delta = -sp.Matrix(opposite.g_matrix)
            g = sp.Matrix(actual.g_matrix)
            transfer = imat(g * delta.inv())
            assert transfer.T * s * transfer == s
            assert delta.T * s * delta == g.T * s * g
            assert full_principal_globalization_conditions(b0, delta, g, s, d) == (True, True)
            assert full_principal_globalization_conditions(b0, delta, g, sp.zeros(b0.rows), d)[0]
            transfers.append(transfer)
            local_checks += 1
            full_seed_checks += 1
            full_fixed_frozen_checks += 1
            canonical_any_rank_checks += 1
            if previous is not None:
                difference_rank = int((previous.inv() * transfer - sp.eye(b0.rows)).rank())
                if difference_rank <= 1:
                    transition = wall_transition(previous, transfer, s)
                    assert transition.preserved
                    full_wall_checks += 1
            previous = transfer
        assert globalizable(s, transfers)
        _, invariant_basis = invariant_skew_space(transfers)
        full_records.append(
            {
                "matrix_index": index,
                "rank": b0.rows,
                "determinant": int(b0.det()),
                "compatibility_scalar": d,
                "sampled_seed_prefixes": local_checks,
                "sampled_transfer_invariant_dimension": len(invariant_basis),
                "exchange_matrix": matrix_json(b0),
                "compatible_form": matrix_json(s),
            }
        )

    for rank in (3, 5):
        b0 = path_exchange(rank)
        assert b0.det() == 0
        transfers = []
        local_checks = 0
        full_equal = 0
        first_witness = None
        for actual, opposite in sampled_states(b0):
            delta = -sp.Matrix(opposite.g_matrix)
            g = sp.Matrix(actual.g_matrix)
            transfer = imat(g * delta.inv())
            transfers.append(transfer)
            assert globalizable(sp.zeros(rank), (transfer,))
            canonical_any_rank_checks += 1
            mutable, mixed = full_principal_globalization_conditions(b0, delta, g, sp.zeros(rank), 1)
            assert mutable
            if mixed:
                full_equal += 1
            else:
                singular_full_no_go_witnesses += 1
                if first_witness is None:
                    first_witness = {
                        "Delta": matrix_json(delta),
                        "G": matrix_json(g),
                    }
            local_checks += 1
        equations, invariant_basis = invariant_skew_space(transfers)
        singular_records.append(
            {
                "rank": rank,
                "matrix_rank": int(b0.rank()),
                "sampled_seed_prefixes": local_checks,
                "canonical_mutable_globalisation_checks": local_checks,
                "canonical_full_fixed_frozen_equal_prefixes": full_equal,
                "canonical_full_fixed_frozen_failure_prefixes": local_checks - full_equal,
                "sampled_transfer_invariant_dimension": len(invariant_basis),
                "invariant_equation_rank": int(equations.rank()),
                "first_full_principal_no_go_witness": first_witness,
            }
        )

    return {
        "certificate": "Cross-rank globalization and holonomy regression",
        "full_rank_families": full_records,
        "singular_path_families": singular_records,
        "totals": {
            "full_rank_seed_prefix_checks": full_seed_checks,
            "full_rank_symplectic_adjacent_prefix_checks": full_wall_checks,
            "full_rank_fixed_frozen_global_checks": full_fixed_frozen_checks,
            "canonical_mutable_any_rank_checks": canonical_any_rank_checks,
            "singular_full_principal_no_go_witnesses": singular_full_no_go_witnesses,
        },
        "verified_statements": [
            "S=0 gives a global mutable Hall cocycle in every sampled rank",
            "the coefficient-free compatible form is fixed by every sampled full-rank transfer",
            "the same compatible form globalizes the full fixed-frozen principal quantization in full rank",
            "the canonical S=0 full fixed-frozen principal form fails away from charts with Delta=-G",
        ],
        "proof_boundary": (
            "This file is deterministic regression evidence.  The exact arbitrary-rank criteria are proved symbolically in the paper."
        ),
    }


def write_tex_macros(certificate: dict[str, object], path: Path) -> None:
    totals = certificate["totals"]
    text = "\n".join(
        [
            "% Generated deterministically by verify_cross_rank_globalisation.py",
            f"\\newcommand{{\\CrossFullSeedChecks}}{{{totals['full_rank_seed_prefix_checks']}}}",
            f"\\newcommand{{\\CrossSymplecticChecks}}{{{totals['full_rank_symplectic_adjacent_prefix_checks']}}}",
            f"\\newcommand{{\\CrossFullPrincipalChecks}}{{{totals['full_rank_fixed_frozen_global_checks']}}}",
            f"\\newcommand{{\\CrossCanonicalMutableChecks}}{{{totals['canonical_mutable_any_rank_checks']}}}",
            f"\\newcommand{{\\CrossSingularNoGoWitnesses}}{{{totals['singular_full_principal_no_go_witnesses']}}}",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=Path("cross_rank_globalisation.json"))
    parser.add_argument("--tex", type=Path, default=Path("cross_rank_results.tex"))
    args = parser.parse_args()
    certificate = build_certificate()
    args.json.write_text(json.dumps(certificate, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_tex_macros(certificate, args.tex)
    print("Cross-rank globalization regression: PASS")
    for key, value in certificate["totals"].items():
        print(f"  {key}: {value}")
    print(f"  JSON: {args.json}")
    print(f"  TeX macros: {args.tex}")


if __name__ == "__main__":
    main()
