#!/usr/bin/env python3
"""Complete A4 certificate for globalizable Hall quantisations and fan holonomy."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from globalisation_holonomy import (
    chart_key,
    enumerate_paired_atlas,
    form_coordinates,
    full_principal_globalization_conditions,
    globalizable,
    imat,
    invariant_skew_space,
    local_mutable_form,
    matrix_json,
    principal_quantum_form,
    reduced_term_lift,
    skew_basis,
    term_restriction,
    wall_transition,
)

B0 = imat(
    [
        [0, -1, 0, 0],
        [1, 0, 1, -1],
        [0, -1, 0, 1],
        [0, 1, -1, 0],
    ]
)
LAMBDA0 = imat(
    [
        [0, -1, -1, -1],
        [1, 0, 0, 0],
        [1, 0, 0, 1],
        [1, 0, -1, 0],
    ]
)
N = 4
D = 1


def _tuple_matrix(m: sp.Matrix) -> tuple[int, ...]:
    return tuple(int(x) for x in m)


def _primitive_solution(vector: sp.Matrix) -> sp.Matrix:
    den = 1
    for x in vector:
        den = sp.ilcm(den, int(sp.Rational(x).q))
    return sp.Matrix([int(sp.Rational(x) * den) for x in vector])


def full_principal_parameter_space(charts) -> dict[str, Any]:
    basis = skew_basis(N)
    rows: list[list[int]] = []
    rhs: list[int] = []
    for state in charts:
        delta = state.delta
        g = state.g
        # Mutable invariance equations.
        transformed = [delta.T * e * delta - g.T * e * g for e in basis]
        for i in range(N):
            for j in range(i + 1, N):
                rows.append([int(x[i, j]) for x in transformed])
                rhs.append(0)
        # Fixed-frozen mixed-block equations.
        a = delta.T + g.T
        base = a * (D * sp.eye(N))
        coeffs = [a * e * B0 for e in basis]
        for i in range(N):
            for j in range(N):
                rows.append([int(x[i, j]) for x in coeffs])
                rhs.append(-int(base[i, j]))
    matrix = sp.Matrix(rows)
    vector = sp.Matrix(rhs)
    rank = int(matrix.rank())
    aug_rank = int(matrix.row_join(vector).rank())
    consistent = rank == aug_rank
    dimension = len(basis) - rank if consistent else -1
    solution_matrix = None
    solution_form = None
    if consistent:
        solution_set = sp.linsolve((matrix, vector))
        solutions = list(solution_set)
        if len(solutions) == 1 and not any(x.free_symbols for x in solutions[0]):
            solution_matrix = sp.Matrix(solutions[0])
            solution_form = sp.zeros(N)
            for coefficient, e in zip(solution_matrix, basis):
                solution_form += coefficient * e
            solution_form = imat(solution_form)
    return {
        "unknowns": len(basis),
        "equations": len(rows),
        "rank": rank,
        "augmented_rank": aug_rank,
        "consistent": consistent,
        "affine_dimension": dimension,
        "unique_solution_coordinates": [str(x) for x in solution_matrix] if solution_matrix is not None else None,
        "unique_solution_form": matrix_json(solution_form) if solution_form is not None else None,
    }


def build_certificate() -> dict[str, Any]:
    assert B0.T * LAMBDA0 == sp.eye(N)
    atlas = enumerate_paired_atlas(B0)
    assert len(atlas.states) == 1008
    assert len(atlas.chart_representatives) == 42
    assert len(atlas.oriented_edges) == 4032

    transfers_by_key: dict[tuple[int, ...], sp.Matrix] = {}
    for state in atlas.chart_representatives:
        transfer = state.transfer
        transfers_by_key.setdefault(_tuple_matrix(transfer), transfer)
    transfers = tuple(transfers_by_key.values())
    assert len(transfers) == 22

    equations, invariant_basis = invariant_skew_space(transfers)
    lambda_coordinates = form_coordinates(LAMBDA0, invariant_basis)
    assert lambda_coordinates is not None
    assert globalizable(sp.zeros(N), transfers)
    assert globalizable(LAMBDA0, transfers)

    # Verify the two equivalent chamber restrictions for every chart.
    mutable_globalisation_checks = 0
    canonical_zero_checks = 0
    lambda_checks = 0
    reduced_lift_checks = 0
    full_principal_canonical_equal_charts = 0
    full_principal_lambda_equal_charts = 0
    chart_digest = hashlib.sha256()
    witness = None
    target_delta = sp.Matrix([0, 1, 0, -1])
    target_g = sp.Matrix([1, -1, 1, 0])

    for state in atlas.chart_representatives:
        delta = state.delta
        g = state.g
        transfer = state.transfer
        lift = reduced_term_lift(delta)
        assert imat(sp.Matrix.hstack(-sp.eye(N), sp.eye(N)) * lift) == delta
        reduced_lift_checks += 1

        assert term_restriction(delta, sp.zeros(N)) == local_mutable_form(g, sp.zeros(N))
        canonical_zero_checks += 1
        assert term_restriction(delta, LAMBDA0) == local_mutable_form(g, LAMBDA0)
        lambda_checks += 1
        assert transfer.T * LAMBDA0 * transfer == LAMBDA0
        mutable_globalisation_checks += 1

        canonical_conditions = full_principal_globalization_conditions(B0, delta, g, sp.zeros(N), D)
        lambda_conditions = full_principal_globalization_conditions(B0, delta, g, LAMBDA0, D)
        if all(canonical_conditions):
            full_principal_canonical_equal_charts += 1
        if all(lambda_conditions):
            full_principal_lambda_equal_charts += 1

        for j in range(N):
            if delta[:, j] == target_delta and g[:, j] == target_g:
                witness = {
                    "split_index": [int(x) for x in target_delta],
                    "cluster_g_vector": [int(x) for x in target_g],
                    "canonical_mixed_condition": canonical_conditions[1],
                    "lambda_parameter_mixed_condition": lambda_conditions[1],
                    "chart_path": [x + 1 for x in state.path],
                }

        chart_digest.update(
            repr((chart_key(delta, g), _tuple_matrix(transfer))).encode("ascii")
        )

    assert witness is not None
    assert full_principal_canonical_equal_charts < len(atlas.chart_representatives)
    assert full_principal_lambda_equal_charts == len(atlas.chart_representatives)

    # Labeled wall transitions.  The full-rank form makes every nonzero wall a
    # symplectic transvection.  Deduplicate by unordered chart pair for the
    # 84-edge exchange graph.
    labeled_zero = 0
    labeled_rank_one = 0
    labeled_symplectic = 0
    scalar_histogram: dict[str, int] = {}
    wall_digest = hashlib.sha256()
    unique_edges: dict[tuple[object, object], tuple[int, int | None]] = {}
    for source_index, target_index, direction in atlas.oriented_edges:
        source = atlas.states[source_index]
        target = atlas.states[target_index]
        transition = wall_transition(source.transfer, target.transfer, LAMBDA0)
        rank = int((sp.Matrix(transition.relative) - sp.eye(N)).rank())
        if rank == 0:
            labeled_zero += 1
        else:
            assert rank == 1
            assert transition.preserved
            assert int(sp.Matrix(transition.relative).det()) == 1
            assert (sp.Matrix(transition.covector).T * sp.Matrix(transition.vector))[0] == 0
            labeled_rank_one += 1
            labeled_symplectic += 1
            key = str(transition.scalar)
            scalar_histogram[key] = scalar_histogram.get(key, 0) + 1
        ck1 = chart_key(source.delta, source.g)
        ck2 = chart_key(target.delta, target.g)
        edge = tuple(sorted((ck1, ck2)))
        old = unique_edges.get(edge)
        current = (rank, transition.scalar)
        if old is not None:
            assert old[0] == rank
        else:
            unique_edges[edge] = current
        wall_digest.update(
            repr((source_index, target_index, direction, rank, transition.scalar)).encode("ascii")
        )
    assert labeled_zero + labeled_rank_one == 4032
    assert len(unique_edges) == 84
    unlabeled_zero = sum(rank == 0 for rank, _ in unique_edges.values())
    unlabeled_rank_one = sum(rank == 1 for rank, _ in unique_edges.values())
    assert (unlabeled_zero, unlabeled_rank_one) == (22, 62)

    full_space = full_principal_parameter_space(atlas.chart_representatives)
    assert full_space["consistent"]
    assert full_space["affine_dimension"] == 0
    assert full_space["unique_solution_form"] == matrix_json(LAMBDA0)

    # Produce one explicit non-invariant parameter witness.
    noninvariant = None
    for index, e in enumerate(skew_basis(N), start=1):
        failures = [
            i for i, transfer in enumerate(transfers)
            if transfer.T * e * transfer != e
        ]
        if failures:
            noninvariant = {
                "standard_skew_basis_index": index,
                "form": matrix_json(e),
                "failing_transfer_count": len(failures),
            }
            break
    assert noninvariant is not None

    return {
        "certificate": "Globalizable Hall quantisations and symplectic fan holonomy for the nonhereditary A4 atlas",
        "initial_data": {
            "B0": matrix_json(B0),
            "Lambda0": matrix_json(LAMBDA0),
            "compatibility": matrix_json(B0.T * LAMBDA0),
        },
        "paired_fan": {
            "labeled_states": len(atlas.states),
            "unlabeled_chambers": len(atlas.chart_representatives),
            "oriented_mutations": len(atlas.oriented_edges),
            "distinct_transfer_matrices": len(transfers),
            "chart_sha256": chart_digest.hexdigest(),
        },
        "globalizable_mutable_parameters": {
            "ambient_skew_dimension": len(skew_basis(N)),
            "equation_rank": int(equations.rank()),
            "invariant_dimension": len(invariant_basis),
            "primitive_invariant_basis": [matrix_json(x) for x in invariant_basis],
            "Lambda0_coordinates_in_invariant_basis": [str(x) for x in lambda_coordinates],
            "canonical_principal_parameter_zero_globalizes": True,
            "coefficient_free_parameter_Lambda0_globalizes": True,
            "noninvariant_witness": noninvariant,
            "chart_restriction_checks": mutable_globalisation_checks,
            "canonical_zero_checks": canonical_zero_checks,
            "Lambda0_checks": lambda_checks,
            "reduced_term_lift_checks": reduced_lift_checks,
        },
        "wall_holonomy": {
            "labeled_zero_shears": labeled_zero,
            "labeled_rank_one_transvections": labeled_rank_one,
            "labeled_symplectic_transvection_checks": labeled_symplectic,
            "unlabeled_zero_shear_edges": unlabeled_zero,
            "unlabeled_rank_one_edges": unlabeled_rank_one,
            "transvection_scalar_histogram": scalar_histogram,
            "wall_sha256": wall_digest.hexdigest(),
        },
        "fixed_frozen_full_principal": {
            "canonical_parameter_zero_equal_chambers": full_principal_canonical_equal_charts,
            "Lambda0_parameter_equal_chambers": full_principal_lambda_equal_charts,
            "complete_affine_system": full_space,
            "explicit_C_S2_witness": witness,
            "interpretation": (
                "the mutable canonical principal cocycle globalizes in every rank, but the full fixed-frozen canonical principal form does not; "
                "for this full-rank A4 pattern the unique full fixed-frozen global parameter is Lambda0"
            ),
        },
        "proof_boundary": {
            "proved_in_paper": (
                "globalizable parameters are exactly the skew forms fixed by every transfer matrix; "
                "rank-one walls are symplectic precisely when S v is an integral multiple of the primitive wall covector; "
                "the full fixed-frozen principal form is governed by an additional mixed-block equation"
            ),
            "computed_here": (
                "the complete 1008-state/42-chamber A4 paired fan, its invariant skew lattice, all 4032 wall transitions, "
                "and the complete fixed-frozen affine parameter system"
            ),
            "not_claimed": (
                "a global Hall-to-theta homomorphism on incompatible products or a direct exact-to-critical Hall comparison"
            ),
        },
    }


def write_tex_macros(certificate: dict[str, Any], path: Path) -> None:
    fan = certificate["paired_fan"]
    inv = certificate["globalizable_mutable_parameters"]
    wall = certificate["wall_holonomy"]
    full = certificate["fixed_frozen_full_principal"]
    system = full["complete_affine_system"]
    text = "\n".join(
        [
            "% Generated deterministically by verify_A4_globalisation.py",
            f"\\newcommand{{\\AIVLabeledStates}}{{{fan['labeled_states']}}}",
            f"\\newcommand{{\\AIVChambers}}{{{fan['unlabeled_chambers']}}}",
            f"\\newcommand{{\\AIVOrientedMutations}}{{{fan['oriented_mutations']}}}",
            f"\\newcommand{{\\AIVTransferMatrices}}{{{fan['distinct_transfer_matrices']}}}",
            f"\\newcommand{{\\AIVInvariantDimension}}{{{inv['invariant_dimension']}}}",
            f"\\newcommand{{\\AIVLabeledZeroWalls}}{{{wall['labeled_zero_shears']}}}",
            f"\\newcommand{{\\AIVLabeledTransvections}}{{{wall['labeled_rank_one_transvections']}}}",
            f"\\newcommand{{\\AIVUnlabeledZeroWalls}}{{{wall['unlabeled_zero_shear_edges']}}}",
            f"\\newcommand{{\\AIVUnlabeledTransvections}}{{{wall['unlabeled_rank_one_edges']}}}",
            f"\\newcommand{{\\AIVCanonicalFullCharts}}{{{full['canonical_parameter_zero_equal_chambers']}}}",
            f"\\newcommand{{\\AIVLambdaFullCharts}}{{{full['Lambda0_parameter_equal_chambers']}}}",
            f"\\newcommand{{\\AIVFullParameterRank}}{{{system['rank']}}}",
            f"\\newcommand{{\\AIVFullParameterDimension}}{{{system['affine_dimension']}}}",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=Path("A4_globalisation_certificate.json"))
    parser.add_argument("--tex", type=Path, default=Path("A4_globalisation_results.tex"))
    args = parser.parse_args()
    certificate = build_certificate()
    args.json.write_text(json.dumps(certificate, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_tex_macros(certificate, args.tex)
    print("A4 globalisation/holonomy certificate: PASS")
    print(f"  labeled states: {certificate['paired_fan']['labeled_states']}")
    print(f"  rigid chambers: {certificate['paired_fan']['unlabeled_chambers']}")
    print(f"  invariant skew dimension: {certificate['globalizable_mutable_parameters']['invariant_dimension']}")
    print(f"  rank-one symplectic walls: {certificate['wall_holonomy']['labeled_rank_one_transvections']}")
    print(f"  fixed-frozen parameter dimension: {certificate['fixed_frozen_full_principal']['complete_affine_system']['affine_dimension']}")
    print(f"  JSON: {args.json}")
    print(f"  TeX macros: {args.tex}")


if __name__ == "__main__":
    main()
