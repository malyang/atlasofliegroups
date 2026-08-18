#!/usr/bin/env python3
"""Exact combinatorial certificates for polygon Hall--skein straightening.

The library deliberately separates the theorem from the finite audit.  It
implements the type-A polygon combinatorics used in the paper:

* diagonals and boundary arcs of a convex m-gon;
* crossing detection and the two Ptolemy/Kauffman smoothings;
* recursive Laurent-polynomial straightening with coefficients w and w^{-1};
* confluence checks for different first-crossing choices;
* products of straightened multicurves;
* enumeration of triangulations.

It does not attempt to compute Hall structure constants of an arbitrary bound
quiver.  Those constants enter the paper through the exact Hall filtration and
its triangular comparison with the skein straightening system.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import comb
from typing import Callable, Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

import sympy as sp

Arc = Tuple[int, int]
Curve = Tuple[Arc, ...]
Expansion = Dict[Curve, sp.Expr]

W = sp.Symbol("w", nonzero=True)


def normalize_arc(a: int, b: int, m: int) -> Arc:
    if not (0 <= a < m and 0 <= b < m):
        raise ValueError(f"arc endpoints {(a, b)} outside 0..{m-1}")
    if a == b:
        raise ValueError("an arc needs two distinct endpoints")
    return (a, b) if a < b else (b, a)


def canonical_curve(arcs: Iterable[Arc], m: int) -> Curve:
    return tuple(sorted(normalize_arc(a, b, m) for a, b in arcs))


def is_boundary(arc: Arc, m: int) -> bool:
    a, b = normalize_arc(*arc, m)
    return b == a + 1 or (a == 0 and b == m - 1)


def is_diagonal(arc: Arc, m: int) -> bool:
    return not is_boundary(arc, m)


def all_boundary_arcs(m: int) -> Tuple[Arc, ...]:
    if m < 3:
        raise ValueError("a polygon has at least three vertices")
    return tuple((i, i + 1) for i in range(m - 1)) + ((0, m - 1),)


def all_diagonals(m: int) -> Tuple[Arc, ...]:
    return tuple(
        (a, b)
        for a in range(m)
        for b in range(a + 1, m)
        if is_diagonal((a, b), m)
    )


def crossing(a: Arc, b: Arc, m: int) -> bool:
    """Return whether two straight chords cross in the polygon interior."""
    a = normalize_arc(*a, m)
    b = normalize_arc(*b, m)
    if len(set(a + b)) < 4:
        return False
    e0, e1, e2, e3 = sorted(a + b)
    return {a, b} == {(e0, e2), (e1, e3)}


def crossing_pairs(curve: Curve, m: int) -> Tuple[Tuple[int, int], ...]:
    curve = canonical_curve(curve, m)
    return tuple(
        (i, j)
        for i in range(len(curve))
        for j in range(i + 1, len(curve))
        if crossing(curve[i], curve[j], m)
    )


def crossing_number(curve: Curve, m: int) -> int:
    return len(crossing_pairs(curve, m))


def smoothings(a: Arc, b: Arc, m: int) -> Tuple[Tuple[Arc, Arc], Tuple[Arc, Arc]]:
    """The two non-crossing endpoint pairings for a crossing chord pair."""
    a = normalize_arc(*a, m)
    b = normalize_arc(*b, m)
    if not crossing(a, b, m):
        raise ValueError(f"arcs {a} and {b} do not cross")
    e0, e1, e2, e3 = sorted(a + b)
    s0 = (normalize_arc(e0, e1, m), normalize_arc(e2, e3, m))
    s1 = (normalize_arc(e0, e3, m), normalize_arc(e1, e2, m))
    return s0, s1


def _clean(expansion: Mapping[Curve, sp.Expr]) -> Expansion:
    result: Expansion = {}
    for curve, coefficient in expansion.items():
        value = sp.expand(coefficient)
        if value != 0:
            result[curve] = value
    return dict(sorted(result.items()))


def add_expansions(*expansions: Mapping[Curve, sp.Expr]) -> Expansion:
    result: Expansion = {}
    for expansion in expansions:
        for curve, coefficient in expansion.items():
            result[curve] = result.get(curve, sp.Integer(0)) + coefficient
    return _clean(result)


def scale_expansion(expansion: Mapping[Curve, sp.Expr], scalar: sp.Expr) -> Expansion:
    return _clean({curve: scalar * coefficient for curve, coefficient in expansion.items()})


def _resolve_with_first(curve: Curve, m: int, pair: Tuple[int, int]) -> Expansion:
    curve = canonical_curve(curve, m)
    i, j = pair
    if not (0 <= i < j < len(curve)) or not crossing(curve[i], curve[j], m):
        raise ValueError("the requested first pair is not a crossing")
    remainder = [arc for k, arc in enumerate(curve) if k not in (i, j)]
    s0, s1 = smoothings(curve[i], curve[j], m)
    c0 = canonical_curve(remainder + list(s0), m)
    c1 = canonical_curve(remainder + list(s1), m)
    return add_expansions(
        scale_expansion(resolve(c0, m), W),
        scale_expansion(resolve(c1, m), W ** -1),
    )


@lru_cache(maxsize=None)
def resolve(curve: Curve, m: int) -> Expansion:
    """Straighten a chord multidiagram into non-crossing multicurves.

    The first crossing in lexicographic order is used.  Confluence is audited
    independently by ``assert_first_choice_confluence``.
    """
    curve = canonical_curve(curve, m)
    pairs = crossing_pairs(curve, m)
    if not pairs:
        return {curve: sp.Integer(1)}
    return _resolve_with_first(curve, m, pairs[0])


def equivalent_expansions(left: Mapping[Curve, sp.Expr], right: Mapping[Curve, sp.Expr]) -> bool:
    keys = set(left) | set(right)
    return all(sp.expand(left.get(k, 0) - right.get(k, 0)) == 0 for k in keys)


def assert_first_choice_confluence(curve: Curve, m: int) -> int:
    curve = canonical_curve(curve, m)
    pairs = crossing_pairs(curve, m)
    if len(pairs) <= 1:
        return len(pairs)
    reference = _resolve_with_first(curve, m, pairs[0])
    for pair in pairs[1:]:
        candidate = _resolve_with_first(curve, m, pair)
        if not equivalent_expansions(reference, candidate):
            raise AssertionError(
                f"non-confluent first crossing for m={m}, curve={curve}, "
                f"reference={pairs[0]}, candidate={pair}"
            )
    return len(pairs)


def multiply_expansions(left: Mapping[Curve, sp.Expr], right: Mapping[Curve, sp.Expr], m: int) -> Expansion:
    result: Expansion = {}
    for c1, a1 in left.items():
        for c2, a2 in right.items():
            resolved = resolve(canonical_curve(c1 + c2, m), m)
            for curve, coefficient in resolved.items():
                result[curve] = result.get(curve, 0) + a1 * a2 * coefficient
    return _clean(result)


def multiply_curves(left: Curve, right: Curve, m: int) -> Expansion:
    return resolve(canonical_curve(left + right, m), m)


def assert_associative(a: Curve, b: Curve, c: Curve, m: int) -> None:
    left = multiply_expansions(multiply_curves(a, b, m), {canonical_curve(c, m): 1}, m)
    right = multiply_expansions({canonical_curve(a, m): 1}, multiply_curves(b, c, m), m)
    if not equivalent_expansions(left, right):
        raise AssertionError(f"associativity failed for m={m}: {a}, {b}, {c}")


def pairwise_compatible(arcs: Sequence[Arc], m: int) -> bool:
    return all(not crossing(a, b, m) for a, b in combinations(arcs, 2))


def triangulations(m: int) -> Tuple[Tuple[Arc, ...], ...]:
    """Enumerate maximal pairwise non-crossing diagonal sets."""
    diagonals = all_diagonals(m)
    target = m - 3
    result = [
        tuple(choice)
        for choice in combinations(diagonals, target)
        if pairwise_compatible(choice, m)
    ]
    return tuple(result)


def catalan(n: int) -> int:
    if n < 0:
        raise ValueError("Catalan index must be nonnegative")
    return comb(2 * n, n) // (n + 1)


def expected_polygon_counts(m: int) -> Dict[str, int]:
    return {
        "vertices": m,
        "rank": m - 3,
        "diagonals": m * (m - 3) // 2,
        "crossing_pairs": comb(m, 4),
        "triangulations": catalan(m - 2),
    }


def curve_multisets(arcs: Sequence[Arc], degree: int) -> Iterator[Curve]:
    for indices in combinations_with_replacement_indices(len(arcs), degree):
        yield tuple(arcs[i] for i in indices)


def combinations_with_replacement_indices(n: int, r: int) -> Iterator[Tuple[int, ...]]:
    if r < 0:
        return
    if r == 0:
        yield ()
        return

    def rec(start: int, depth: int, acc: List[int]) -> Iterator[Tuple[int, ...]]:
        if depth == r:
            yield tuple(acc)
            return
        for i in range(start, n):
            acc.append(i)
            yield from rec(i, depth + 1, acc)
            acc.pop()

    yield from rec(0, 0, [])


def expansion_signature(expansion: Mapping[Curve, sp.Expr]) -> Tuple[Tuple[Curve, str], ...]:
    return tuple((curve, str(sp.expand(coeff))) for curve, coeff in sorted(expansion.items()))


@dataclass(frozen=True)
class PolygonAudit:
    vertices: int
    diagonals: int
    crossing_pairs: int
    triangulations: int
    triple_associativity_checks: int
    confluence_checks: int
    degree_four_curves: int


def audit_polygon(m: int, *, degree_four_limit: int | None = None) -> PolygonAudit:
    diagonals = all_diagonals(m)
    expected = expected_polygon_counts(m)
    crossing_count = sum(1 for a, b in combinations(diagonals, 2) if crossing(a, b, m))
    tris = triangulations(m)
    if len(diagonals) != expected["diagonals"]:
        raise AssertionError("diagonal count mismatch")
    if crossing_count != expected["crossing_pairs"]:
        raise AssertionError("crossing-pair count mismatch")
    if len(tris) != expected["triangulations"]:
        raise AssertionError("triangulation count mismatch")

    singleton = [(arc,) for arc in diagonals]
    assoc = 0
    for a in singleton:
        for b in singleton:
            for c in singleton:
                assert_associative(a, b, c, m)
                assoc += 1

    conf = 0
    degree_four = 0
    for curve in curve_multisets(diagonals, 4):
        degree_four += 1
        if degree_four_limit is not None and degree_four > degree_four_limit:
            break
        conf += assert_first_choice_confluence(curve, m)

    return PolygonAudit(
        vertices=m,
        diagonals=len(diagonals),
        crossing_pairs=crossing_count,
        triangulations=len(tris),
        triple_associativity_checks=assoc,
        confluence_checks=conf,
        degree_four_curves=min(degree_four, degree_four_limit or degree_four),
    )
