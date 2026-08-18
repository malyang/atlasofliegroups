#!/usr/bin/env python3
"""Exact tools for globalizing Hall cocycles on rigid presentation fans.

The module is self-contained apart from SymPy.  It implements:

* skew-symmetric cluster mutation and tropical C/G-matrices;
* paired mutation of B and -B, producing presentation matrices Delta_R and
  cluster g-matrices G_R;
* the chamber transfer L_R = G_R Delta_R^{-1};
* the invariant lattice of skew forms S satisfying L_R^T S L_R = S;
* rank-one wall factorization and the symplectic-transvection criterion;
* mutable and full-principal globalization tests.

All arithmetic is exact.  No floating point or network access is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections import deque
from functools import reduce
from math import gcd
from typing import Iterable, Sequence

import sympy as sp


def imat(values: Sequence[Sequence[int]] | sp.Matrix) -> sp.Matrix:
    m = sp.Matrix(values)
    for x in m:
        if sp.Rational(x).q != 1:
            raise ValueError("nonintegral matrix entry")
    return m.applyfunc(int)


def require_skew(m: sp.Matrix) -> None:
    if m.rows != m.cols or m.T != -m:
        raise ValueError("matrix must be square and skew-symmetric")


def positive_part(x: int) -> int:
    return max(int(x), 0)


def mutate_exchange(b: sp.Matrix, k: int) -> sp.Matrix:
    """Fomin--Zelevinsky mutation of an m-by-n extended exchange matrix."""
    b = imat(b)
    m, n = b.shape
    if m < n or not 0 <= k < n:
        raise ValueError("invalid mutation direction")
    out = sp.zeros(m, n)
    for i in range(m):
        for j in range(n):
            if i == k or j == k:
                out[i, j] = -b[i, j]
            else:
                out[i, j] = (
                    b[i, j]
                    + positive_part(b[i, k]) * positive_part(b[k, j])
                    - positive_part(-b[i, k]) * positive_part(-b[k, j])
                )
    return imat(out)


def c_vector_sign(c: sp.Matrix, k: int) -> int:
    col = [int(c[i, k]) for i in range(c.rows)]
    if all(x >= 0 for x in col):
        return 1
    if all(x <= 0 for x in col):
        return -1
    raise ValueError(f"c-vector is not sign coherent: {col}")


def g_mutation_matrix(b: sp.Matrix, k: int, epsilon: int) -> sp.Matrix:
    b = imat(b)
    n = b.rows
    if b.cols != n or epsilon not in (-1, 1):
        raise ValueError("invalid g-mutation data")
    e = sp.eye(n)
    for i in range(n):
        e[i, k] = -1 if i == k else positive_part(-epsilon * b[i, k])
    return imat(e)


@dataclass(frozen=True)
class TropicalSeed:
    exchange: sp.ImmutableMatrix
    c_matrix: sp.ImmutableMatrix
    g_matrix: sp.ImmutableMatrix

    @classmethod
    def initial(cls, b0: sp.Matrix) -> "TropicalSeed":
        b0 = imat(b0)
        require_skew(b0)
        n = b0.rows
        return cls(sp.ImmutableMatrix(b0), sp.ImmutableMatrix(sp.eye(n)), sp.ImmutableMatrix(sp.eye(n)))

    @property
    def rank(self) -> int:
        return self.exchange.rows

    def extended_exchange(self) -> sp.Matrix:
        return sp.Matrix.vstack(sp.Matrix(self.exchange), sp.Matrix(self.c_matrix))

    def mutate(self, k: int) -> "TropicalSeed":
        b = sp.Matrix(self.exchange)
        c = sp.Matrix(self.c_matrix)
        g = sp.Matrix(self.g_matrix)
        eps = c_vector_sign(c, k)
        e = g_mutation_matrix(b, k, eps)
        ext = mutate_exchange(sp.Matrix.vstack(b, c), k)
        n = b.rows
        return TropicalSeed(
            sp.ImmutableMatrix(ext[:n, :]),
            sp.ImmutableMatrix(ext[n:, :]),
            sp.ImmutableMatrix(imat(g * e)),
        )

    def key(self) -> tuple[int, ...]:
        return tuple(int(x) for x in self.exchange) + tuple(int(x) for x in self.c_matrix)


@dataclass(frozen=True)
class PairedState:
    actual: TropicalSeed
    opposite: TropicalSeed
    path: tuple[int, ...]

    @property
    def delta(self) -> sp.Matrix:
        return -sp.Matrix(self.opposite.g_matrix)

    @property
    def g(self) -> sp.Matrix:
        return sp.Matrix(self.actual.g_matrix)

    @property
    def transfer(self) -> sp.Matrix:
        return imat(self.g * self.delta.inv())

    @property
    def local_exchange(self) -> sp.Matrix:
        return sp.Matrix(self.actual.exchange)


@dataclass(frozen=True)
class PairedAtlas:
    states: tuple[PairedState, ...]
    chart_representatives: tuple[PairedState, ...]
    oriented_edges: tuple[tuple[int, int, int], ...]


def _column_tuple(m: sp.Matrix, j: int) -> tuple[int, ...]:
    return tuple(int(m[i, j]) for i in range(m.rows))


def chart_key(delta: sp.Matrix, g: sp.Matrix) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    return tuple(sorted((_column_tuple(delta, j), _column_tuple(g, j)) for j in range(delta.cols)))


def enumerate_paired_atlas(b0: sp.Matrix, max_states: int | None = None) -> PairedAtlas:
    """Enumerate the finite paired pattern for B0 and -B0.

    ``max_states`` may be used for deterministic finite regressions in infinite
    type; when omitted, enumeration continues until the queue is exhausted.
    """
    b0 = imat(b0)
    require_skew(b0)
    actual0 = TropicalSeed.initial(b0)
    opposite0 = TropicalSeed.initial(-b0)
    queue: deque[PairedState] = deque([PairedState(actual0, opposite0, tuple())])
    seen: dict[tuple[int, ...], int] = {}
    states: list[PairedState] = []
    edges: list[tuple[int, int, int]] = []
    chart_map: dict[tuple[tuple[tuple[int, ...], tuple[int, ...]], ...], PairedState] = {}

    while queue:
        state = queue.popleft()
        key = state.actual.key()
        if key in seen:
            continue
        if sp.Matrix(state.opposite.exchange) != -sp.Matrix(state.actual.exchange):
            raise AssertionError("opposite exchange patterns diverged")
        idx = len(states)
        seen[key] = idx
        states.append(state)
        ck = chart_key(state.delta, state.g)
        old = chart_map.get(ck)
        if old is None or len(state.path) < len(old.path):
            chart_map[ck] = state
        if max_states is not None and len(states) >= max_states:
            continue
        for k in range(b0.rows):
            target = PairedState(state.actual.mutate(k), state.opposite.mutate(k), state.path + (k,))
            tkey = target.actual.key()
            if tkey not in seen:
                queue.append(target)

    # Record edges after all retained states are known.  For truncated patterns,
    # only retain edges whose endpoints were enumerated.
    for i, state in enumerate(states):
        for k in range(b0.rows):
            target_key = state.actual.mutate(k).key()
            if target_key in seen:
                edges.append((i, seen[target_key], k))

    charts = tuple(sorted(chart_map.values(), key=lambda s: chart_key(s.delta, s.g)))
    return PairedAtlas(tuple(states), charts, tuple(edges))


def skew_basis(n: int) -> tuple[sp.Matrix, ...]:
    basis: list[sp.Matrix] = []
    for i in range(n):
        for j in range(i + 1, n):
            e = sp.zeros(n)
            e[i, j] = 1
            e[j, i] = -1
            basis.append(e)
    return tuple(basis)


def primitive_integer_vector(v: sp.Matrix) -> sp.Matrix:
    v = sp.Matrix(v)
    denoms = [int(sp.Rational(x).q) for x in v]
    lcm = 1
    for d in denoms:
        lcm = sp.ilcm(lcm, d)
    ints = [int(sp.Rational(x) * lcm) for x in v]
    nonzero = [abs(x) for x in ints if x]
    if not nonzero:
        return sp.zeros(v.rows, v.cols)
    g = reduce(gcd, nonzero)
    ints = [x // g for x in ints]
    first = next(x for x in ints if x)
    if first < 0:
        ints = [-x for x in ints]
    return sp.Matrix(v.rows, v.cols, ints)


def invariant_skew_space(transfers: Iterable[sp.Matrix]) -> tuple[sp.Matrix, tuple[sp.Matrix, ...]]:
    transfers = tuple(imat(x) for x in transfers)
    if not transfers:
        raise ValueError("at least one transfer is required")
    n = transfers[0].rows
    if any(t.shape != (n, n) or abs(int(t.det())) != 1 for t in transfers):
        raise ValueError("transfers must be unimodular square matrices")
    basis = skew_basis(n)
    rows: list[list[int]] = []
    for t in transfers:
        transformed = [t.T * e * t - e for e in basis]
        for i in range(n):
            for j in range(i + 1, n):
                rows.append([int(m[i, j]) for m in transformed])
    equations = sp.Matrix(rows)
    null = equations.nullspace()
    invariant_forms: list[sp.Matrix] = []
    for vector in null:
        primitive = primitive_integer_vector(vector)
        form = sp.zeros(n)
        for coefficient, e in zip(primitive, basis):
            form += int(coefficient) * e
        invariant_forms.append(imat(form))
    return equations, tuple(invariant_forms)


def form_coordinates(form: sp.Matrix, basis: Sequence[sp.Matrix]) -> sp.Matrix | None:
    form = imat(form)
    if not basis:
        return None
    cols = []
    n = form.rows
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for b in basis:
        cols.append(sp.Matrix([b[i, j] for i, j in pairs]))
    a = sp.Matrix.hstack(*cols)
    y = sp.Matrix([form[i, j] for i, j in pairs])
    solution = sp.linsolve((a, y))
    if solution is sp.EmptySet:
        return None
    values = list(solution)
    if not values:
        return None
    vector = sp.Matrix(values[0])
    if any(x.free_symbols for x in vector):
        return None
    return vector


def factor_rank_one(matrix: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix]:
    """Factor a nonzero integral rank-one matrix as v ell^T over Z."""
    m = imat(matrix)
    if m.rank() != 1:
        raise ValueError("matrix must have rank one")
    j0 = next(j for j in range(m.cols) if any(m[i, j] for i in range(m.rows)))
    col = sp.Matrix(m[:, j0])
    gcol = reduce(gcd, [abs(int(x)) for x in col if x])
    v = col.applyfunc(lambda x: int(x) // gcol)
    i0 = next(i for i in range(v.rows) if v[i, 0])
    ell_entries = []
    for j in range(m.cols):
        q = sp.Rational(m[i0, j], v[i0, 0])
        if q.q != 1 or m[:, j] != int(q) * v:
            raise AssertionError("rank-one factorization is not integral")
        ell_entries.append(int(q))
    ell = sp.Matrix(ell_entries)
    # Move common content from ell into v, leaving ell primitive.
    nonzero = [abs(int(x)) for x in ell if x]
    gell = reduce(gcd, nonzero) if nonzero else 1
    if gell > 1:
        ell = ell.applyfunc(lambda x: int(x) // gell)
        v *= gell
    first = next(int(x) for x in ell if x)
    if first < 0:
        ell = -ell
        v = -v
    if m != v * ell.T:
        raise AssertionError("rank-one reconstruction failed")
    return imat(v), imat(ell)


def proportional_integer(left: sp.Matrix, right_primitive: sp.Matrix) -> tuple[bool, int | None]:
    left = imat(left)
    right = imat(right_primitive)
    if left.shape != right.shape or right.cols != 1:
        raise ValueError("column-vector mismatch")
    if all(x == 0 for x in right):
        raise ValueError("right vector must be nonzero")
    i0 = next(i for i in range(right.rows) if right[i, 0])
    q = sp.Rational(left[i0, 0], right[i0, 0])
    if q.q != 1:
        return False, None
    c = int(q)
    return left == c * right, c


@dataclass(frozen=True)
class WallTransition:
    relative: sp.ImmutableMatrix
    vector: sp.ImmutableMatrix
    covector: sp.ImmutableMatrix
    scalar: int | None
    preserved: bool


def wall_transition(left_transfer: sp.Matrix, right_transfer: sp.Matrix, form: sp.Matrix) -> WallTransition:
    l = imat(left_transfer)
    r = imat(right_transfer)
    s = imat(form)
    relative = imat(l.inv() * r)
    difference = relative - sp.eye(relative.rows)
    if difference == sp.zeros(*difference.shape):
        return WallTransition(
            sp.ImmutableMatrix(relative), sp.ImmutableMatrix(sp.zeros(relative.rows, 1)),
            sp.ImmutableMatrix(sp.zeros(relative.rows, 1)), 0, relative.T * s * relative == s,
        )
    v, ell = factor_rank_one(difference)
    preserved = relative.T * s * relative == s
    proportional, scalar = proportional_integer(s * v, ell)
    if preserved != proportional:
        raise AssertionError("rank-one symplectic criterion failed")
    return WallTransition(
        sp.ImmutableMatrix(relative), sp.ImmutableMatrix(v), sp.ImmutableMatrix(ell), scalar, preserved
    )


def split_difference_matrix(n: int) -> sp.Matrix:
    return sp.Matrix.hstack(-sp.eye(n), sp.eye(n))


def reduced_term_lift(delta: sp.Matrix) -> sp.Matrix:
    delta = imat(delta)
    n, m = delta.shape
    columns = []
    for j in range(m):
        column = sp.zeros(2 * n, 1)
        for i in range(n):
            x = int(delta[i, j])
            column[i, 0] = max(-x, 0)
            column[n + i, 0] = max(x, 0)
        columns.append(column)
    return imat(sp.Matrix.hstack(*columns))


def global_term_form(parameter: sp.Matrix) -> sp.Matrix:
    s = imat(parameter)
    require_skew(s)
    dmap = split_difference_matrix(s.rows)
    return imat(dmap.T * s * dmap)


def globalizable(parameter: sp.Matrix, transfers: Iterable[sp.Matrix]) -> bool:
    s = imat(parameter)
    require_skew(s)
    return all(imat(l).T * s * imat(l) == s for l in transfers)


def local_mutable_form(g: sp.Matrix, parameter: sp.Matrix) -> sp.Matrix:
    g = imat(g)
    s = imat(parameter)
    return imat(g.T * s * g)


def term_restriction(delta: sp.Matrix, parameter: sp.Matrix) -> sp.Matrix:
    delta = imat(delta)
    s = imat(parameter)
    return imat(delta.T * s * delta)


def principal_quantum_form(b0: sp.Matrix, g: sp.Matrix, parameter: sp.Matrix, d: int = 1) -> sp.Matrix:
    """Canonical principal form plus the complete affine S-parameter correction."""
    b0 = imat(b0)
    g = imat(g)
    s = imat(parameter)
    n = b0.rows
    canonical = int(d) * sp.Matrix.vstack(
        sp.Matrix.hstack(sp.zeros(n), -g.T),
        sp.Matrix.hstack(g, b0.T),
    )
    k = sp.Matrix.vstack(g.T, -b0.T)
    return imat(canonical + k * s * k.T)


def quotient_full_principal_form(b0: sp.Matrix, parameter: sp.Matrix, d: int = 1) -> sp.Matrix:
    """The quotient form forced by the initial chart Delta_0=-I, G_0=I."""
    b0 = imat(b0)
    s = imat(parameter)
    n = b0.rows
    c = int(d) * sp.eye(n) + s * b0
    bottom = int(d) * b0.T + b0.T * s * b0
    return imat(sp.Matrix.vstack(sp.Matrix.hstack(s, c), sp.Matrix.hstack(-c.T, bottom)))


def full_principal_global_pullback(
    b0: sp.Matrix, delta: sp.Matrix, parameter: sp.Matrix, d: int = 1
) -> sp.Matrix:
    delta = imat(delta)
    qbar = quotient_full_principal_form(b0, parameter, d)
    n = delta.rows
    j = sp.diag(1, 1)  # overwritten to avoid SymPy's empty-block corner cases
    j = sp.Matrix.vstack(
        sp.Matrix.hstack(delta, sp.zeros(n)),
        sp.Matrix.hstack(sp.zeros(n), sp.eye(n)),
    )
    return imat(j.T * qbar * j)


def full_principal_globalization_conditions(
    b0: sp.Matrix, delta: sp.Matrix, g: sp.Matrix, parameter: sp.Matrix, d: int = 1
) -> tuple[bool, bool]:
    b0 = imat(b0)
    delta = imat(delta)
    g = imat(g)
    s = imat(parameter)
    mutable = delta.T * s * delta == g.T * s * g
    mixed = (delta.T + g.T) * (int(d) * sp.eye(b0.rows) + s * b0) == sp.zeros(b0.rows)
    desired = principal_quantum_form(b0, g, s, d)
    actual = full_principal_global_pullback(b0, delta, s, d)
    if (actual == desired) != (mutable and mixed):
        raise AssertionError("full-principal criterion does not match block equality")
    return mutable, mixed


def matrix_json(m: sp.Matrix) -> list[list[int]]:
    return [[int(m[i, j]) for j in range(m.cols)] for i in range(m.rows)]


__all__ = [
    "PairedAtlas", "PairedState", "TropicalSeed", "WallTransition",
    "chart_key", "enumerate_paired_atlas", "factor_rank_one",
    "form_coordinates", "full_principal_global_pullback",
    "full_principal_globalization_conditions", "global_term_form",
    "globalizable", "imat", "invariant_skew_space", "local_mutable_form",
    "matrix_json", "mutate_exchange", "principal_quantum_form",
    "quotient_full_principal_form", "reduced_term_lift",
    "skew_basis", "split_difference_matrix", "term_restriction",
    "wall_transition",
]
