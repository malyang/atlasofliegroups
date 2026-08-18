#!/usr/bin/env python3
"""Exact symbolic identities behind the whole-Hall type-A quotient.

The module keeps the two ingredients separate:

1. the one-dimensional-extension calculation in the automorphism-normalized
   Ringel--Hall basis, followed by the self-Euler normalization and a skew
   bicharacter twist;
2. elimination between the two ordered Berenstein--Zelevinsky exchange
   products in the canonical based quantum torus normalization.

All calculations are exact SymPy expressions. The code deliberately does not
pretend to compute Hall numbers for an arbitrary bound quiver.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Dict, Sequence, Tuple

import sympy as sp

v = sp.Symbol("v", nonzero=True)
q = v**2


@dataclass(frozen=True)
class HallOneExtIdentity:
    hom_xy: int
    hom_yx: int
    exact_skew: int
    target_skew: int
    forward_split_exponent: int
    reverse_split_exponent: int
    commutator_exponent: int


def verify_one_ext_identity(h: int, hp: int, lam: int) -> HallOneExtIdentity:
    """Verify the normalized Hall identity when Ext^1(X,Y)=F_q and Ext^1(Y,X)=0."""
    if h < 0 or hp < 0:
        raise ValueError("Hom dimensions must be nonnegative")
    mu = hp - h + 1
    untwisted_forward = hp - h - 1
    untwisted_reverse = h - hp - 1
    assert untwisted_forward == mu - 2
    assert untwisted_reverse == -mu
    omega = lam - mu
    forward = untwisted_forward + omega
    reverse = untwisted_reverse - omega
    assert forward == lam - 2
    assert reverse == -lam
    split, middle = sp.symbols("split middle")
    lhs_forward = v**forward * (split + (v**2 - 1) * middle)
    lhs_reverse = v**reverse * split
    comm_exp = 2 * lam - 2
    residual = sp.expand(lhs_forward - v**comm_exp * lhs_reverse)
    expected = sp.expand((v**lam - v**(lam - 2)) * middle)
    assert sp.simplify(residual - expected) == 0
    return HallOneExtIdentity(h, hp, mu, lam, forward, reverse, comm_exp)


def verify_bz_exchange_elimination(lam: int) -> bool:
    """Verify elimination of one monomial from two ordered BZ products."""
    mp, mm = sp.symbols("Mplus Mminus")
    forward = v**lam * mp + v**(lam - 1) * mm
    reverse = v**(-lam) * mp + v**(1 - lam) * mm
    residual = sp.expand(forward - v**(2 * lam - 2) * reverse)
    expected = sp.expand((v**lam - v**(lam - 2)) * mp)
    return sp.simplify(residual - expected) == 0


def verify_a4_ordered_exchange() -> Dict[str, str]:
    x2, x3 = sp.symbols("X2 X3")
    forward = x2 + v**-1 * x3
    reverse = x2 + v * x3
    residual = sp.expand(forward - v**-2 * reverse)
    expected = sp.expand((1 - v**-2) * x2)
    assert sp.simplify(residual - expected) == 0
    return {
        "forward": str(forward),
        "reverse": str(reverse),
        "q_commutator": str(residual),
        "isolated_exact_smoothing": str(x2),
    }

Arc = Tuple[int, int]


def norm_arc(a: int, b: int) -> Arc:
    if a == b:
        raise ValueError("arc endpoints must be distinct")
    return (a, b) if a < b else (b, a)


def boundary(arc: Arc, m: int) -> bool:
    a, b = norm_arc(*arc)
    return b == a + 1 or (a == 0 and b == m - 1)


def diagonals(m: int) -> Tuple[Arc, ...]:
    return tuple((a, b) for a in range(m) for b in range(a + 1, m)
                 if not boundary((a, b), m))


def crossing(a: Arc, b: Arc) -> bool:
    a = norm_arc(*a)
    b = norm_arc(*b)
    if len(set(a + b)) != 4:
        return False
    e0, e1, e2, e3 = sorted(a + b)
    return {a, b} == {(e0, e2), (e1, e3)}


def compatible(arcs: Sequence[Arc]) -> bool:
    return all(not crossing(a, b) for a, b in combinations(arcs, 2))


def triangulations(m: int) -> Tuple[Tuple[Arc, ...], ...]:
    ds = diagonals(m)
    return tuple(tuple(c) for c in combinations(ds, m - 3) if compatible(c))


def exchange_pair_witness(a: Arc, b: Arc, m: int,
                          tris: Sequence[Tuple[Arc, ...]]) -> bool:
    if not crossing(a, b):
        return False
    for tri in tris:
        if a not in tri:
            continue
        replaced = tuple(sorted((set(tri) - {a}) | {b}))
        if len(replaced) == m - 3 and compatible(replaced):
            return True
    return False


def common_cluster_witness(a: Arc, b: Arc,
                           tris: Sequence[Tuple[Arc, ...]]) -> bool:
    if crossing(a, b):
        return False
    return any(a in tri and b in tri for tri in tris)


def catalan(n: int) -> int:
    return comb(2*n, n) // (n+1)


def audit_type_a(rank: int) -> Dict[str, int]:
    m = rank + 3
    ds = diagonals(m)
    tris = triangulations(m)
    cross = [(a, b) for a, b in combinations(ds, 2) if crossing(a, b)]
    noncross = [(a, b) for a, b in combinations(ds, 2) if not crossing(a, b)]
    assert len(ds) == m * (m - 3) // 2
    assert len(cross) == comb(m, 4)
    assert len(tris) == catalan(m - 2)
    assert all(exchange_pair_witness(a, b, m, tris) for a, b in cross)
    assert all(common_cluster_witness(a, b, tris) for a, b in noncross)
    local_hall_checks = 0
    for h in range(4):
        for hp in range(4):
            for lam in range(-3, 4):
                verify_one_ext_identity(h, hp, lam)
                assert verify_bz_exchange_elimination(lam)
                local_hall_checks += 1
    return {
        "rank": rank,
        "vertices": m,
        "diagonals": len(ds),
        "crossing_pairs": len(cross),
        "compatible_pairs": len(noncross),
        "triangulations": len(tris),
        "exchange_pair_witnesses": len(cross),
        "common_cluster_witnesses": len(noncross),
        "symbolic_local_identity_checks": local_hall_checks,
    }
