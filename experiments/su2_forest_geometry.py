#!/usr/bin/env python3
"""Finite-period geometry audit for the pinned 4Y-M forest gauge.

Source: creelie/4Y-M, strong-coupling/forest/main.tex,
commit b2d359f1197f9f50aae8f6890edb6d9106531fa3:
https://github.com/creelie/4Y-M/blob/b2d359f1197f9f50aae8f6890edb6d9106531fa3/strong-coupling/forest/main.tex

This checks the tabulated period-two line field and finite periodic cells.
It does not assert periodic line components are forests or supply an
infinite-volume Gibbs adapter.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path


SIGMA_ROWS = {
    (0, 0, 0, 0): (0, 0, 2, 3, 1),
    (1, 0, 0, 0): (0, 0, 3, 1, 2),
    (0, 1, 0, 0): (2, 3, 0, 2, 1),
    (1, 1, 0, 0): (3, 2, 0, 1, 3),
    (0, 0, 1, 0): (3, 1, 2, 0, 3),
    (1, 0, 1, 0): (1, 3, 1, 0, 2),
    (0, 1, 1, 0): (2, 1, 3, 2, 0),
    (1, 1, 1, 0): (1, 2, 1, 3, 0),
    (0, 0, 0, 1): (1, 2, 1, 3, 0),
    (1, 0, 0, 1): (2, 1, 3, 2, 0),
    (0, 1, 0, 1): (1, 3, 1, 0, 2),
    (1, 1, 0, 1): (3, 1, 2, 0, 3),
    (0, 0, 1, 1): (3, 2, 0, 1, 3),
    (1, 0, 1, 1): (2, 3, 0, 2, 1),
    (0, 1, 1, 1): (0, 0, 3, 1, 2),
    (1, 1, 1, 1): (0, 0, 2, 3, 1),
}

BLOCK_BASES = {
    (0, 1): {(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 0, 1, 1),
             (1, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 1, 1)},
    (2, 3): {(0, 0, 0, 0), (0, 0, 1, 1), (0, 1, 0, 0), (0, 1, 1, 1),
             (1, 0, 0, 0), (1, 0, 1, 1), (1, 1, 0, 0), (1, 1, 1, 1)},
}


def add(x, mu, step=1, p=None):
    if p is None:
        return tuple(a + (step if i == mu else 0) for i, a in enumerate(x))
    return tuple((a + (step if i == mu else 0)) % p for i, a in enumerate(x))


def sig(x):
    return SIGMA_ROWS[tuple(v % 2 for v in x)][0]


def plaquette_edges(x, mu, nu, p):
    return ((x, mu), (add(x, mu, p=p), nu),
            (add(x, nu, p=p), mu), (x, nu))


def all_plaquettes(p):
    return [(x, mu, nu) for x in product(range(p), repeat=4)
            for mu in range(4) for nu in range(mu + 1, 4)]


def is_block_plaquette(x, mu, nu):
    r = tuple(v % 2 for v in x)
    return r in BLOCK_BASES.get((mu, nu), set())


def check_table_parity():
    """Check the 16-row table's local proof obligations once, exactly."""
    for parity, row in SIGMA_ROWS.items():
        s, *pi = row
        assert len(row) == 5
        assert pi[s] == s
        other = set(range(4)) - {s}
        assert {pi[mu] for mu in other} == other
        assert all(pi[mu] != mu for mu in other)
        for mu in range(4):
            for nu in range(mu + 1, 4):
                incidence = int(s in (mu, nu)) + int(pi[mu] == nu) + int(pi[nu] == mu)
                assert incidence == 1, (parity, mu, nu, row, incidence)
    return {"parity_rows": len(SIGMA_ROWS), "plaquette_incidence_cases": 16 * 6,
            "exactly_one_fixed_edge_per_case": True}


def check_finite_free_boxes():
    """Run finite path-segment checks once; they do not depend on torus size."""
    box_cases = 0
    for lengths in product(range(1, 5), repeat=4):
        for parity_origin in product(range(2), repeat=4):
            box_cases += 1
            verts = {tuple(parity_origin[d] + i for d, i in enumerate(offset))
                     for offset in product(*(range(n) for n in lengths))}
            fedges = {(x, sig(x)) for x in verts if add(x, sig(x)) in verts}
            adj = defaultdict(set)
            for x, mu in fedges:
                y = add(x, mu)
                adj[x].add(y)
                adj[y].add(x)
            unseen = set(adj)
            while unseen:
                root = next(iter(unseen))
                stack = [root]
                comp = set()
                while stack:
                    v = stack.pop()
                    if v in comp:
                        continue
                    comp.add(v)
                    unseen.discard(v)
                    stack.extend(adj[v] - comp)
                comp_edges = sum(len(adj[v] & comp) for v in comp) // 2
                assert comp_edges == len(comp) - 1
                assert max(map(len, (adj[v] for v in comp))) <= 2
                dirs = {mu for x, mu in fedges if x in comp and add(x, mu) in comp}
                assert len(dirs) == 1
    return {"cases": box_cases, "side_lengths_per_axis": [1, 2, 3, 4],
            "origins_per_residue": 16, "all_components_are_path_segments": True}


def check_period(p):
    assert p >= 4 and p % 2 == 0
    sites = list(product(range(p), repeat=4))
    links = {(x, mu) for x in sites for mu in range(4)}
    fixed = {(x, mu) for x in sites if (mu := sig(x)) >= 0}
    active = links - fixed
    plaquettes = all_plaquettes(p)

    block_faces = set()
    owner = defaultdict(list)
    active_per_block = {}
    incident_faces = defaultdict(list)
    for x, mu, nu in plaquettes:
        edges = plaquette_edges(x, mu, nu, p)
        face = (x, mu, nu)
        for edge in edges:
            incident_faces[edge].append(face)
        fixed_edges = [e for e in edges if e in fixed]
        assert len(fixed_edges) == 1, (p, x, mu, nu, fixed_edges)
        if is_block_plaquette(x, mu, nu):
            block_faces.add((x, mu, nu))
            active_edges = [e for e in edges if e not in fixed]
            assert len(active_edges) == 3
            active_per_block[(x, mu, nu)] = active_edges
            for edge in active_edges:
                owner[edge].append((x, mu, nu))

    assert len(fixed) == p**4
    assert len(active) == 3 * p**4
    assert len(block_faces) == p**4
    assert set(owner) == active
    assert all(len(v) == 1 for v in owner.values())

    nonblock_active_owner_counts = Counter()
    outside_faces_per_active = Counter()
    for x, mu, nu in plaquettes:
        face = (x, mu, nu)
        if face in block_faces:
            continue
        active_edges = [e for e in plaquette_edges(x, mu, nu, p) if e in active]
        assert len(active_edges) == 3
        block_owners = [owner[e][0] for e in active_edges]
        assert len(set(block_owners)) == 3, (p, face, active_edges, block_owners)
        nonblock_active_owner_counts[face] = len(set(block_owners))
        for e in active_edges:
            outside_faces_per_active[e] += 1

    assert len(nonblock_active_owner_counts) == len(plaquettes) - len(block_faces)
    assert set(outside_faces_per_active.values()) == {5}

    total_degree = {}
    for e in active:
        block = owner[e][0]
        block_partners = set(active_per_block[block]) - {e}
        ext_partners = set()
        incident = []
        for face in incident_faces[e]:
            edges = plaquette_edges(face[0], face[1], face[2], p)
            incident.append(face)
            if face not in block_faces:
                ext_partners.update({q for q in edges if q in active and q != e})
        assert len(incident) == 6
        assert sum(face in block_faces for face in incident) == 1
        assert len(block_partners) == 2
        assert len(ext_partners) == 10, (p, e, len(ext_partners))
        assert not (block_partners & ext_partners)
        total_degree[e] = len(block_partners | ext_partners)
        assert total_degree[e] == 12

    # Period-two direction-field identities and periodic line structure.
    for x in sites:
        s = sig(x)
        pi = SIGMA_ROWS[tuple(v % 2 for v in x)][1:]
        assert all(sig(add(x, mu, p=p)) == pi[mu] for mu in range(4))
        assert pi[s] == s
        rest = [mu for mu in range(4) if mu != s]
        assert sorted(pi[mu] for mu in rest) == sorted(rest)
        assert len({(mu, pi[mu]) for mu in rest}) == 3
        assert all(pi[mu] != mu for mu in rest)
        assert sig(add(x, s, p=p)) == s
        assert sig(add(x, s, step=-1, p=p)) == s

    return {
        "side": p,
        "vertices": len(sites),
        "positive_links": len(links),
        "fixed_links": len(fixed),
        "active_links": len(active),
        "plaquettes": len(plaquettes),
        "block_plaquettes": len(block_faces),
        "nonblock_plaquettes": len(plaquettes) - len(block_faces),
        "fixed_links_per_plaquette": [1],
        "active_link_block_owners": sorted(set(len(v) for v in owner.values())),
        "nonblock_active_links": sorted(set(nonblock_active_owner_counts.values())),
        "nonblock_distinct_block_owners": sorted(set(nonblock_active_owner_counts.values())),
        "outside_plaquettes_per_active_link": sorted(set(outside_faces_per_active.values())),
        "block_partners_per_active_link": [2],
        "external_partners_per_active_link": [10],
        "total_active_partners_per_link": sorted(set(total_degree.values())),
    }


def main():
    table = check_table_parity()
    free_boxes = check_finite_free_boxes()
    reports = [check_period(p) for p in (4, 6, 8)]
    out = {"source": "creelie/4Y-M strong-coupling/forest/main.tex at b2d359f1197f9f50aae8f6890edb6d9106531fa3",
           "scope": "exact finite periodic geometry audit; no global Gibbs claim",
           "table_parity_proof_checks": table,
           "finite_free_box_check": free_boxes,
           "reports": reports}
    rendered = json.dumps(out, indent=2, sort_keys=True) + "\n"
    output = Path(".artifacts/su2-forest-geometry.json")
    output.parent.mkdir(exist_ok=True)
    output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
