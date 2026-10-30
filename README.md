# SU(2) mixing by exact color elimination

A reproducible analytic and rational-arithmetic certificate for
four-dimensional SU(2) Wilson lattice gauge theory at
**0 <= beta_W <= 77/200 = 0.385**.

After a specified forest gauge and exact integration of one independent
link color, the retained singleton chordal Wasserstein influence matrix
has both row and column sums below **0.990737 < 1**. The computed rational
upper endpoint is approximately **0.9907366241076274**. The exact fraction,
not this decimal, decides success.

This release contains a written analytic argument and a reproducible
computer-assisted bound. It is not a Lean proof, peer-reviewed result,
claim of priority, or solution of the continuum Yang--Mills mass-gap
problem. The adopted upstream lemmas and original contributions are
identified in [PROVENANCE.md](PROVENANCE.md).

## The statement and its scope

Use the Wilson weight exp(beta_W sum_p Re Tr(U_p)/2) and normalized Haar
measure. The metric is Euclidean chordal distance on unit S^3. The
certificate is uniform over the retained-link boundary configurations
of the exactly projected specification. It implies a unique projected
Gibbs state and exponential local Lipschitz covariance decay.

The same estimate covers partial stars in finite free boxes. Exact
conditional lifting and the finite forest gauge identity transfer decay
to local gauge-invariant Lipschitz observables in the original Wilson
free-box limiting state. It does not certify arbitrary frozen original
boundaries: those can retain bare plaquette terms after elimination.

The separate [physical transfer corollary](docs/transfer-corollary.md)
derives a_t Delta E > 0.0046315 for the selected fixed-spacing physical
vacuum, conditional on the analytic mixing theorem and lift. The
arithmetic verifier does not formalize that reconstruction argument.
No uniform continuum limit, area law, full-measure LSI, or periodic/twist
sector theorem is asserted here.

## Reproduce

Python 3.11 or newer, standard library only; no pip installation or
downloaded inputs are needed. From this directory:

```sh
python verify.py
python -m unittest discover -s tests -v
python verification/independent_check.py
```

Use `python3` if that is your interpreter name. The verifier recomputes
the derivative enclosures and geometry from source. The 19 focused tests
exercise signed interval operations, positive-series tails, exact
quaternion identities, geometry and success/failure controls. Saved
verification results are conveniences, not trusted inputs.
The independent checker imports none of the certificate implementation;
it proves contraction with a coarser rational bound below 0.997328724.
The sharper 0.990737 bound is verified by `verify.py`.

For the larger raw arithmetic report:

```sh
python -m experiments.su2_decimation_certificate --output .artifacts/decimation.json
python -m experiments.su2_color_geometry --output .artifacts/geometry.json
```

The success at 0.385 and failure of this sufficient expression at 0.390
are both checked. Failure of an upper bound does not prove failure of
mixing. The interval 0 <= beta_W <= 0.385 is covered by analytic
monotonicity with fixed field-domain constants, not by sampling beta.

## Read the argument

1. [Theorem and exact elimination](docs/theorem.md).
2. [Adopted sharp one-link lemmas](docs/analytic-inputs.md).
3. [Tilt estimates and rigorous interval bounds](docs/tilt-and-intervals.md).
4. [Verification scope and Gibbs comparison](docs/verification-scope.md).
5. [Physical transfer corollary](docs/transfer-corollary.md).

The geometry checker uses periodic cells to verify local incidences.
Those cells do not justify gauge fixing wrapped lines on a torus. The
physical lift here uses finite free rectangles, where the fixed graph
is a forest.

## Reuse and review

See [provenance](PROVENANCE.md) and [third-party notices](THIRD_PARTY_NOTICES.md)
before reusing adapted material. A separate license for original additions
has not yet been selected. Mathematical corrections should identify the
specific lemma, scope assumption, or reproducible failing check.
