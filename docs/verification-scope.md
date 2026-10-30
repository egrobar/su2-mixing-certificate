# Verification scope and the Gibbs comparison

The certificate combines a written analytic proof with finite exact
arithmetic and geometry. Passing the verifier checks the latter; it
does not turn the analytic proof into a kernel-checked formal theorem.

## What is checked computationally

`python verify.py` recomputes the positive-series derivative enclosures
on 4096 closed cells covering 0 <= r <= 12/5. It checks the rational
endpoint at beta_W=77/200, the field-domain guard, the square-root
majorant, the 6+60 occurrence formula, both row and column counts,
all 96 parity face types, and the interaction range using unwrapped
integer coordinates. Its compact outward summary gives

    alpha <= 247684156026907/250000000000000
          = 0.990736624107628 < 0.990737 < 1.

The internal exact rational endpoint is slightly smaller. The verifier
recomputes it rather than trusting this summary or a saved JSON file.
Floating-point fields are display-only. Runtime information is recorded
separately from proof parameters; Python patch versions and test timings
are not mathematical inputs.

The 19 tests additionally cover interval signs and division, positive
series against longer independent sums, rejection of invalid domains,
the exact success/failure controls, quaternion conformality and geometry
corruption. These are software checks, not a formalization of calculus,
optimal coupling, infinite-volume limits or reflection positivity.

## Forest gauge and finite boxes

For the tabulated direction field sigma, the 16 parity cases obey
sigma(x+e_sigma(x))=sigma(x-e_sigma(x))=sigma(x). Consequently the fixed
edges extend in straight lines. If a fixed incoming edge from a different
direction entered x, the same identity applied at its starting vertex
would force sigma(x) to be that other direction, a contradiction. The
fixed graph therefore consists of disjoint straight lines. Intersecting
it with any finite rectangular free box gives a forest.

Choose one root in each finite component and recursively use vertex
gauge transformations to set the forest links to identity. Haar
invariance makes the remaining link variables product Haar; the
normalized gauge variables integrate to one. The Wilson action and
gauge-invariant observables are unchanged. This proves the finite-box
gauge identity used by the theorem. It does not gauge-fix an infinite
configuration by a nonlocal transformation.

Periodic geometry checks validate local counts, not acyclicity of the
wrapped graph. This is why the physical statement uses free rectangles.

## A finite-volume comparison proof

Let delta_i(f) be the single-coordinate Lipschitz seminorm of f for the
unit-S^3 chord. On a finite retained product space, suppose the reference
law mu has influence matrix C with row sums at most alpha<1. Let nu be
a second law, and b_i a uniform upper bound on the chordal W1 discrepancy
between its i-conditional and that of mu at the same outside coordinates.

Choose a coupling of mu and nu minimizing the sum of expected coordinate
chords; compactness ensures a minimizer. Write d_i for these expected
chords. Resampling the i-coordinate pair by an optimal coupling of its
two conditional laws preserves both full marginals. Optimality, the
triangle inequality and the influence bound give

    d_i <= sum_j C_ij d_j + b_i.

Measurable optimal coupling selection on compact spaces justifies this
resampling. Iterating the inequality and using the row contraction gives

    d <= D b,       D=(I-C)^(-1)=sum_{n>=0} C^n,
    |mu(f)-nu(f)| <= sum_{i,j} delta_i(f) D_ij b_j.       (1)

Only the reference influence matrix is used. The perturbed law need
not obey its own contraction bound.

## Covariance decay, with an explicit bound

Take nu proportional to exp(tg) mu for a real local Lipschitz g. For a
chordal 1-Lipschitz test h on one coordinate, osc(h)<=2 and
osc_i(g)<=2 delta_i(g). Under every intermediate one-coordinate tilt,
the elementary range covariance inequality gives

    |Cov(h,g)| <= osc(h) osc_i(g)/4 <= delta_i(g).

Integrating the tilt parameter and using Kantorovich duality therefore
allows b_i<=|t|delta_i(g) in (1). Divide by |t| and let t tend to zero:

    |Cov_mu(f,g)| <= sum_{i,j} delta_i(f) D_ij delta_j(g). (2)

All derivatives exist for finite-volume bounded functions. Complex
observables follow by real and imaginary parts, with a harmless change
in prefactor.

The matrix is supported on pairs at midpoint lattice l1 distance at
most two. For 0<lambda<-log(alpha_0)/2, alpha_0=990737/1000000, set
q=alpha_0 exp(2lambda)<1. The triangle inequality and positivity yield

    sup_i sum_j exp(lambda d(i,j)) D_ij <= 1/(1-q).

If the supports of f and g are at distance at least r, (2) gives

    |Cov_mu(f,g)| <= exp(-lambda r)
                    ||delta f||_1 ||delta g||_infinity /(1-q). (3)

The transpose estimate follows from the identical column bound. The
same estimates apply to the partial-star free-box laws: no omitted
plaquette can increase the field norms or occurrence counts used by the
majorant. The local constants and exponent are uniform in box size.

## Infinite volume and the physical lift

Compactness gives subsequential local limits. Every fixed interior
conditional eventually has the full projected local factors. The
comparison argument applied to larger finite sets shows that the effect
of different exteriors on a fixed local observable tends to zero: each
influence path from its support to the exterior has length tending to
infinity, and its total weight is bounded by a geometric tail. This
gives uniqueness of the projected Gibbs state and identifies the
free-box limit. Equation (3) passes to that limit.

Given retained links, the eliminated links are independent vMF variables.
Conditionally integrating a local observable enlarges its support by
only one plaquette neighborhood. Varying a retained coordinate changes
only finitely many vMF parameters; the quarter transport bound gives
a finite Lipschitz constant for the averaged observable. For sufficiently
separated supports, the two eliminated sets are conditionally independent,
so their covariance is exactly that of their retained averages.

Combine this local lift with the finite gauge identity above to obtain
the stated covariance decay for gauge-invariant local observables in
the original Wilson free-box limiting state. This does not require,
or assert, uniqueness for every original frozen-link boundary.

Translation covariance of original finite free-box expectations, together
with the unique projected limit for every translated exhaustion, makes
the original physical limiting state translation invariant. Symmetric
exhaustions also give reflection invariance. The period-two proof gauge
does not introduce a physical two-step time unit.

## Remaining mathematical trust boundary

The written proof uses compact-space optimal transport, the adopted vMF
lemmas, differentiation under finite integrals, and standard compactness
and spectral arguments. These have been checked in project audits but
not externally refereed or formalized in Lean. The separate
[transfer note](transfer-corollary.md) also requires the original Wilson
reflection-positive reconstruction and its dense physical observable
algebra. The arithmetic script does not prove those steps.

The package makes no continuum construction assumption. Its endpoint
is a fixed-lattice strong-coupling statement. The methods and numerical
threshold do not supply a bound along the continuum limit.
