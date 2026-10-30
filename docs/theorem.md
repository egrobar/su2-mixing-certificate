# Exact color elimination: a SU(2) mixing certificate

For the forest-gauge SU(2) model, integrate one of three active-link colors
exactly. The projected specification on the two remaining colors has
uniform chordal single-link row and column sums bounded by
\[
 \boxed{\alpha_{\rm dec}<0.990737<1
        \quad\text{for }0\le\beta_W\le77/200=0.385.}
\]
At \(0.38\), the stronger margin is
\(\alpha_{\rm dec}<0.963503\).
The bounds are ordinary analytic arguments combined with exact finite
geometry and outward rational arithmetic. They are not Lean theorems or
externally reviewed results. Novelty is unassessed.

This is a structural change from the three-link proof: an entire independent
color class is integrated out, and the resulting interactions are kept
exactly. It remains a **fixed-lattice strong-coupling estimate**.
There is no iteration of a renormalization map, weak-coupling control, or
continuum physical-gap result.

The forest and sharp vMF inputs are adopted from the pinned upstream
manuscript and independently checked in [analytic inputs](analytic-inputs.md).
See [provenance](../PROVENANCE.md) and the full upstream license notice.
The tilt lemmas and derivative enclosure are in
[tilt and interval bounds](tilt-and-intervals.md).

Define C_ij as the supremum of the chordal W1 distance between the retained
single-site conditional laws at i, divided by the input chord at j, over
all pairs of retained configurations differing only at j. Diagonal entries
are zero. The theorem bounds both sup_i sum_j C_ij and sup_j sum_i C_ij.
The input and output chord is Euclidean distance on the unit sphere S^3.
No group-radius or geodesic normalization is implicit.

## Exact finite-range projected specification

The period-two forest fixes exactly one edge of every plaquette.
Its three active edges admit a period-two coloring with labels 0, 1, 2
such that every plaquette contains each color once.
The explicit 16-row table and checks are in
[color geometry](../experiments/su2_color_geometry.py).
Checking all 96 parity/orientation faces proves the infinite-lattice
coloring by periodicity; larger periodic cells independently check counts.

Retain colors 0 and 1, and write \(T\) for color 2.
Conditional on retained links, the variables in \(T\) are independent.
For each \(e\in T\),
\[
 c_e=b\sum_{p\ni e}v_{p,e},\qquad |v_{p,e}|=1,\qquad b=\beta_W,
\]
and its conditional law is \(\nu_{c_e}\).
Integrating these variables produces the finite-range potential
\[
 H_{\rm eff}=\sum_{e\in T}g(c_e),\qquad
 g(c)=\log\int_{S^3}e^{c\cdot x}\,d\sigma(x).
\]
This notation denotes local finite factors, not a convergent infinite
partition integral. Each factor depends on twelve retained links:
two in each of the six incident plaquettes.

The one-site conditional for a retained link \(x\) therefore has exponent
\[
 \sum_{e:\,x\in\operatorname{star}(e)}
       g(c_{0,e}+b L_e x).
\]
There are six terms. In each, \(|c_{0,e}|\le5b\), and \(L_e\) is an
orthogonal quaternion multiplication/inversion map.

Boundary scope matters. For a finite set \(T_0\) of integrated links with
other color-2 links frozen, the exact action is instead
\[
 b\sum_{p:\ e(p)\notin T_0}\tau_p+\sum_{e\in T_0}g(c_e),
\]
where \(e(p)\) is the unique color-2 edge of the included plaquette.
The first sum can retain order-\(b\) interactions. Our certificate concerns
the full projected bulk specification; it does **not** cover those
arbitrary frozen-original-boundary conditionals.

In finite free boxes, every included plaquette's color-2 edge is inside
the box and is integrated. The effective law has only \(g\)-factors, with
partial stars near the boundary. Thus free-box limits can be handled
without discarding an original boundary term.

## Uniform one-site tilt bound

Use the audited global inputs
\[
 W_1(\nu_c,\nu_d)\le|c-d|/4,\quad
 \operatorname{Var}_{\nu_c}f\le\tfrac13\nu_c|\nabla f|^2,\quad
 D^2g\le\tau I,\quad \tau=1/4.
\]
On \(|c|\le6b\), let \(D^2g\ge\lambda_* I\) and \(\|D^3g\|\le M\),
and put \(\Delta=\tau-\lambda_*\).
Taylor expansion of each of the six factors at \(c_{0,e}\) writes the
one-site law as a vMF law tilted by a remainder whose oscillation is at most
\[
 w=3b^2\Delta.
\]
Consequently, for chordal 1-Lipschitz \(f\),
\[
 |\operatorname{Cov}_\kappa(f,h\cdot x)|\le s|h|,
 \qquad
 s=\tfrac14+(e^w-1)e^w\sqrt{\tau/3},
\]
and for a chordal Lipschitz function \(N\),
\[
 |\operatorname{Cov}_\kappa(f,N)|\le(e^w/3)\operatorname{Lip}(N).
\]
These follow from the [tilt/Poincare lemmas](tilt-and-intervals.md),
now applied to six remainder terms.
They bound a single retained-link conditional, not the variance of the
entire retained field.

For finite endpoint changes of another link \(y\), use the straight segment
\(y(t)=(1-t)y_0+ty_1\) in quaternion space. Its norm is at most one.
Before interpolating, write each inverse of a unit quaternion as its
conjugate and retain that linear conjugation off the unit sphere; do not
use the reciprocal of a nonunit quaternion.
The intermediate law need not be a physical boundary law; it is a smooth
positive interpolation with the required endpoints.
All intermediate arguments remain in \(|c|\le6b\).

One detail is essential. The affected map is conformal:
\[
 L_y^{\,T}L_y=|y|^2 I.
\]
Thus \(|L_yx|\) is constant on unit \(x\), even at intermediate \(y(t)\).
The Taylor remainder of that term has oscillation at most
\(b^2|y(t)|^2\Delta/2\), so the bound \(w\) survives.
An operator-norm bound alone would not justify subtracting \(\lambda_*\).

## Two types of factor response

A retained link appears only once in the six staples of a fixed eliminated
link. This follows from the distinct-link plaquette incidence: two distinct
lattice links share at most one plaquette.
For input \(y\), only factors containing both \(x\) and \(y\) matter;
the others contribute constants to the one-site log-density derivative.

Let \(A(r)=I_2(r)/I_1(r)\), so \(|\nabla g(c)|=A(|c|)\).
It is increasing since \(A'(r)=\operatorname{Var}_{\nu_{re_1}}x_1\ge0\).

**Same plaquette.** Suppose \(x,y\) occur in the same staple of factor \(e\).
Then \(c_0\) is independent of \(y\).
Normalize the input direction to unit length and write \(B=\partial_yL_y\);
its operator norm is one. The source is
\[
 b\,\nabla g(c_0+bL_yx)\cdot Bx.
\]
Set \(E(c)=\nabla g(c)-\tau c\). Its Lipschitz constant on the field ball
is at most \(\Delta\). The source splits into:
\[
 b\,\nabla g(c_0)\cdot Bx
 + b^2\tau\,L_yx\cdot Bx
 + b\,[E(c_0+bL_yx)-E(c_0)]\cdot Bx.
\]
The middle term is constant in unit \(x\). Indeed, differentiating
\(|L_yx|^2=|y|^2|x|^2\) gives
\(L_yx\cdot Bx=y\cdot h\).
It therefore vanishes from covariance.

The first term has linear coefficient at most \(b A(5b)\).
The last has chordal Lipschitz constant at most \(2b^2\Delta\):
one factor changes by at most \(b\Delta|x-x'|\), while the other
has norm at most one; its value also has norm at most \(b\Delta\).
Thus the response per such factor occurrence is at most
\[
 u_{\rm same}=s\,b A(5b)+(e^w/3)\,2b^2\Delta.
\]

**Different plaquettes of one eliminated star.**
Here \(L\) is fixed and \(c_0\) changes in a direction \(b v\), \(|v|\le1\).
The source is
\[
 b\,\nabla g(c_0+bLx)\cdot v.
\]
The constant term drops out; its linear part has coefficient at most
\(b^2\tau\).
The gradient Taylor remainder has Lipschitz constant at most \(b^2M\)
as a function of \(x\). After the additional source factor \(b\), the
response is at most
\[
 u_{\rm cross}=s\,b^2\tau+(e^w/3)\,b^3M.
\]

All covariance bounds hold along the full input segment. Integration
therefore bounds the finite-endpoint chordal transport ratio.
No geodesic/chord conversion or small-input assumption is used.

## Count occurrences before distinct neighbors

For each retained input link \(y\), there are six eliminated stars.
In each, one of the other retained links shares \(y\)'s plaquette and
ten lie in the other five plaquettes.
Hence the column sum is at most
\[
 \boxed{\alpha_{\rm dec}
 =s\{6b A(5b)+60b^2\tau\}
 +(e^w/3)\{12b^2\Delta+60b^3M\}.}
\]
The same count bounds every row.
Some pairs occur in multiple stars; their contributions must be added,
rather than treating the distinct-neighbor count as the occurrence count.
The exact geometry checker verifies both.
On the period-six cell it finds 52 distinct retained neighbors per link;
the response estimate correctly uses 66 factor occurrences instead.

This is an ordinary Wasserstein Dobrushin estimate on the projected product
space. It need not be relabeled as the ungauged cube coefficient
\(\kappa_L\), or as the three-link forest block coefficient.

## Rational evaluation

Run:
```sh
python3 -m experiments.su2_decimation_certificate
python3 -m experiments.su2_color_geometry
```
The full response report can be generated at
`.artifacts/su2-decimation-certificate.json`. The standalone release has
19 focused tests; no saved research inputs are required. `python verify.py`
checks the endpoint and geometry together.
It reuses the independently checked integer interval machinery, with
4,096 subintervals covering \(0\le r\le2.4\).
For this larger ball, the safe displayed bounds are
\[
 \lambda_*>0.1393407619,\qquad M<0.0617748750.
\]
The actual program uses the full rational endpoints, not these displays.

The endpoint mean is bounded directly by the positive-series quotient
\(A(r)=2r Z'(r^2)/Z(r^2)\).
The exponential uses a positive rational series and geometric remainder;
\(\sqrt{1/12}\) is bounded by \(289/1000\).
No floating-point value is used in a proof decision.

| Wilson coupling | Certified column-sum upper bound |
| --- | ---: |
| \(1/3\) | \(<0.730171\) |
| \(0.348\) | \(<0.799526\) |
| \(0.380\) | \(<0.963503\) |
| \(0.385\) | \(<0.990737\) |
| \(0.390\) | expression exceeds 1; inconclusive |
| \(0.400\) | expression exceeds 1; inconclusive |

For all \(b\le0.385\), use the same domain constants and the endpoint
bound on \(A(5b)\); monotonicity of the true \(A\), and positivity of the
other factors, extend the endpoint certificate across the interval.
This does not rely on floating-point monotonicity or a grid in \(b\).
The field guard holds since \(6(0.385)=2.31<2.4\).

## Gibbs and physical-observable scope

The projected potential is finite range with smooth positive conditionals
on a compact product space. The row/column contraction gives a unique
projected Gibbs state and exponential local response. Two retained links
sharing an eliminated star have midpoint lattice distance at most two,
so weighted comparison allows \(\lambda<-\frac12\log\alpha_{\rm dec}\).

There is a local correspondence with the full forest-gauge specification.
Given retained links, lift by independently sampling every eliminated
link from \(\nu_{c_e}\).
For a finite retained set, only finitely many eliminated stars meet it.
The projected conditional density contains their \(Z(c_e)\) factors;
conditioning additionally on those eliminated links contributes the
densities \(e^{U_e\cdot c_e}/Z(c_e)\), cancelling the normalizers.
This recovers the original local Wilson weights. Conversely, integrating
those same neighboring eliminated variables in the original finite
conditional gives the projected factors. This establishes the local
specification correspondence; it does not invoke a nonlocal gauge
transformation of an infinite configuration.

Finite free boxes integrate all their color-2 variables and yield partial
stars at the boundary. Deleting plaquettes can only reduce the number of
one-site factors (at most six), the residual staple norm (at most \(5b\)),
and the two occurrence counts (at most six and sixty). The same local
majorant therefore covers these partial-star free-box laws uniformly.
Every fixed interior conditional eventually has
the full projected factors. Compactness and projected uniqueness identify
their limit. The finite forest gauge identity then gives the same
gauge-invariant free-box expectations for the ungauged Wilson theory.

A local observable involving eliminated links can be conditionally averaged
over the independent vMF variables. Its support enlarges by only one
plaquette neighborhood, and the quarter transport bound gives finite
Lipschitz constants for that average.
For sufficiently separated observables, their eliminated-link sets are
conditionally independent, so their covariance equals the covariance of
the retained averages. Thus the exponential bound transfers to the local
gauge-invariant Lipschitz observables of the free-box limiting state.

This is not a statement about every originally frozen boundary and does
not establish uniqueness of all ungauged Gibbs states.
A physical transfer-matrix spectral statement requires reflection-positive
reconstruction and a dense physical observable family. The
[transfer corollary](transfer-corollary.md) supplies those analytic
steps for this selected free-box limit, conditional on the certificate
and lift above. No Lean formalization of those steps is supplied.
Nothing here provides a positive rate uniform along the continuum trajectory.

The boundary and observable restrictions also prevent importing an
arbitrary-original-boundary area-law theorem without a separate adapter.
No Wilson-loop area-law claim is part of this release.

The proof-to-checker boundary and the finite/infinite-volume comparison
are detailed in [verification scope](verification-scope.md).

## Group dependence and a negative control

The coloring and the conformal norm cancellation also exist for U(1).
The quantitative S^3 transport and variance inputs do not transfer to
the circle. At circle Haar measure a coordinate has variance 1/2,
and testing transport response with that coordinate gives sensitivity
at least 1/2. Thus a circle version of this same nonnegative majorant,
with uniform constants tau and s, must have tau>=1/2 and s>=1/2.
Its cross-factor term alone is at least

    60 beta_W^2 s tau >= 15 beta_W^2.

Such a certificate is necessarily inconclusive for beta_W>=1/sqrt(15).
This is a lower bound on the sufficient expression, not on the true
U(1) influence matrix or a proof of a phase transition. It shows why
using the SU(2) quarter and Poincare constants in an abelian test would
be invalid; this release supplies no weak-coupling mechanism for either
group.
