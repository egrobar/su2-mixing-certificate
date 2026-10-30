# Fixed-spacing physical transfer corollary

This is a separate analytic implication of the mixing theorem, not a
result checked by the arithmetic program. The state, physical observable
algebra, boundary exhaustion and transfer normalization are specified
below. It concerns the fixed-spacing Wilson theory and does not rely on
an announced continuum construction or a Kogut--Susskind limiting step.

## 1. Conditional theorem, with normalization and state specified

Use four-dimensional SU(2) Wilson weights

\[
 \exp\!\left(\beta_W\sum_p{\operatorname{ReTr}U_p\over2}\right),
 \qquad 0<\beta_W\le {77\over200}=0.385.
\]

Assume the analytic influence estimate and physical free-box lift in
[the mixing theorem](theorem.md). In particular:

* the projected chordal W1 matrix has row and column sums at most
  alpha_0=990737/1000000;
* an influence step spans at most two midpoint lattice-distance units;
* the same majorant controls partial stars in finite free rectangles;
* conditional lifting and the exact finite forest-gauge identity
  transfer local gauge-invariant Lipschitz covariance bounds to the
  original Wilson free-box limiting state omega.

Let T be the normalized positive physical transfer operator in the
Osterwalder–Schrader (OS) reconstruction of this selected state.
The vacuum is Omega and T Omega=Omega. Then

\[
 \boxed{\quad
 \|T|_{\Omega^\perp}\|\le\sqrt{\alpha_0},\qquad
 a_t\Delta E\ge -{1\over2}\log\alpha_0
 >{9263\over2000000}=0.0046315 .
 \quad}
 \tag{1}
\]

Here a_t is one Euclidean time step in physical units and
H_T=-a_t^(-1)log T is the transfer Hamiltonian relative to the vacuum.
The logarithmic rate is approximately 0.004653084185 per lattice step;
the rational number in (1) is a conservative lower bound.

This is an analytic corollary conditional on the mixing theorem and lift for the **infinite-volume
free-box physical vacuum at fixed lattice spacing**. It is not
a new coupling threshold, an identification of the actual glueball
mass, an externally reviewed theorem, or a Lean proof.
The beta_W=0 case is the degenerate product limit and is not needed
for the displayed positive-coupling statement.

## 2. The existing covariance estimate supplies a common exponent

Write

\[
 r=-{1\over2}\log\alpha_0.
\]

The range-two influence bound supports weighted comparison for every
0<lambda<r, since alpha_0 exp(2 lambda)<1. Thus the projected field
has exponential local Lipschitz covariance bounds at every such rate.

For an observable involving eliminated links, conditionally integrate
those links using their independent vMF laws. Its retained average
has a finite Lipschitz norm and support enlarged by only one plaquette
neighborhood. For sufficiently separated observables the eliminated
sets are conditionally independent, so covariance equals covariance
of the retained averages. The exact finite forest-gauge identity
then identifies the corresponding original physical expectation.

In particular, if F is any fixed smooth local gauge-invariant
positive-time cylinder observable, its reflected support and its
translate by n time steps are separated by at least n-c_F.
For all sufficiently large n,

\[
 |\operatorname{Cov}_\omega(\Theta F,\tau_nF)|
 \le K_{F,\lambda}e^{-\lambda n}.
 \tag{2}
\]

The support offset is absorbed into K. This prefactor may depend
on F and lambda, but the exponent does not depend on F or the
spatial or temporal box extents. The period-two proof gauge does not
change that uniform exponent: the finite local averaging bounds hold
for every retained configuration and every translated support.

## 3. Translation invariance and free-boundary reflection positivity

Keep the period-two forest and coloring fixed globally. Every
translated rectangular exhaustion has the same projected bulk
specification in its interior. Compactness and projected uniqueness
identify its limit; the local conditional lift identifies the
forest-gauge limit. The exact finite gauge identity then gives a
common limit for the original gauge-invariant expectations.

The original finite-volume Wilson laws satisfy exact translation
covariance:
omega_Lambda(tau_a F)=omega_(Lambda-a)(F).
Taking the two exhaustions proves omega(tau_a F)=omega(F).
This does not require the proof forest to be invariant under unit
translation. Symmetric exhaustions also give reflection invariance.
No uniqueness claim for arbitrary original frozen-link boundaries
is used.

Reflection positivity can be checked on finite free cylinders
directly, so periodic spatial boundary conditions need not be imported.
For a finite free spatial graph Gamma, let

\[
 M(U)=\exp\!\left({\beta_W\over2}\sum_{p\subset\Gamma}q(U_p)\right),
 \quad K(U,V)=\prod_{e\subset\Gamma}k_{\beta_W}(U_eV_e^{-1}),
 \quad T_\Gamma=MKM,
\]

where q=Re Tr/2 and k_beta is the normalized Wilson one-link density.
The nonnegative character coefficients of this density make K a
positive operator; hence T_Gamma is positive on any such finite graph.
It commutes with the spatial gauge projector G.

After fixing the temporal-link forest in an open cylinder with
slices -N,...,N, its partition function is proportional to

\[
 \langle M,T_\Gamma^{2N}M\rangle.
\]

Each interior slice has the full spatial factor M^2. The boundary
vectors supply the missing half-weight at both free ends.

Integrate an entire positive-time gauge-invariant F over its half
of this cylinder. The resulting boundary vector v_F is invariant
under the common residual spatial gauge transformations. A cut at
a site reflection gives a positive multiple of ||v_F||^2; a cut
between slices gives a positive multiple of
<v_F,T_Gamma v_F>. Both are nonnegative.

For a multi-time observable, its individual time-slice factors need
not be gauge invariant. Use the unprojected positive kernel while
integrating the whole half-path; restrict to physical vectors at
the cut. Inserting separate gauge projectors inside the support
of such an observable is not justified by this argument.

These site and link reflection inequalities pass to the local weak
limit because their observables are bounded cylinders. They give
the usual positive physical Wilson transfer reconstruction.
The operator construction and Gauss-law sector are classical;
the free-end boundary-vector check above specifies our exhaustion.
[Lüscher, positive Wilson transfer construction](https://bib-pubdb1.desy.de/record/396349/files/7611148.pdf?version=1).

The forest used to prove the covariance estimate is not used to
assert reflection positivity of a gauge-fixed action. Reflection
positivity belongs to the original Wilson state.

## 4. Density and the spectral argument

Define the physical OS space by completing the positive-time
gauge-invariant local cylinder algebra modulo its null space, with

\[
 \langle[F],[G]\rangle_{\rm OS}=\omega((\Theta F)G).
\]

Here Theta includes reflection and complex conjugation. On each
finite set of links, Peter–Weyl polynomials uniformly approximate
continuous functions. Gauge averaging preserves locality and gives
smooth invariant approximants. Since
||[F]||_OS<=||F||_infinity, these approximations also converge in
OS norm. Centering by omega(F) yields a dense family in Omega-perp.
This uses all smooth physical positive-time cylinders, not a finite
collection of Wilson loops or one selected spectral channel.

For centered F and the normalized positive contraction T,

\[
 0\le\langle[F],T^n[F]\rangle_{\rm OS}
 =\omega((\Theta F)\tau_nF)
 \le K_{F,\lambda}e^{-\lambda n}.
 \tag{3}
\]

Let nu_F be the positive spectral measure of [F] for T. If it had
positive mass on [s,1] for any s>exp(-lambda), then the left side
of (3) would be at least nu_F([s,1])s^n. This contradicts (3) as
n tends to infinity. Thus [F] has no spectral support above
exp(-lambda).

The spectral projection onto (exp(-lambda),1] therefore vanishes on
a dense centered family and hence on all of Omega-perp. No bound
uniform in F on K_(F,lambda) is needed. Letting lambda increase
to r gives ||T|Omega-perp||<=exp(-r)=sqrt(alpha_0).
Spectral calculus proves (1).

Finally,

\[
 -\log(1-x)\ge x,\qquad
 x=1-\alpha_0={9263\over1000000},
\]

gives the rational lower bound in (1). The decimal logarithm is
only a display; it is not used in the proof decision.

## Scope

The corollary concerns the selected infinite-volume free-box physical
vacuum at fixed lattice spacing. It proves no uniform periodic/twisted
finite-volume statement, all-boundary uniqueness, full-measure LSI,
continuum mass gap or positive bound along a continuum trajectory.
The argument is written for mathematical review and has not been
formalized in Lean or externally peer reviewed.
