# Tilt lemmas and rigorous derivative enclosures

The [analytic inputs](analytic-inputs.md) supply the vMF quarter transport,
Poincare constant 1/3, and global Hessian upper bound tau=1/4. This note
supplies the perturbation estimates used in the main theorem and explains
the integer interval calculation. The decimation checker uses field
radius 12/5, 4096 closed cells, polynomial degree 24 and scale 2^120.
These are the release's 0.385 parameters, not the smaller field domain
used by the earlier three-link experiment.

## The variance sharpening

For a probability law \(\nu\) and its tilt
\(\kappa=e^R\nu/\nu(e^R)\), suppose \(\operatorname{osc}R\le w\).
Writing \(\rho=d\kappa/d\nu\), \(\delta=\|\rho-1\|_\infty\le e^w-1\),
direct centering gives
\[
 |\operatorname{Cov}_\kappa(f,h)-\operatorname{Cov}_\nu(f,h)|
 \le(\delta+\delta^2)\operatorname{sd}_\nu(f)\operatorname{sd}_\nu(h).
\]
If \(\nu\) is a vMF law and \(f\) is chordal 1-Lipschitz, its tangential
gradient is bounded by one, so \(\operatorname{sd}_\nu(f)\le1/\sqrt3\).
For \(h(x)=v\cdot x\), its standard deviation is at most \(\sqrt\tau|v|\).
This yields the sharper constant
\[
 s=\tfrac14+(e^w-1)e^w\sqrt{\tau/3}.
\]
The source used the diameter bound \(\operatorname{sd}_\nu(f)\le1\).

There is a second improvement. If \(a\le\rho\le d\), then \(d/a\le e^w\).
The variational formula for variance and the vMF Poincare bound give
\[
 \operatorname{Var}_\kappa f
 \le d\,\operatorname{Var}_\nu f
 \le {d\over3}\nu|\nabla f|^2
 \le {e^w\over3}\kappa|\nabla f|^2.
\]
Consequently the covariance of two chordal Lipschitz functions with
constants \(L_1,L_2\) is at most \(e^w L_1L_2/3\).
This derivation makes the perturbation factor explicit; it does not assume
a spectral gap for the multi-link block.

## Interval evaluation

The four positive series derivatives are evaluated with integer interval
endpoints, rounding every division and multiplication outward. The
denominator interval must stay strictly positive. All derivative cells
include both endpoints and cover the entire field interval, including
zero. The exponential upper bound uses positive terms and a geometric
tail. No floating-point comparison decides success.

Writing \(s_r=r^2\), the exact series is
\[
 Z(s_r)=\sum_{m\ge0}{s_r^m\over4^m m!(m+1)!},\quad
 Z^{(j)}(s_r)=\sum_{m\ge0}
 {s_r^m\over4^{m+j}m!(m+j+1)!},\quad 0\le j\le3.
\]
Twenty-five terms (degrees 0 through 24) and a positive geometric tail
bound enclose each derivative. Interval quotient formulas give the
derivatives of \(\varphi=\log Z\), including the origin without division
by \(r\).

The Hessian eigenvalues are \(2\varphi'\) and
\(2\varphi'+4r^2\varphi''\).
For a unit direction, the cubic derivative is \(At+\Gamma t^3\), with
\(A=12r\varphi''\), \(\Gamma=8r^3\varphi'''\), \(|t|\le1\).
Its norm is **at most**
\(\max(|A+\Gamma|,2|A|/3)\): endpoint values and interior critical values
give that bound. The manuscript calls this an equality, which is generally
too strong; only the valid upper bound is used here. The equality of the
symmetric trilinear norm with its diagonal norm on a real Hilbert space
is the standard Banach theorem.

The exponential is bounded by another positive rational series and
geometric remainder. To avoid a square-root routine we use
\(289/1000>\sqrt{1/12}\), verified by squaring. This loses a small amount.

The complete code is [su2_vmf_certificate.py](../experiments/su2_vmf_certificate.py).
Its earlier three-link command remains available as a dependency check;
the release decision is made by [verify.py](../verify.py) using the radius
12/5 decimation calculation. Endpoint monotonicity in beta is justified
analytically in [the theorem](theorem.md), not inferred from sampled beta
values or interval output monotonicity.
