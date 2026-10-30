# Sharp single-link transport and Poincare inputs

This exposition adopts the pinned-source argument. The project checked the harmonic reduction, Bessel identities,
coefficient signs, Riccati barriers, and auxiliary Poincare estimate in the
[PRD manuscript at commit b2d359f](https://github.com/creelie/4Y-M/blob/b2d359f1197f9f50aae8f6890edb6d9106531fa3/strong-coupling/prd/main.tex).
The argument below is an adopted and checked source argument, not an
originality claim, peer review, or a Lean formalization.

For normalized Haar measure on the unit sphere \(S^3\), set
\[
 \nu_c(dx)=Z(c)^{-1}e^{c\cdot x}\,d\sigma(x).
\]
In Euclidean chordal transport,
\[
 \boxed{W_1(\nu_c,\nu_d)\le |c-d|/4.}
\]
The constant is sharp at Haar. This is not the corresponding theorem for
the intrinsic geodesic transport metric.

## Covariance and harmonic reduction

For each chordal 1-Lipschitz \(f\),
\[
 {d\over dt}\nu_{c+th}(f)=\operatorname{Cov}_{\nu_{c+th}}(f,h\cdot x).
\]
Thus a uniform covariance bound \(1/4\) for unit \(h\) gives finite endpoint
transport by integration along the straight Euclidean parameter segment.
No conversion from an input geodesic to a chord is used.

Let \(u=(h\cdot x-\nu_c(h\cdot x))\,d\nu_c/d\sigma\), and let \(u_n\)
be its degree-\(n\) spherical harmonic component. Define
\[
 Q(u)=\sum_{n\ge1}{\|u_n\|_2^2\over n}.
\]
Extend \(f\) to a Euclidean 1-Lipschitz function on the unit four-ball,
and use the harmonic function with normal derivative \(u\).
Green's identity and Cauchy--Schwarz, using
\(|B^4|/|S^3|=1/4\), give
\[
 |\operatorname{Cov}_{\nu_c}(f,h\cdot x)|\le\tfrac12\sqrt{Q(u)}.
\]
The exponential density is analytic, so the harmonic series and the
required derivatives converge. Lipschitz extensions have ambient gradient
bounded by one almost everywhere; integration by parts follows by
approximation.

Rotate \(c=k e_1\), \(h=\cos\phi\,e_1+\sin\phi\,e_2\), and put
\(d_n=(n+1)I_{n+1}(k)/I_1(k)\).
Chebyshev orthogonality on \(S^3\), and odd/even parity in \(x_2\), yield
\[
 Q(u)=\cos^2\phi\,Q_z+\sin^2\phi\,Q_x,\quad
 Q_z=\sum_{n\ge1}{d_n'^2\over n},\quad
 Q_x={1\over3k^2}\sum_{n\ge1}(n+2)d_n^2.
\]
Both limits at \(k=0\) are \(1/4\). At that point \(f=h\cdot x\)
attains the resulting covariance \(1/4\).

## The transverse energy, with all coefficients controlled

Parseval and the Bessel recurrence give
\[
 \sum_{j\ge1}j^2 I_j^2={k\over4}I_1(2k),\qquad
 \sum_{j\ge1}j^3 I_j^2={k^2\over4}(I_0^2+I_1^2),\qquad
 I_1(2k)=2\sum_{j\ge1}I_{j-1}I_j.
\]
The second identity follows by substituting
\(jI_j=k(I_{j-1}-I_{j+1})/2\) and telescoping; all series converge
absolutely. It follows that \(Q_x\le1/4\) is equivalent to
\[
 D(k)=(8+2k^2)I_1^2-k^2I_0^2-kI_1(2k)\ge0.
\]
For \(m\ge1\), the coefficient of \(k^{2m+2}\) is
\[
 {N_m\over4^m m!(m+1)!(m+1)(m+2)},\quad
 N_m={2m\choose m}\{4(2m+1)+(m-1)(m+1)(m+2)\}
       -4^m(m+1)(m+2).
\]
This identity follows by convolution of the positive Bessel power series
and Vandermonde's identity. The earlier coefficients vanish.
One finds \(N_1=N_2=0\), \(N_3=80\), \(N_4=1140\), \(N_5=10416\).
For \(m\ge6\), induction gives
\(16^m\le4m{2m\choose m}^2\), and \((m-1)^2\ge4m\).
Consequently \((m-1){2m\choose m}\ge4^m\), which makes \(N_m>0\).
This proves the entire-series inequality; it is not a finite numerical
check of a few coefficients.

## The zonal energy and large-parameter bound

Write \(A=I_2/I_1\). Integration by parts gives
\(A'=1-A^2-3A/k=\operatorname{Var}(x_1)\).
Comparison for this Riccati equation gives
\[
 {k\over2+\sqrt{k^2+4}}\le A\le
 {k\over1+\sqrt{k^2+9}}\le{k\over4},
 \qquad A'\le A/k,
\]
and the stronger lower barrier
\[
 A_\flat={k\over3/2+\sqrt{k^2+25/4}}\le A.
\]
These comparisons are valid at the singular endpoint: for a barrier
\(\varphi\), multiply the difference equation by
\(k^3\exp\int_0^k(A+\varphi)\); the initial product tends to zero, and
the explicitly signed barrier residual determines the sign.

The positive-series ratio argument gives
\(I_{n+2}/I_{n+1}\le A\). Therefore
\(|d_n'|\le n d_n/k\) for \(n\ge2\), \(k\le4\).
The preceding Bessel sums, retaining three positive terms in \(I_1(2k)\),
give, for \(0<k\le\sqrt{88}/3\),
\[
 Q_z\le\tfrac14+4A'^2+T(A),\quad
 T(a)=a^2/4-a/(2k)-2a^2/k^2.
\]
Here \(T\) decreases on \([0,k/4]\), since \(k^2<12\).
Put \(s=\sqrt{k^2+25/4}\) and \(\alpha=16/(2s+3)^2\).
Then \(0\le A'\le\alpha\), and
\[
 4\alpha^2+T(A_\flat)
 ={(2s-5)(8s^3+28s^2-98s-695)\over4(2s+3)^4}<0.
\]
Indeed \(5/2<s\le\sqrt{577}/6<4.004\); the cubic increases there
and is still below \(-120\) at \(4.004\).
Together with the transverse estimate, this proves the quarter on
the stated compact parameter range.

For larger \(k\), use the **auxiliary** Poincare bound
\(\operatorname{Var}_{\nu_c}f\le\frac13\nu_c|\nabla f|^2\).
Here is the checked proof input. In polar angle \(\theta\), the generator is
\[
 L=\partial_\theta^2+(2\cot\theta-k\sin\theta)\partial_\theta
       +\sin^{-2}\theta\,\Delta_{S^2}.
\]
In the first nonzonal sector, the positive trial function
\(\psi=\sin\theta\,e^{a k\cos\theta}\), \(a=-1/5\), gives
\(-L_1\psi/\psi=3+4k^2\sin^2\theta/25\ge3\).
Higher angular sectors have larger potential. Differentiating a
nonconstant zonal eigenfunction gives the potential
\(2/\sin^2\theta+k\cos\theta\); using \(a=-2/5\) gives
\(3+6k^2\sin^2\theta/25\).
The ground-state integration by parts has vanishing endpoint terms:
smooth sector functions and zonal derivatives vanish at least linearly,
and the polar density is \(\sin^2\theta e^{k\cos\theta}\).
Compact ellipticity supplies a discrete spectrum.
Thus every nonzero eigenvalue is at least three, uniformly in \(c\).

For chordal 1-Lipschitz \(f\), its tangential gradient is bounded by one.
The largest eigenvalue of the coordinate covariance is \(A/k\), so
\[
 |\operatorname{Cov}(f,h\cdot x)|
 \le [3(1+\sqrt{k^2+9})]^{-1/2}\le1/4
 \quad(k\ge\sqrt{88}/3).
\]
This completes the all-parameter proof. The auxiliary spectral gap used
here estimates conditional variances; it is not the physical Yang--Mills
Hamiltonian gap.
