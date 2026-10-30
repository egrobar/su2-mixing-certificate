# Release verification record

Date: 2026-10-29. This records project checks, not independent human
refereeing or formal proof verification.

Three separately tasked smaller-model agents audited the analytic proof,
reproduced arithmetic/geometry, and checked provenance and scope. The
primary agent checked the physical lift and transfer argument, assembled
this standalone package, and reviewed the audit conclusions.

The analytic audit found no fatal gap in the stated projected/free-box
claim. It independently checked the conformal cancellation, the vMF
Poincare trial residuals, the full chord interpolation, occurrence
counting, and the Dobrushin covariance implication. The explicit
comparison proof is in [verification scope](../docs/verification-scope.md).
The transfer implication remains a written analytic argument.

The 4096-cell rational calculation reproduces an upper endpoint of
approximately 0.9907366241076274 at beta_W=77/200. All 19 standalone
tests pass. An independently written Fraction implementation, without
importing the certificate modules, also verifies contraction with a
coarser 1024-cell estimate below 0.997328724. That independent calculation
checks contraction; the sharper displayed 0.990737 bound comes from the
included 4096-cell implementation and its exact arithmetic tests.

Both wrapped local counts and an unwrapped parity enumeration confirm
the 6+60 occurrence bound and midpoint l1 range at most two. The forest
argument for arbitrary finite free rectangles is analytic and stated
separately; periodic cells are not treated as gauge-fixable forests.

The arithmetic/geometry code and four focused test modules were extracted
from the research implementation, with line endings normalized for this
release. Original-source and release-file hashes are recorded separately
in `provenance-sha256.json`; current code hashes are also recorded by the
verifier. No broader research notes or saved field configurations are
needed to run the package.

The original research suite passed 97 tests during extraction. The
public package deliberately includes the 19 tests relevant to this
certificate; its test count is not 97.

No Lean build of this theorem was performed. Upstream MIT notices and
immutable source pins were verified. Publication does not establish
novelty, priority, external peer review, or a continuum gap.

The GitHub Actions workflow is provided for reproducibility. A local
successful run is not evidence that hosted CI has run; its status should
be checked on the published repository.
