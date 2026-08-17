# Next step: history-dependent survival and MFPT

The current `ComputeSegmentSurvivalFunction()` is still an **instantaneous geometric occupancy** observable. It can leave the tube and later return.

For mean-first-passage analysis, do not infer first passage from the already averaged instantaneous curve. The history has to be tracked at the trajectory/time-origin level.

Recommended next implementation:

1. Keep the current instantaneous function as a regression/reference observable.
2. Add a history-dependent first-passage segment survival function that scans every stored intermediate frame for each time origin.
3. Define an absorbing first-passage time for each reference contour sample and tube diameter, e.g.

   T_FP(s,a) = min{ t : escape criterion is first satisfied }.

4. The first-passage survival curve is

   S_FP(s,t;a) = P[T_FP(s,a) > t].

5. The per-contour mean first-passage time is the area under that absorbing survival curve:

   MFPT(s,a) = integral_0^infinity S_FP(s,t;a) dt.

   With a finite trajectory, use a restricted MFPT (RMFPT) up to the observation horizon unless a tail model is explicitly fitted.

6. The contour-integrated tube survival

   mu_FP(t;a) = integral_0^1 S_FP(s,t;a) ds

   remains a survival probability curve. Its time integral is the contour-averaged MFPT. `TubeSurvivalFunction` itself should not be relabeled as an MFPT scalar.

7. The affine correction must be applied at every intermediate frame using the same time-origin box:

   future/intermediate frame box -> box at t0

   before the escape/first-passage criterion is evaluated.

This separation keeps the old nonabsorbing result available for comparison and avoids silently changing the meaning of existing CSVs/tests.
