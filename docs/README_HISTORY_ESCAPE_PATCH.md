# History-dependent permanent tube escape patch

This patch adds a history-dependent tube survival estimator on top of the
existing affine-corrected geometry.

## Rule

For tube diameter `a` and affine-corrected distance `d` from a reference
contour sample to the future primitive path:

- `d <= a/2`: inside the inner tube.
- `a/2 < d < a`: temporarily outside; recovery is allowed.
- `d >= a`: permanently escaped for that time origin.

The returned `HistoryDependentSurvivalResult` contains:

1. `stephanou_survival`
   - 1 only when the sample is currently inside `a/2`
   - and has never previously reached `d >= a`.
   - This is the curve used as the third member of the smoke-test comparison.

2. `permanent_escape_survival`
   - pure absorbing outer-boundary survival
   - 1 while the running maximum distance remains below `a`
   - 0 forever after `d >= a`.

Every saved intermediate frame is scanned, even when it is not a requested
output lag, so an escape between output lags cannot be missed.

## Files

Replace the project-root files:

- `Tube_Survival.hpp`
- `Tube_Survival.cpp`

Add:

- `tests/History_Dependent_Escape_Smoke_Test.cpp`
- `run_history_escape_smoke_test.sh`

## Run

From the project root:

```bash
BOX_CENTER_Z=39.22658609485464 \
bash run_history_escape_smoke_test.sh Z1+SP.dat 20 5.0
```

The real-data smoke test defaults to `CHAIN_INDEX=0` and lags through frame 80
(1600 tau_LJ for 20 tau_LJ spacing) so the history scan remains modest.

Choose another chain with:

```bash
CHAIN_INDEX=10 \
BOX_CENTER_Z=39.22658609485464 \
bash run_history_escape_smoke_test.sh Z1+SP.dat 20 5.0
```

## Main three-way comparison

The smoke test compares:

1. raw instantaneous tube survival;
2. affine-corrected instantaneous tube survival;
3. affine-corrected history-dependent survival with permanent kill at `d >= a`.

It also writes the pure absorbing outer-boundary survival as a separate
diagnostic.

## Regression tests

The fast synthetic regression uses distance history

`0, 7, 3, 11, 2` with `a=10`.

Expected history-dependent inner-tube sequence:

`1, 0, 1, 0, 0`

Expected pure permanent-escape sequence:

`1, 1, 1, 0, 0`

This verifies temporary recovery before the outer escape boundary and
irreversibility after permanent escape.
