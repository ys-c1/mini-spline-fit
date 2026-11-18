# mini-spline-fit

C11 natural cubic spline interpolation with no third-party dependencies.

## API

- `int spline_natural(const double *x, const double *y, int n, double *m)`
  computes the knot second derivatives `m` using a tridiagonal Thomas solve
  with natural boundary conditions (`m[0] = m[n-1] = 0`). It returns `-1`
  for `n < 3` or invalid arguments.
- `double spline_eval(const double *x, const double *y, const double *m, int n, double t)`
  evaluates the piecewise cubic at `t`.

## Build and test

```sh
make test
```

A correctly implemented build reports the following summary:

```text
checks: 77, failures: 0
```

Only `src/spline.c` and `src/spline.h` are implementation files. The tests in
`tests/test_spline.c` define the expected behavior and are not part of the
implementation.

This seed repository intentionally contains algorithmic defects in
`src/spline.c`; the final checked-in state is the failing variant and currently
reports `checks: 77, failures: 16`. The test file must not be modified to make
the suite pass.
