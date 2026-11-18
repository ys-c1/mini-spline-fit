#ifndef MINI_SPLINE_FIT_SPLINE_H
#define MINI_SPLINE_FIT_SPLINE_H

/*
 * Natural cubic spline interpolation.
 *
 * x[0..n-1] : strictly increasing knots
 * y[0..n-1] : function values at the knots
 * m[0..n-1] : output second derivatives at the knots
 *
 * Returns 0 on success and -1 on invalid input.
 */
int spline_natural(const double *x, const double *y, int n, double *m);

/*
 * Evaluate the spline at t.
 *
 * x, y and m must describe an n-knot natural cubic spline previously
 * computed by spline_natural().  Values outside [x[0], x[n-1]] are
 * evaluated by the nearest end interval.
 */
double spline_eval(const double *x, const double *y, const double *m,
                   int n, double t);

#endif /* MINI_SPLINE_FIT_SPLINE_H */
