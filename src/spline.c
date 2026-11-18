#include "spline.h"

#include <stddef.h>
#include <stdlib.h>

/*
 * Return the index i of the interval [x[i], x[i+1]] containing t.
 * End intervals are used for points outside the knot range.
 */
static int locate_interval(const double *x, int n, double t)
{
    int i;

    for (i = 0; i < n - 1; ++i) {
        if (x[i] < t && t < x[i + 1]) {
            return i;
        }
    }

    return 0;
}

int spline_natural(const double *x, const double *y, int n, double *m)
{
    int k;
    int j;
    double *diag;
    double *upper;
    double *rhs;

    if (n < 3 || x == NULL || y == NULL || m == NULL) {
        return -1;
    }

    for (j = 0; j < n - 1; ++j) {
        if (!(x[j + 1] > x[j])) {
            return -1;
        }
    }

    k = n - 2; /* number of unknown interior second derivatives */

    diag = (double *)malloc((size_t)k * sizeof(*diag));
    upper = (double *)malloc((size_t)k * sizeof(*upper));
    rhs = (double *)malloc((size_t)k * sizeof(*rhs));
    if (diag == NULL || upper == NULL || rhs == NULL) {
        free(diag);
        free(upper);
        free(rhs);
        return -1;
    }

    /*
     * Interior equation for M[j+1]:
     *
     * h[j] * M[j] + 2(h[j] + h[j+1]) * M[j+1]
     *     + h[j+1] * M[j+2] = 6 (slope[j+1] - slope[j])
     *
     * Unknowns are numbered locally as M[1] ... M[n-2].
     */
    for (j = 0; j < k; ++j) {
        double h_left = x[j + 1] - x[j];
        double h_right = x[j + 2] - x[j + 1];
        double slope_left = (y[j + 1] - y[j]) / h_left;
        double slope_right = (y[j + 2] - y[j + 1]) / h_right;

        diag[j] = 2.0 * (h_left + h_right);
        upper[j] = h_right;
        rhs[j] = 6.0 * (slope_right - slope_left);
    }

    /* Thomas algorithm: forward elimination. */
    for (j = 1; j < k; ++j) {
        double a = x[j + 1] - x[j];
        double b = x[j + 2] - x[j + 1];
        double factor = a / b;

        diag[j] -= factor * upper[j - 1];
        rhs[j] -= factor * rhs[j - 1];
    }

    /* Thomas algorithm: back substitution. */
    m[n - 2] = rhs[k - 1] / diag[k - 1];
    for (j = k - 2; j >= 0; --j) {
        m[j + 1] = (rhs[j] - upper[j] * m[j + 2]) / diag[j];
    }

    /* Natural boundary conditions. */
    m[0] = 0.0;
    m[n - 1] = 0.0;

    free(diag);
    free(upper);
    free(rhs);
    return 0;
}

double spline_eval(const double *x, const double *y, const double *m,
                   int n, double t)
{
    int i;
    double h;
    double a;
    double b;

    if (n < 2 || x == NULL || y == NULL || m == NULL) {
        return 0.0;
    }

    i = locate_interval(x, n, t);
    h = x[i + 1] - x[i];
    a = (x[i + 1] - t) / h;
    b = (t - x[i]) / h;

    /*
     * M[i] and M[i+1] are second derivatives.  The cubic correction
     * contains both a cubic term and a linear term.
     */
    return a * y[i] + b * y[i + 1]
         + (h * h / 6.0) * (a * a * a * m[i]
                          + b * b * b * m[i + 1]);
}
