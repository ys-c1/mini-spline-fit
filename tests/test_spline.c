#include "spline.h"

#include <math.h>
#include <stdio.h>

static int checks = 0;
static int failures = 0;

#define CHECK(condition, label)                                               \
    do {                                                                      \
        ++checks;                                                             \
        if (!(condition)) {                                                   \
            ++failures;                                                       \
            printf("FAIL [%s] %s\n", (label), #condition);                    \
        }                                                                     \
    } while (0)

#define CHECK_NEAR(actual, expected, tolerance, label)                        \
    do {                                                                      \
        double check_actual = (actual);                                       \
        double check_expected = (expected);                                   \
        double check_tolerance = (tolerance);                                 \
        double check_error = fabs(check_actual - check_expected);             \
        ++checks;                                                             \
        if (!(check_error <= check_tolerance)) {                              \
            ++failures;                                                       \
            printf("FAIL [%s] actual=%.17g expected=%.17g error=%.6e\n",      \
                   (label), check_actual, check_expected, check_error);       \
        }                                                                     \
    } while (0)

static double derivative_from_left(const double *x, const double *y,
                                   const double *m, int i)
{
    double h = x[i] - x[i - 1];

    return (y[i] - y[i - 1]) / h + h * (2.0 * m[i] + m[i - 1]) / 6.0;
}

static double derivative_from_right(const double *x, const double *y,
                                    const double *m, int i)
{
    double h = x[i + 1] - x[i];

    return (y[i + 1] - y[i]) / h - h * (2.0 * m[i] + m[i + 1]) / 6.0;
}

static void test_three_points(void)
{
    const double x[3] = {0.0, 1.0, 2.0};
    const double y[3] = {0.0, 1.0, 0.0};
    double m[3];
    int i;

    CHECK(spline_natural(x, y, 3, m) == 0, "three-point solve");

    for (i = 0; i < 3; ++i) {
        CHECK_NEAR(spline_eval(x, y, m, 3, x[i]), y[i], 1e-12,
                   "three-point knot value");
    }

    CHECK_NEAR(spline_eval(x, y, m, 3, 0.5), 0.6875, 1e-12,
               "three-point first midpoint");
    CHECK_NEAR(spline_eval(x, y, m, 3, 1.5), 0.6875, 1e-12,
               "three-point second midpoint");
}

static void test_quadratic(void)
{
    const double h = 1e-4;
    const double x[4] = {0.0, h, 2.0 * h, 3.0 * h};
    double y[4];
    double m[4];
    double t = 1.5 * h;
    int i;

    for (i = 0; i < 4; ++i) {
        y[i] = x[i] * x[i];
    }

    CHECK(spline_natural(x, y, 4, m) == 0, "quadratic solve");

    for (i = 0; i < 4; ++i) {
        CHECK_NEAR(spline_eval(x, y, m, 4, x[i]), y[i], 1e-12,
                   "quadratic knot value");
    }

    CHECK_NEAR(spline_eval(x, y, m, 4, t), t * t, 1e-9,
               "quadratic midpoint");
}

static void test_natural_boundary(void)
{
    const double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    const double y[5] = {0.0, 1.0, 0.0, 1.0, 0.0};
    double m[5];

    CHECK(spline_natural(x, y, 5, m) == 0, "natural-boundary solve");
    CHECK_NEAR(m[0], 0.0, 1e-15, "left natural boundary");
    CHECK_NEAR(m[4], 0.0, 1e-15, "right natural boundary");
}

static void test_known_second_derivatives(void)
{
    const double x[4] = {0.0, 1.0, 2.0, 3.0};
    const double y[4] = {0.0, 1.0, 4.0, 9.0};
    double m[4];

    CHECK(spline_natural(x, y, 4, m) == 0, "known-M solve");
    CHECK_NEAR(m[1], 2.4, 1e-12, "known first interior second derivative");
    CHECK_NEAR(m[2], 2.4, 1e-12, "known second interior second derivative");
}

static void test_smoothness(void)
{
    const double x[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const double y[6] = {0.0, 1.15, 2.60, 4.35, 6.40, 8.75};
    double m[6];
    int i;

    CHECK(spline_natural(x, y, 6, m) == 0, "smoothness solve");

    for (i = 1; i <= 4; ++i) {
        double left = derivative_from_left(x, y, m, i);
        double right = derivative_from_right(x, y, m, i);

        CHECK_NEAR(left, right, 1e-9, "first derivative continuity");
    }
}

static void test_monotonicity(void)
{
    const double x[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const double y[6] = {0.0, 1.15, 2.60, 4.35, 6.40, 8.75};
    double m[6];
    double previous;
    int sample;

    CHECK(spline_natural(x, y, 6, m) == 0, "monotonicity solve");

    previous = spline_eval(x, y, m, 6, 0.0);
    for (sample = 1; sample <= 50; ++sample) {
        double t = (double)sample * 0.1;
        double current = spline_eval(x, y, m, 6, t);

        CHECK(current >= previous - 1e-12, "sampled monotonicity");
        previous = current;
    }
}

static void test_short_inputs(void)
{
    const double x[2] = {0.0, 1.0};
    const double y[2] = {0.0, 1.0};
    double m[2];

    CHECK(spline_natural(x, y, 2, m) == -1, "n=2 rejected");
    CHECK(spline_natural(x, y, 1, m) == -1, "n=1 rejected");
    CHECK(spline_natural(x, y, 0, m) == -1, "n=0 rejected");
}

int main(void)
{
    test_three_points();
    test_quadratic();
    test_natural_boundary();
    test_known_second_derivatives();
    test_smoothness();
    test_monotonicity();
    test_short_inputs();

    printf("checks: %d, failures: %d\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
