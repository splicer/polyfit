#include "polyfit.h"
#include <stdio.h>


static void print_weights(double *w, unsigned count)
{
    for (unsigned i = 0; i < count; i++) {
        printf(" %f", w[i]);
    }
    putchar('\n');
}


static double f(double x)
{
    return 3.2 - 12.5 * x + 0.223 * x * x * x;
}


int main(void)
{
    polyfit_t *cubic_fit = polyfit_create(3);

    for (unsigned i = 0; i < 10; i++) {
        polyfit_add_point(cubic_fit, i, f(i));

        double weights[POLYFIT_NUM_WEIGHTS(3)];
        polyfit_get_weights(cubic_fit, weights);
        print_weights(weights, POLYFIT_NUM_WEIGHTS(3));
    }

    printf("f(%f) = %f\n", 3.5, f(3.5));
    printf("f(%f) ~ %f\n", 3.5, polyfit_estimate_y(cubic_fit, 3.5));

    double cubic_fit_archive[POLYFIT_ARCHIVE_LEN(3)];
    polyfit_archive(cubic_fit, cubic_fit_archive);
    polyfit_destroy(cubic_fit);

    polyfit_t *new_cubic_fit = polyfit_create_from_archive(3, cubic_fit_archive);
    printf("f(%f) ~ %f\n", 3.5, polyfit_estimate_y(new_cubic_fit, 3.5));

    return 0;
}
