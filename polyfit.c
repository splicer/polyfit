#include "polyfit.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define SMALL_VALUE 1.0E-32
#define FORGETTING_FACTOR (1.0 - 1.0E-11)

struct polyfit_t {
    unsigned rows;
    unsigned cols;
    bool weights_need_update;
    double *weights;
    double cells[];
};


// givens generation
static void boundary_cell(double *restrict cell,
                          double *restrict c,
                          double *restrict s,
                          double in)
{
    if (fabs(in) < SMALL_VALUE) {
        // close enough to zero
        *c = 1;
        *s = 0;
    } else {
        double norm = sqrt(*cell * *cell + in * in);
        *c = *cell / norm;
        *s = in / norm;
        *cell = FORGETTING_FACTOR * norm;
    }
}


// givens rotation
static double internal_cell(double *cell, double c, double s, double in)
{
    double out = c * in - FORGETTING_FACTOR * s * *cell;
    *cell = s * in + FORGETTING_FACTOR * c * *cell;
    return out;
}


polyfit_t * polyfit_create(unsigned degree)
{
    const unsigned rows = degree + 1;
    const unsigned cols = degree + 2;

    polyfit_t *p = malloc(sizeof(polyfit_t) + sizeof(double) * rows * (cols + 1));
    if (!p) return NULL;

    p->rows = rows;
    p->cols = cols;
    p->weights_need_update = true;
    p->weights = &p->cells[rows * cols];

    for (unsigned i = 0; i < p->rows * p->cols; i++) {
        p->cells[i] = SMALL_VALUE;
    }

    return p;
}


polyfit_t * polyfit_create_from_archive(unsigned degree, const double *archive)
{
    polyfit_t *p = polyfit_create(degree);
    if (!p) return NULL;

    for (unsigned i = 0, arch_index = 0; i < p->rows; i++) {
        for (unsigned j = 0; j < p->cols; j++) {
            if (j < i) continue; // skip unused cells
            p->cells[i * p->cols + j] = archive[arch_index++];
        }
    }

    return p;
}


void polyfit_archive(const polyfit_t *p, double *archive)
{
    for (unsigned i = 0, arch_index = 0; i < p->rows; i++) {
        for (unsigned j = 0; j < p->cols; j++) {
            if (j < i) continue; // skip unused cells
            archive[arch_index++] = p->cells[i * p->cols + j];
        }
    }
}


void polyfit_destroy(polyfit_t *p)
{
    free(p);
}


void polyfit_add_point(polyfit_t *p, double x, double y)
{
    double in[p->cols];

    in[0] = 1;
    for (unsigned j = 1; j < p->cols - 1; j++) {
        in[j] = in[j - 1] * x;
    }
    in[p->cols - 1] = y;

    for (unsigned i = 0; i < p->rows; i++) {
        double c, s;
        boundary_cell(&p->cells[i * p->cols + i], &c, &s, in[i]);
        for (unsigned j = i + 1; j < p->cols; j++) {
            double out = internal_cell(&p->cells[i * p->cols + j], c, s, in[j]);
            if (i < p->rows - 1) in[j] = out;
        }
    }

    p->weights_need_update = true;
}


static void compute_weights(polyfit_t *p)
{
    if (!p->weights_need_update) return;

    for (int i = p->rows - 1; i >= 0; i--) {
        p->weights[i] = p->cells[i * p->cols + p->cols - 1];
        for (unsigned j = i + 1; j < p->cols - 1; j++) {
            p->weights[i] -= p->cells[i * p->cols + j] * p->weights[j];
        }
        p->weights[i] /= p->cells[i * p->cols + i];
    }

    p->weights_need_update = false;
}


void polyfit_get_weights(polyfit_t *p, double *weights)
{
    compute_weights(p);
    for (unsigned i = 0; i < p->rows; i++) {
        weights[i] = p->weights[i];
    }
}


double polyfit_estimate_y(polyfit_t *p, double x)
{
    unsigned i = p->rows - 1;
    double y;

    compute_weights(p);

    y = p->weights[i];
    while (i-- > 0) {
        y = y * x + p->weights[i];
    }

    return y;
}
