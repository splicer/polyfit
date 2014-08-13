#ifndef POLYFIT_H
#define POLYFIT_H

#define POLYFIT_NUM_WEIGHTS(deg) ((deg) + 1)
#define POLYFIT_ARCHIVE_LEN(deg) (((deg) + 4) * ((deg) + 1) / 2)

typedef struct polyfit_t polyfit_t;

/* Create a polyfit_t object (on the heap) of given degree.
 */
polyfit_t * polyfit_create(unsigned degree);

/* Create a polyfit_t object (on the heap) of given degree using the values
 * found in an archive produced by polyfit_archive.
 */
polyfit_t * polyfit_create_from_archive(unsigned degree, const double *archive);

/* Outputs vector of values that can be used to recreate a polyfit_t object
 * using polyfit_create_from_archive. This is useful if you'd like to save the
 * state of a polyfit_t object to a file in a portable fashion. Note that for a
 * polynomial of degree deg, archive must contain at least
 * POLYFIT_ARCHIVE_LEN(deg) elements.
 */
void polyfit_archive(const polyfit_t *p, double *archive);

/* Frees the memory associated with object p.
 */
void polyfit_destroy(polyfit_t *p);

/* Update object p's estimate of the underlying polynomial with an x-y pair.
 */
void polyfit_add_point(polyfit_t *p, double x, double y);

/* Copies the coefficient vector represented by object p into array weights.
 * Note that for a polynomial of degree deg, weights must contain at least
 * POLYFIT_NUM_WEIGHTS(deg) elements. weights are ordered by increasing power
 * (i.e. weights[2] is the coefficient for x^2.
 */
void polyfit_get_weights(polyfit_t *p, double *weights);

/* Let P be the polynomial represented by the object p. This function returns
 * P(x).
 */
double polyfit_estimate_y(polyfit_t *p, double x);

#endif
