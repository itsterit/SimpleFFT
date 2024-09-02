#include "../../../SimpleFFT.h"
#include "../conversion.h"

/***********************************************************
 *  void fht_dit_iter(double *x, unsigned long n)
 *
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using an iterative radix 2
 *   decimation-in-time FHT.  n must be a power of 2.
 *   Entering this function, x[] must be in bit-reversed
 *   order.  On return, x[] contains the Hartley transform,
 *   returned to normal order.
 *
 * Arguments:
 *   double *x       - array of n doubles, representing n
 *                     real numbers
 *   unsigned long n - dimension of x, must be a power of 2
 ************************************************************/

void fht_dit_iter(fht_dit_type *x, len_type n)
{
#if SIMPLE_FFT__BUILD__FHT_DIT_ITER
    unsigned long m;

    for (m = 2; m <= n; m <<= 1)
    {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        c = 1.0;
        s = 0.0;
        for (j = 1, k = mh - 1; j < mq; ++j, --k)
        {
            double tmp;
            fht_dit_type *xj, *xk;
            xj = x + j + mh;
            xk = x + k + mh;
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            for (i = 0; i < n; i += m)
            {
                double u, v;
                u = xj[i];
                v = xk[i];
                xj[i] = (fht_dit_type)(u * c + v * s);
                xk[i] = (fht_dit_type)(u * s - v * c);
            }
        }
        for (i = 0; i < n; i += m)
        {
            fht_dit_type*xp;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k)
            {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = (fht_dit_type)(u + v);
                xp[k] = (fht_dit_type)(u - v);
            }
        }
    }
#endif /* SIMPLE_FFT__BUILD__FHT_DIT_ITER */
    return;
}

/***********************************************************
 *  void fht_dit_iter_seq(double *x, unsigned long n)
 *
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using an iterative radix 2
 *   decimation-in-time FHT.  n must be a power of 2.
 *   Entering this function, x[] must be in bit-reversed
 *   order.  On return, x[] contains the Hartley transform,
 *   returned to normal order.
 *
 *   The two inner loops of the FHT computation are ordered
 *   to favor sequential memory access at the expense of
 *   redundant trig computations.  See J. Arndt, "Algorithms
 *   for Programmers," online at http://www.jjj.de/fxt/.
 *
 * Arguments:
 *   double *x       - array of n doubles, representing n
 *                     real numbers
 *   unsigned long n - dimension of x, must be a power of 2
 ************************************************************/

void fht_dit_iter_seq(fht_dit_type *x, len_type n)
{
#if SIMPLE_FFT__BUILD__FHT_DIT_ITER_SEQ
    unsigned long m;

    for (m = 2; m <= n; m <<= 1)
    {
        double a, b, t;
        unsigned long i, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m)
        {
            double c, s;
            fht_dit_type *xp;
            unsigned long j, k;
            xp = x + i + mh;
            c = 1.0;
            s = 0.0;
            for (j = 1, k = mh - 1; j < mq; ++j, --k)
            {
                double tmp, u, v;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
                u = xp[j];
                v = xp[k];
                xp[j] = (fht_dit_type)(u * c + v * s);
                xp[k] = (fht_dit_type)(u * s - v * c);
            }
            xp -= mh;
            for (j = 0, k = mh; j < mh; ++j, ++k)
            {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = (fht_dit_type)(u + v);
                xp[k] = (fht_dit_type)(u - v);
            }
        }
    }
#endif /* SIMPLE_FFT__BUILD__FHT_DIT_ITER_SEQ */
    return;
}

/***********************************************************
 *  void fht_dit_rec(double *x, unsigned long n, int nbranch)
 *
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using a radix 2 decimation-in-time
 *   FHT.  If the computation is small enough to fit in
 *   cache, it is done iteratively.  Otherwise, it is done
 *   recursively until the recursion descends to cache-sized
 *   chunks.
 *
 *   n must be a power of 2.  Entering this function, x[]
 *   must be in bit-reversed order.  On return, x[] contains
 *   the Hartley transform, returned to normal order.
 *
 *   To support OpenMP parallelism, nbranch keeps track of
 *   the number of active transforms at a given recursion
 *   level. On the first call to this function, nbranch
 *   should be 1.  It is then doubled for each recursion.
 *
 * Arguments:
 *   double *x       - array of n doubles, representing n
 *                     real numbers
 *   unsigned long n - dimension of x, must be a power of 2
 *   int nbranch     - number of transforms at this recursion
 *                     level
 ************************************************************/

void fht_dit_rec(fht_dit_type *x, len_type n, int nbranch)
{
#if SIMPLE_FFT__BUILD__FHT_DIT_REC
    double a, b, c, s, t;
    uint32_t j, jmax, k, nh, nq;

    if (n == 1)
        return;
    if (n <= (uint32_t)(L1_CACHE_BYTES / sizeof(double)))
    {
        if (FHT_UNIT_STRIDE)
            fht_dit_iter_seq(x, n);
        else
            fht_dit_iter(x, n);
        return;
    }
    nh = (uint32_t)n >> 1;
    nq = nh >> 1;
    nbranch <<= 1;
#if (_OPENMP >= 200203)
    if (nbranch > omp_get_max_threads())
    {
        fht_dit_rec(x, nh, nbranch);
        fht_dit_rec(x + nh, nh, nbranch);
    }
    else
    {
#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            fht_dit_rec(x, nh, nbranch);
#pragma omp section
            fht_dit_rec(x + nh, nh, nbranch);
        }
    }
#else
    fht_dit_rec(x, nh, nbranch);
    fht_dit_rec(x + nh, nh, nbranch);
#endif
    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    jmax = nq + nh;
    c = 1.0;
    s = 0.0;
    for (j = nh + 1, k = (uint32_t)n - 1; j < jmax; ++j, --k)
    {
        double tmp, u, v;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = (fht_dit_type)(u * c + v * s);
        x[k] = (fht_dit_type)(u * s - v * c);
    }
    for (j = 0, k = nh; j < nh; ++j, ++k)
    {
        double u, v;
        u = x[j];
        v = x[k];
        x[j] = (fht_dit_type)(u + v);
        x[k] = (fht_dit_type)(u - v);
    }
#endif /* SIMPLE_FFT__BUILD__FHT_DIT_REC */
    return;
}
