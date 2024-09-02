#include "../../../SimpleFFT.h"
#include "../conversion.h"

/***********************************************************
 *  void fht_dif_iter(double *x, unsigned long n)
 *
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using an iterative radix 2
 *   decimation-in-frequency FHT.  n must be a power of 2.
 *   Entering this function, x[] is in normal order.  On
 *   return, x[] contains the Hartley transform, stored in
 *   bit-reversed order.
 *
 * Arguments:
 *   double *x       - array of n doubles, representing n
 *                     real numbers
 *   unsigned long n - dimension of x, must be a power of 2
 ************************************************************/

void fht_dif_iter(fht_dif_type *x, len_type n)
{
#if SIMPLE_FFT__BUILD__FHT_DIF_ITER
    unsigned long m;

    for (m = n; m > 1; m >>= 1)
    {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m)
        {
            double *xp;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k)
            {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
        }
        c = 1.0;
        s = 0.0;
        for (j = 1, k = mh - 1; j < mq; ++j, --k)
        {
            double tmp;
            double *xj, *xk;
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
                xj[i] = u * c + v * s;
                xk[i] = u * s - v * c;
            }
        }
    }
#endif /* SIMPLE_FFT__BUILD__FHT_DIF_ITER */
    return;
}

/***********************************************************
 *  void fht_dif_iter_seq(double *x, unsigned long n)
 *
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using an iterative radix 2
 *   decimation-in-frequency FHT.  n must be a power of 2.
 *   Entering this function, x[] is in normal order.  On
 *   return, x[] contains the Hartley transform, stored in
 *   bit-reversed order.
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

void fht_dif_iter_seq(fht_dif_type *x, len_type n)
{
#if SIMPLE_FFT__BUILD__FHT_DIF_ITER_SEQ
    unsigned long m;

    for (m = n; m > 1; m >>= 1)
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
            double *xp;
            unsigned long j, k;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k)
            {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
            xp += mh;
            c = 1.0;
            s = 0.0;
            for (j = 1, k = mh - 1; j < mq; ++j, --k)
            {
                double u, v, tmp;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
                u = xp[j];
                v = xp[k];
                xp[j] = u * c + v * s;
                xp[k] = u * s - v * c;
            }
        }
    }
#endif /* SIMPLE_FFT__BUILD__FHT_DIF_ITER_SEQ */
    return;
}

/***********************************************************
 *  void fht_dif_rec(double *x, unsigned long n, int nbranch)
 *
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using a radix 2 decimation-in-
 *   frequency FHT.  If the computation is small enough to
 *   fit in cache, it is done iteratively.  Otherwise, it is
 *   done recursively until the recursion descends to
 *   cache-sized chunks.
 *
 *   n must be a power of 2.  Entering this function, x[]
 *   must be in normal order.  On return, x[] contains the
 *   Hartley transform, in bit-reversed order.
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

void fht_dif_rec(fht_dif_type *x, len_type n, int nbranch)
{
#if SIMPLE_FFT__BUILD__FHT_DIF_REC
    double a, b, c, s, t;
    unsigned long j, jmax, k, nh, nq;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / sizeof(double)))
    {
        if (FHT_UNIT_STRIDE)
            fht_dif_iter_seq(x, n);
        else
            fht_dif_iter(x, n);
        return;
    }
    nh = n >> 1;
    nq = nh >> 1;
    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    for (j = 0, k = nh; j < nh; ++j, ++k)
    {
        double u, v;
        u = x[j];
        v = x[k];
        x[j] = u + v;
        x[k] = u - v;
    }
    c = 1.0;
    s = 0.0;
    jmax = nq + nh;
    for (j = nh + 1, k = n - 1; j < jmax; ++j, --k)
    {
        double u, v, tmp;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = u * c + v * s;
        x[k] = u * s - v * c;
    }
    nbranch <<= 1;
#if (_OPENMP >= 200203)
    if (nbranch > omp_get_max_threads())
    {
        fht_dif_rec(x, nh, nbranch);
        fht_dif_rec(x + nh, nh, nbranch);
    }
    else
    {
#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            fht_dif_rec(x, nh, nbranch);
#pragma omp section
            fht_dif_rec(x + nh, nh, nbranch);
        }
    }
#else
    fht_dif_rec(x, nh, nbranch);
    fht_dif_rec(x + nh, nh, nbranch);
#endif
#endif /* SIMPLE_FFT__BUILD__FHT_DIF_REC */
    return;
}
