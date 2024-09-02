#include "../../../SimpleFFT.h"
#include "../conversion.h"

/***********************************************************
 *  void fft_dif_iter(double *z, unsigned long n)
 *
 * Purpose:
 *   Computes the forward discrete Fourier transform of a
 *   complex sequence z, using an iterative radix 2
 *   decimation-in-frequency FFT.  z[0..2n-1] is an array of
 *   n complex numbers, stored in the usual way with real
 *   elements in z[0,2,..2n-2] and imaginary elements in
 *   z[1,3,..2n-1].
 *
 *   Entering this function, z[] should be in normal order.
 *   On return, the FT is stored in bit-reversed order.
 *
 *   n must be a power of 2.
 *
 * Arguments:
 *   double *z       - array of 2n doubles, representing n
 *                     complex numbers
 *   unsigned long n - dimension of z, must be a power of 2
 ************************************************************/

void fft_dif_iter(f_fft_dif_type *z, len_type n)
{
#if SIMPLE_FFT__BUILD__FFT_DIF_ITER
    unsigned long i, n2;

    n2 = n << 1;
    for (i = n; i > 1; i >>= 1)
    {
        double a, b, c, s, t;
        unsigned long i2, j;
        i2 = i << 1;
        t = TWOPI / (double)i;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        c = 1.0;
        s = 0.0;
        for (j = 0; j < i; j += 2)
        {
            double tmp;
            unsigned long kr, kmax;
            kmax = n2 + j;
            for (kr = j; kr < kmax; kr += i2)
            {
                double ur, ui;
                unsigned long ki, mr, mi;
                ki = kr + 1;
                mr = kr + i;
                mi = mr + 1;
                ur = z[kr];
                ui = z[ki];
                z[kr] = ur + z[mr];
                z[ki] = ui + z[mi];
                ur -= z[mr];
                ui -= z[mi];
                z[mr] = ur * c - ui * s;
                z[mi] = ur * s + ui * c;
            }
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
        }
    }
#endif /* SIMPLE_FFT__BUILD__FFT_DIF_ITER */
    return;
}

/***********************************************************
 *  void fft_dif_iter_seq(double *z, unsigned long n)
 *
 * Purpose:
 *   Computes the forward discrete Fourier transform of a
 *   complex sequence z, using an iterative radix 2
 *   decimation-in-frequency FFT.  Compared with
 *   fft_dif_iter(), above, the inner loops have been swapped
 *   to obtain sequential memory access at the expense of
 *   more trig overhead.
 *
 *   z[0..2n-1] is an array of n complex numbers, stored in
 *   the usual way with real elements in z[0,2,..2n-2] and
 *   imaginary elements in z[1,3,..2n-1].
 *
 *   Entering this function, z[] should be in normal order.
 *   On return, the FT is stored in bit-reversed order.
 *
 *   n must be a power of 2.
 *
 * Arguments:
 *   double *z       - array of 2n doubles, representing n
 *                     complex numbers
 *   unsigned long n - dimension of z, must be a power of 2
 ************************************************************/

void fft_dif_iter_seq(f_fft_dif_type *z, len_type n)
{
#if SIMPLE_FFT__BUILD__FFT_DIF_ITER_SEQ
    unsigned long i, n2;

    n2 = n << 1;
    for (i = n; i > 1; i >>= 1)
    {
        double a, b, t;
        unsigned long i2, k;
        i2 = i << 1;
        t = TWOPI / (double)i;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (k = 0; k < n2; k += i2)
        {
            double c, s;
            unsigned long j;
            c = 1.0;
            s = 0.0;
            for (j = 0; j < i; j += 2)
            {
                double ur, ui, tmp;
                unsigned long kr, ki, mr, mi;
                kr = k + j;
                ki = kr + 1;
                mr = kr + i;
                mi = mr + 1;
                ur = z[kr];
                ui = z[ki];
                z[kr] = ur + z[mr];
                z[ki] = ui + z[mi];
                ur -= z[mr];
                ui -= z[mi];
                z[mr] = ur * c - ui * s;
                z[mi] = ur * s + ui * c;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
            }
        }
    }
#endif /* SIMPLE_FFT__BUILD__FFT_DIF_ITER_SEQ */
    return;
}

/***********************************************************
*  void fft_dif_rec(double *z, unsigned long n, int nbranch)
*
* Purpose:
*   Computes the forward discrete Fourier transform of a
*   complex sequence z, using a radix 2 decimation-in-
*   frequency FFT.  If the computation is small enough to
*   fit in cache, it is done iteratively.  Otherwise, it is
*   done recursively until the recursion descends to
*   cache-sized chunks.
*
*   z[0..2n-1] is an array of n complex numbers, stored in
*   the usual way with real elements in z[0,2,..2n-2] and
*   imaginary elements in z[1,3,..2n-1].
*
*   Entering this function, z[] should be in normal order.
*   On return, the FT is stored in bit-reversed order.
*
*   n must be a power of 2.

*   To support OpenMP parallelism, nbranch keeps track of
*   the number of active transforms at a given recursion
*   level. On the first call to this function, nbranch
*   should be 1.  It is then doubled for each recursion.
*
* Arguments:
*   double *z       - array of 2n doubles, representing n
*                     complex numbers
*   unsigned long n - dimension of z, must be a power of 2
*   int nbranch     - number of transforms at this recursion
*                     level
************************************************************/
void fft_dif_rec(f_fft_dif_type *z, len_type n, int nbranch)
{
#if SIMPLE_FFT__BUILD__FFT_DIF_REC
    double a, b, c, s, t;
    unsigned long nh, kr;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / (2 * sizeof(double))))
    {
        if (FFT_UNIT_STRIDE)
            fft_dif_iter_seq(z, n);
        else
            fft_dif_iter(z, n);
        return;
    }
    t = TWOPI / (double)n;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    c = 1.0;
    s = 0.0;
    for (kr = 0; kr < n; kr += 2)
    {
        double ur, ui, tmp;
        unsigned long ki, mr, mi;
        ki = kr + 1;
        mr = kr + n;
        mi = mr + 1;
        ur = z[kr];
        ui = z[ki];
        z[kr] = ur + z[mr];
        z[ki] = ui + z[mi];
        ur -= z[mr];
        ui -= z[mi];
        z[mr] = ur * c - ui * s;
        z[mi] = ur * s + ui * c;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
    }
    nh = n >> 1;
    nbranch <<= 1;
#if (_OPENMP >= 200203)
    if (nbranch > omp_get_max_threads())
    {
        fft_dif_rec(z, nh, nbranch);
        fft_dif_rec(z + n, nh, nbranch);
    }
    else
    {
#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            fft_dif_rec(z, nh, nbranch);
#pragma omp section
            fft_dif_rec(z + n, nh, nbranch);
        }
    }
#else
    fft_dif_rec(z, nh, nbranch);
    fft_dif_rec(z + n, nh, nbranch);
#endif
#endif /* SIMPLE_FFT__BUILD__FFT_DIF_REC */
    return;
}