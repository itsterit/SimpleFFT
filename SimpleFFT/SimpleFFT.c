#include "SimpleFFT.h"
#include "../SimpleFFT/core/conversion/conversion.h"

/***********************************************************
 * Purpose:
 *   Computes the forward discrete Fourier transform of a
 *   complex sequence z, using a radix 2
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

int forward_fft_dif(f_fft_dif_type *z, len_type len)
{
#if SIMPLE_FFT__BUILD__FFT_DIF_ITER && SIMPLE_FFT__BUILD__FFT_DIF_ITER_SEQ && SIMPLE_FFT__BUILD__FFT_DIF_REC

    if (len > 2 && is_pow_of_two(len) && z != NULL)
    {
#if USE_RECURSIVE_VERSION__FFT_DIF
        fft_dif_rec(z, len, 1);
#else
        fft_dif_iter(z, len);
#endif
        return 1;
    }

#endif
    return 0;
}

/***********************************************************
 * Purpose:
 *   Computes the inverse discrete Fourier transform of a
 *   complex sequence z, using a radix 2 decimation-in-time
 *   FFT.  z[0..2n-1] is an array of n complex numbers,
 *   stored in the usual way with real elements in
 *   z[0,2,..2n-2] and imaginary elements in z[1,3,..2n-1].
 *
 *   Entering this function, z[] should be in bit-reversed
 *   order.  The returned inverse FT is restored to normal
 *   order.
 *
 *   n must be a power of 2.
 *
 * Arguments:
 *   double *z       - array of 2n doubles, representing n
 *                     complex numbers
 *   unsigned long n - dimension of z, must be a power of 2
 ************************************************************/

int inverse_fft_dit(i_fft_dit_type_type *z, len_type len)
{
#if SIMPLE_FFT__BUILD__IFFT_DIT_ITER && SIMPLE_FFT__BUILD__IFFT_DIT_ITER_SEQ && SIMPLE_FFT__BUILD__IFFT_DIT_REC

    if (len > 2 && is_pow_of_two(len) && z != NULL)
    {
#if USE_RECURSIVE_VERSION__IFFT_DIT
        ifft_dit_rec(z, len, 1);
#else
        ifft_dit_iter(z, len);
#endif
        return 1;
    }

#endif
    return 0;
}

/***********************************************************
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using a radix 2 decimation-in-
 *   frequency FHT.  n must be a power of 2.  Entering this
 *   function, x[] is in normal order.  On return, x[]
 *   contains the Hartley transform, stored in bit-reversed
 *   order.
 *
 * Arguments:
 *   double *x       - array of n doubles, representing n
 *                     real numbers
 *   unsigned long n - dimension of x, must be a power of 2
 ************************************************************/

int fht_dif(fht_dif_type *r, len_type len)
{
#if SIMPLE_FFT__BUILD__FHT_DIF_ITER && SIMPLE_FFT__BUILD__FHT_DIF_ITER_SEQ && SIMPLE_FFT__BUILD__FHT_DIF_REC

    if (len > 2 && is_pow_of_two(len) && r != NULL)
    {
#if USE_RECURSIVE_VERSION__FHT_DIF
        fht_dif_rec(r, len, 1);
#else
        fht_dif_iter(r, len);
#endif
        return 1;
    }

#endif
    return 0;
}

/***********************************************************
 * Purpose:
 *   Computes the discrete Hartley transform of a real
 *   sequence x[0..n-1], using a radix 2 decimation-in-time
 *   FHT.  n must be a power of 2.  Entering this function,
 *   x[] must be in bit-reversed order.  On return, x[]
 *   contains the Hartley transform, returned to normal
 *   order.
 *
 * Arguments:
 *   double *x       - array of n doubles, representing n
 *                     real numbers
 *   unsigned long n - dimension of x, must be a power of 2
 ************************************************************/

int fht_dit(fht_dit_type *r, len_type len)
{
#if SIMPLE_FFT__BUILD__FHT_DIT_ITER && SIMPLE_FFT__BUILD__FHT_DIT_ITER_SEQ && SIMPLE_FFT__BUILD__FHT_DIT_REC

    if (len > 2 && is_pow_of_two(len) && r != NULL)
    {
#if USE_RECURSIVE_VERSION__FHT_DIT
        fht_dit_rec(r, len, 1);
#else
        fht_dit_iter(r, len);
#endif
        return 1;
    }

#endif
    return 0;
}

/***********************************************************
 * Purpose:
 *   Checking if a number is a power of two
 *   return true if the number is a power of two
 *
 * Arguments:
 *   double num      - number to check
 ************************************************************/

int is_pow_of_two(int64_t num)
{
    return !((num > 0) && (num & (num - 1)));
}