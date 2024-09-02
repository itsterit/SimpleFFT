#include "../../../SimpleFFT.h"
#include "../conversion.h"

/***********************************************************
 * Purpose:
 *   Given a real sequence, initially stored in the complex
 *   array z, this function computes the corresponding
 *   analytic sequence.  The imaginary part is the Hilbert
 *   transform of the real part.
 *
 * Arguments:
 *   double *z       - array of 2n doubles, representing n
 *                     complex numbers
 *   unsigned long n - dimension of z, must be a power of 2
 ************************************************************/

int hilbert_transform(hilbert_transform_type *z, len_type len)
{
#if SIMPLE_FFT__BUILD__HILBERT_TRANSFORM
	if (len > 2 && is_pow_of_two(len) && z != NULL)
	{
		double x;
		uint32_t i, n2;

		n2 = (uint32_t)len << 1;
		/*
		 * Compute the (bit-reversed) Fourier transform of z.
		 */
		forward_fft_dif(z, len);
		/*
		 * Form the transform of the analytic sequence by zeroing the
		 * transform for negative time, except for the (N/2)th.
		 * element.  Since z is now in bit-reversed order, this means
		 * zeroing every other complex element.  The array indices of
		 * the elements to be zeroed are 6,7,10,11...etc. (The real
		 * and imaginary parts of the (N/2)th element are in z[2] and
		 * z[3], respectively.)
		 */
		for (i = 6; i < n2; i += 4)
		{
			z[i] = 0.;
			z[i + 1] = 0.;
		}
		/*
		 * The 0th and (N/2)th elements get multiplied by 0.5.  Test
		 * for the trivial 1-point transform, just in case.
		 */
		z[0] *= 0.5;
		z[1] *= 0.5;
		if (len > 1)
		{
			z[2] *= 0.5;
			z[3] *= 0.5;
		}
		/*
		 * Compute the inverse transform.
		 */
		// inverse_fft_dit_type(z, len);
		/*
		 * Normalize the array.  The factor of 2 is left over from
		 * forming the transform in the time domain.
		 */
		x = 2. / (double)len;
		for (i = 0; i < n2; ++i)
			z[i] *= x;

		return 1;
	}
#endif /* SIMPLE_FFT__BUILD__HILBERT_TRANSFORM */
	return 0;
}

/***********************************************************
 * Purpose:
 *   Carry out a reverse binary bit-reversal permutation of
 *   the elements of a complex array z, with length len, where
 *   len is a power of 2.
 *
 * Arguments:
 *   double *z       - array of 2n doubles, representing len
 *                     complex numbers
 *   unsigned long len - dimension of z, must be a power of 2
 ************************************************************/

int bit_reverse(bit_reverse_type *z, len_type len)
{
#if SIMPLE_FFT__BUILD__BIT_REVERSE
	if (len > 2 && is_pow_of_two(len) && z != NULL)
	{
		uint32_t i;
		unsigned int ldn = 0;
		unsigned int rshift;

		i = (uint32_t)len;
		while (i >>= 1)
			++ldn;
		rshift = 8 * (unsigned int)sizeof(unsigned long) - ldn;
		for (i = 0; i < len; ++i)
		{
			unsigned long r;
#if (ULONG_MAX == 0xffffffff)
			r = ((i & 0x55555555) << 1) | ((i & ~0x55555555) >> 1);
			r = ((r & 0x33333333) << 2) | ((r & ~0x33333333) >> 2);
			r = ((r & 0x0f0f0f0f) << 4) | ((r & ~0x0f0f0f0f) >> 4);
			r = ((r & 0x00ff00ff) << 8) | ((r & ~0x00ff00ff) >> 8);
			r = (r << 16) | (r >> 16);
#elif (ULONG_MAX == 0xffffffffffffffff)
			r = ((i & 0x5555555555555555) << 1) | ((i & ~0x5555555555555555) >> 1);
			r = ((r & 0x3333333333333333) << 2) | ((r & ~0x3333333333333333) >> 2);
			r = ((r & 0x0f0f0f0f0f0f0f0f) << 4) | ((r & ~0x0f0f0f0f0f0f0f0f) >> 4);
			r = ((r & 0x00ff00ff00ff00ff) << 8) | ((r & ~0x00ff00ff00ff00ff) >> 8);
			r = ((r & 0x0000ffff0000ffff) << 16) |
				((r & ~0x0000ffff0000ffff) >> 16);
			r = (r << 32) | (r >> 32);
#endif
			r >>= rshift;
			if (r > i)
			{
				double tmp;
				unsigned long i2;
				i2 = i << 1;
				r <<= 1;
				tmp = z[i2];
				z[i2] = z[r];
				z[r] = tmp;
				tmp = z[i2 + 1];
				z[i2 + 1] = z[r + 1];
				z[r + 1] = tmp;
			}
		}

		return 1;
	}
#endif /* SIMPLE_FFT__BUILD__BIT_REVERSE */
	return 0;
}

/***********************************************************
 * Purpose:
 *   Carry out a reverse binary bit-reversal permutation of
 *   the elements of a real array x, with length len, where
 *   len is a power of 2.
 *
 * Arguments:
 *   double *x       - array of 2n doubles, representing len
 *                     real numbers
 *   unsigned long len - dimension of x, must be a power of 2
 ************************************************************/

int bit_reverse_real(bit_reverse_real_type *x, len_type len)
{
#if SIMPLE_FFT__BUILD__BIT_REVERSE_REAL
	if (len > 2 && is_pow_of_two(len) && x != NULL)
	{
		volatile uint32_t target = 0;
		for (uint32_t position = 0; position < len; position++)
		{
			if (target > position)
			{
				const bit_reverse_real_type temp = x[target];
				x[target] = x[position];
				x[position] = temp;
			}
			volatile uint32_t mask = (uint32_t)len;
			while (target & (mask >>= 1))
				target &= ~mask;
			target |= mask;
		}
		return 1;
	}
#endif /* SIMPLE_FFT__BUILD__BIT_REVERSE_REAL */
	return 0;
}
