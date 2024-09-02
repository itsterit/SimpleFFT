#ifndef _SIMPLE_FHT_CONFIG_H_
#define _SIMPLE_FHT_CONFIG_H_
#include <stdint.h>

#define USE_BASE_TYPE_IN_ALL_TRANSFORM (1)
typedef int64_t len_type;
typedef int32_t simple_fft_base_type;

#define USE_RECURSIVE_VERSION__FFT_DIF (0)
#define USE_RECURSIVE_VERSION__IFFT_DIT (0)
#define USE_RECURSIVE_VERSION__FHT_DIF (0)
#define USE_RECURSIVE_VERSION__FHT_DIT (0)

#define SIMPLE_FFT__BUILD__HILBERT_TRANSFORM (0)
#define SIMPLE_FFT__BUILD__BIT_REVERSE (0)
#define SIMPLE_FFT__BUILD__BIT_REVERSE_REAL (1)

/**
 * forward_fft_dif
 */
#define SIMPLE_FFT__BUILD__FFT_DIF_ITER (0)
#define SIMPLE_FFT__BUILD__FFT_DIF_ITER_SEQ (0)
#define SIMPLE_FFT__BUILD__FFT_DIF_REC (0)

/**
 * inverse_fft_dit
 */
#define SIMPLE_FFT__BUILD__IFFT_DIT_ITER (0)
#define SIMPLE_FFT__BUILD__IFFT_DIT_ITER_SEQ (0)
#define SIMPLE_FFT__BUILD__IFFT_DIT_REC (0)

/**
 * fht_dif
 */
#define SIMPLE_FFT__BUILD__FHT_DIF_ITER (0)
#define SIMPLE_FFT__BUILD__FHT_DIF_ITER_SEQ (0)
#define SIMPLE_FFT__BUILD__FHT_DIF_REC (0)

/**
 * fht_dit
 */
#define SIMPLE_FFT__BUILD__FHT_DIT_ITER (1)
#define SIMPLE_FFT__BUILD__FHT_DIT_ITER_SEQ (1)
#define SIMPLE_FFT__BUILD__FHT_DIT_REC (1)

/**
 * @brief			Forward discrete Fourier transform
 *
 * @param[in]		z 	normal order.
 * @param[in]		len dimension of z, must be a power of 2.
 * @param[out]		z 	FT is stored in bit-reversed order.
 */
#if USE_BASE_TYPE_IN_ALL_TRANSFORM
typedef simple_fft_base_type f_fft_dif_type;
#else
typedef double f_fft_dif_type;
#endif
int forward_fft_dif(f_fft_dif_type *z, len_type len);

/**
 * @brief			Inverse discrete Fourier transform
 *
 * @param[in]		z 	bit-reversed order.
 * @param[in]		len dimension of z, must be a power of 2.
 * @param[out]		z 	inverse FT is restored to normal order.
 */
#if USE_BASE_TYPE_IN_ALL_TRANSFORM
typedef simple_fft_base_type i_fft_dit_type_type;
#else
typedef double i_fft_dit_type_type;
#endif
int inverse_fft_dit(i_fft_dit_type_type *z, len_type len);

/**
 * @brief			Discrete Hartley transform of a real sequence decimation-in-frequency
 *
 * @param[in]		r 	normal order.
 * @param[in]		len dimension of r, must be a power of 2.
 * @param[out]		r 	Hartley transform, stored in bit-reversed order.
 */
#if USE_BASE_TYPE_IN_ALL_TRANSFORM
typedef simple_fft_base_type fht_dif_type;
#else
typedef double fht_dif_type;
#endif
int fht_dif(fht_dif_type *r, len_type len);

/**
 * @brief			Discrete Hartley transform of a real sequence decimation-in-time
 *
 * @param[in]		r 	bit-reversed order
 * @param[in]		len dimension of r, must be a power of 2.
 * @param[out]		r 	Hartley transform, stored in normal order order.
 */
#if USE_BASE_TYPE_IN_ALL_TRANSFORM
typedef simple_fft_base_type fht_dit_type;
#else
typedef double fht_dit_type;
#endif
int fht_dit(fht_dit_type *r, len_type len);

/**
 * @brief			The imaginary part is the Hilbert transform of the real part.
 *
 * @param[in]		z 	complex array
 * @param[in]		len dimension of z, must be a power of 2.
 * @param[out]		z 	real sequence.
 */
#if USE_BASE_TYPE_IN_ALL_TRANSFORM
typedef simple_fft_base_type hilbert_transform_type;
#else
typedef double hilbert_transform_type;
#endif
int hilbert_transform(hilbert_transform_type *z, len_type len);

/**
 * @brief			Reverse binary bit-reversal permutation of
 *					the elements of a complex array.
 *
 * @param[in]		z 	complex array
 * @param[in]		len dimension of z, must be a power of 2.
 * @param[out]		z 	bit-revers sequence.
 */
#if USE_BASE_TYPE_IN_ALL_TRANSFORM
typedef simple_fft_base_type bit_reverse_type;
#else
typedef double bit_reverse_type;
#endif
int bit_reverse(bit_reverse_type *z, len_type len);

/**
 * @brief			Reverse binary bit-reversal permutation of
 *					the elements of a real array.
 *
 * @param[in]		x 	real array
 * @param[in]		len dimension of x, must be a power of 2.
 * @param[out]		x 	bit-revers sequence.
 */
#if USE_BASE_TYPE_IN_ALL_TRANSFORM
typedef simple_fft_base_type bit_reverse_real_type;
#else
typedef double bit_reverse_real_type;
#endif
int bit_reverse_real(bit_reverse_real_type *x, len_type len);

/**
 * @brief			Checking if a number is a power of two
 *
 * @param[in]		num number to check
 * @param[out]		int return true if the number is a power of two
 */
int is_pow_of_two(int64_t num);

/**
 * @brief			benchmarks sequence
 */
void benchmarks(len_type sampling_rate, len_type pnt_amt);

#endif /* _SIMPLE_FHT_CONFIG_H_ */