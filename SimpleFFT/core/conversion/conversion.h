#ifndef _CONVERSION_H_
#define _CONVERSION_H_
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if (_OPENMP >= 200203)
#include <omp.h>
#endif

#include "../param/am_sysdep.h"
#include "../param/math_const.h"

/*
 * Selects iterative transform types.
 */
#ifndef FFT_UNIT_STRIDE
#define FFT_UNIT_STRIDE 0
#endif
#ifndef FHT_UNIT_STRIDE
#define FHT_UNIT_STRIDE 0
#endif

void fft_dif_iter(f_fft_dif_type *, len_type);
void fft_dif_iter_seq(f_fft_dif_type *, len_type);
void fft_dif_rec(f_fft_dif_type *, len_type, int);

void ifft_dit_iter(i_fft_dit_type_type *, len_type);
void ifft_dit_iter_seq(i_fft_dit_type_type *, len_type);
void ifft_dit_rec(i_fft_dit_type_type *, len_type, int);

void fht_dif_iter(fht_dif_type *, len_type);
void fht_dif_iter_seq(fht_dif_type *, len_type);
void fht_dif_rec(fht_dif_type *, len_type, int);

void fht_dit_iter(fht_dit_type *, len_type);
void fht_dit_iter_seq(fht_dit_type *, len_type);
void fht_dit_rec(fht_dit_type *, len_type, int);
#endif /* _CONVERSION_H_ */