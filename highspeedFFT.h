/**
 * @file mixed_radix_guide.c
 * @brief A comprehensive guide to the Mixed-Radix FFT Algorithm
 * 
 * This guide explains how the mixed-radix Fast Fourier Transform (FFT) works, 
 * based on the implementation in `highspeedFFT.h`. It covers splitting signals 
 * into smaller parts, using twiddle factors, and combining them into frequency 
 * buckets—all in the context of the provided code.
 * 
 * @author Tugbars Heptaskin
 * @date March 06, 2025
 */

/**
 * @mainpage Understanding the Mixed-Radix FFT Algorithm
 * 
 * The Fast Fourier Transform (FFT) transforms a signal into its frequency components efficiently.
 * The **mixed-radix algorithm** is a versatile FFT approach that handles various signal sizes by 
 * splitting them using different radixes (e.g., 2, 3, 4). This guide walks through the process, 
 * focusing on the `mixed_radix_dit_rec` function from the provided code.
 */

/**
 * @section goal_sec What’s the Goal?
 * 
 * The FFT converts a signal (a list of numbers) into **frequency buckets**—complex numbers 
 * representing the strength and phase of frequencies. The mixed-radix method achieves this by:
 * - Splitting the signal into smaller chunks using multiple radixes.
 * - Combining those chunks into frequency buckets with **twiddle factors**.
 * 
 * Why? It reduces computation from \( O(N^2) \) to \( O(N \log N) \), enabling real-time signal processing.
 */

/**
 * @section concepts_sec Key Concepts
 */

/**
 * @subsection signals_subsec Signals and Frequency Buckets
 * - **Input**: A signal of length \( N \), e.g., \( [1, 2, 3, 4] \) for \( N = 4 \).
 * - **Output**: \( N \) complex numbers (buckets), each tied to a frequency:
 *   - Bucket 0: 0 Hz (DC/constant term).
 *   - Bucket 1: \( f_s/N \) Hz (1 cycle), up to \( N/2 \) (highest frequency).
 *   - Each bucket has magnitude (strength) and phase (timing).
 */

/**
 * @subsection twiddles_subsec Twiddle Factors
 * - **Definition**: Complex numbers \( W_N^k = e^{-2\pi i k / N} = \cos(2\pi k / N) - i \sin(2\pi k / N) \).
 * - **Role**: Rotate signal pieces to sort them into frequency buckets.
 * - **In Code**: Precomputed in tables, e.g.:
 *   @code{.c}
 *   static const complex_t twiddle_radix4[] = {{1.0, 0.0}, {0.0, -1.0}, {-1.0, 0.0}, {0.0, 1.0}};
 *   @endcode
 *   - These are \( W_4^0 = 1 \), \( W_4^1 = -i \), etc., for \( N = 4 \).
 */

/**
 * @subsection radixes_subsec Radixes
 * - **Meaning**: Numbers like 2, 3, 4, 5, 7, 8 that define split sizes.
 * - **Mixed-Radix**: Uses multiple radixes (e.g., \( 12 = 4 \times 3 \)) instead of just powers of 2.
 */

/**
 * @section algorithm_sec How the Mixed-Radix Algorithm Works
 * 
 * The mixed-radix FFT uses a **decimation-in-time (DIT)** approach in `mixed_radix_dit_rec`. 
 * It splits the signal recursively, processes small pieces, and combines them into frequency buckets.
 * Let’s explore with \( N = 4 \), signal \( [1, 2, 3, 4] \), radix-2.
 */

/**
 * @subsection step1_subsec Step 1: Initialize the FFT Object
 * - **Function**: `fft_init(signal_length, transform_direction)`
 * - **What It Does**:
 *   - Checks if \( N \) factors into supported primes (2, 3, 4, 5, 7, 8, etc.) via `dividebyN`.
 *   - For \( N = 4 \), factors are \( [2, 2] \) (radix-2 twice).
 *   - Allocates `fft_object` and computes twiddles:
 *     @code{.c}
 *     fft_config->lf = factors(signal_length, fft_config->factors); // [2, 2]
 *     longvectorN(fft_config->twiddle, signal_length, fft_config->factors, fft_config->lf);
 *     @endcode
 *   - Sets `lt = 0` for mixed-radix.
 * - **Twiddles**: For \( N = 4 \), preps \( W_4^0 = 1 \), \( W_4^1 = -i \), stored in `twiddle`.
 */

/**
 * @subsection step2_subsec Step 2: Splitting the Signal
 * - **Function**: `mixed_radix_dit_rec` (recursive)
 * - **How It Works**:
 *   - Splits based on first radix (2):
 *     - Evens: \( [1, 3] \) (indices 0, 2).
 *     - Odds: \( [2, 4] \) (indices 1, 3).
 *   - Recursively splits:
 *     - Evens: \( [1] \), \( [3] \).
 *     - Odds: \( [2] \), \( [4] \).
 *   - **No Twiddles**: Splitting is index-based, using `stride`:
 *     @code{.c}
 *     mixed_radix_dit_rec(output_buffer, input_buffer, fft_obj, transform_sign, 2, 2, factor_index + 1); // Evens
 *     mixed_radix_dit_rec(output_buffer + 2, input_buffer + 1, fft_obj, transform_sign, 2, 2, factor_index + 1); // Odds
 *     @endcode
 */

/**
 * @subsection step3_subsec Step 3: Process Small Pieces (Butterflies)
 * - **When**: Hits smallest size (e.g., \( N = 2 \) or 1).
 * - **What Happens**:
 *   - For \( N = 1 \): Copy data:
 *     @code{.c}
 *     output_buffer[0].re = input_buffer[0].re; // e.g., 1
 *     @endcode
 *   - For \( N = 2 \) (radix-2 butterfly):
 *     - Evens \( [1, 3] \):
 *       - \( X_e[0] = 1 + 3 = 4 \).
 *       - \( X_e[1] = 1 - 3 = -2 \).
 *     - Odds \( [2, 4] \):
 *       - \( X_o[0] = 2 + 4 = 6 \).
 *       - \( X_o[1] = 2 - 4 = -2 \).
 *     @code{.c}
 *     output_buffer[0].re = tau1r + output_buffer[1].re; // 4
 *     output_buffer[1].re = tau1r - output_buffer[1].re; // -2
 *     @endcode
 * - **Twiddles**: Simple (\( W_2^0 = 1 \), \( W_2^1 = -1 \)), often implicit in add/subtract.
 */

/**
 * @subsection step4_subsec Step 4: Combining into Frequency Buckets
 * - **Function**: `mixed_radix_dit_rec` (combining)
 * - **How It Works**:
 *   - Combine \( [4, -2] \) (evens) and \( [6, -2] \) (odds) into \( N = 4 \):
 *     - Twiddles: \( W_4^0 = 1 \), \( W_4^1 = -i \) from `fft_obj->twiddle`.
 *   - Butterfly:
 *     - \( X[0] = X_e[0] + W_4^0 X_o[0] = 4 + 1 \cdot 6 = 10 \).
 *     - \( X[1] = X_e[1] + W_4^1 X_o[1] = -2 + (-i) \cdot (-2) = -2 + 2i \).
 *     - \( X[2] = X_e[0] - W_4^0 X_o[0] = 4 - 1 \cdot 6 = -2 \).
 *     - \( X[3] = X_e[1] - W_4^1 X_o[1] = -2 - (-i) \cdot (-2) = -2 - 2i \).
 *     @code{.c}
 *     tau2r = output_buffer[2].re * twiddle_real - output_buffer[2].im * twiddle_imag; // (-i) * (-2)
 *     output_buffer[0].re = tau1r + tau2r; // 10
 *     output_buffer[1].re = tau1r - tau2r; // -2 - 2i
 *     @endcode
 * - **Twiddles’ Role**:
 *   - **Rotation**: \( W_4^1 = -i \) rotates \( -2 \) to \( 2i \) for frequency 1.
 *   - **Separation**: Add vs. subtract sorts into buckets (0 Hz vs. 2 Hz).
 */

/**
 * @subsection step5_subsec Step 5: Output Frequency Buckets
 * - **Result**: \( [10, -2 + 2i, -2, -2 - 2i] \).
 * - **Buckets**:
 *   - \( X[0] = 10 \): 0 Hz (DC).
 *   - \( X[1] = -2 + 2i \): 1 cycle.
 *   - \( X[2] = -2 \): 2 cycles.
 *   - \( X[3] = -2 - 2i \): 3 cycles (mirror for real signals).
 * - **Execution**: `fft_exec` drives it:
 *   @code{.c}
 *   fft_exec(fft_obj, input_data, output_data);
 *   @endcode
 */

/**
 * @section flexibility_sec Mixed-Radix Flexibility
 * - **Multiple Radixes**: For \( N = 12 \):
 *   - Factors: \( [4, 3] \) (via `factors`).
 *   - Split: \( 12 \to 3 \) groups of 4, then 4 into 1s.
 *   - Combine: Radix-4, then radix-3, using `twiddle_radix4`, `twiddle_radix3`.
 * - **Supported Radixes**: Defined in:
 *   @code{.c}
 *   static const int primes[] = {2, 3, 4, 5, 7, 8, 11, 13, 17, 23, 29, 31, 37, 41, 43, 47, 53};
 *   @endcode
 */

/**
 * @section why_sec Why Split and Combine?
 * - **Speed**:
 *   - DFT: \( N^2 = 16 \) operations for \( N = 4 \).
 *   - FFT: \( N \log N = 4 \cdot 2 = 8 \) operations (2 stages).
 * - **Twiddles Enable It**: Rotate pieces to build frequencies, not the original signal.
 */

/**
 * @section summary_sec Summary
 * 1. **Setup**: `fft_init` factors \( N \) and preps twiddles.
 * 2. **Split**: Recursively divide into radix-sized chunks.
 * 3. **Process**: Butterflies on small pieces.
 * 4. **Combine**: Twiddles rotate and merge into buckets.
 * 5. **Output**: Frequency buckets reveal the signal’s spectrum.
 * 
 * The mixed-radix FFT is like chopping ingredients, mixing with spices (twiddles), 
 * and sorting into flavor jars (buckets)—fast and flexible.
 */


#ifndef HSFFT_H_
#define HSFFT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PI2 6.28318530717958647692528676655900577

#ifndef fft_type
#define fft_type double
#endif


typedef struct fft_t {
  fft_type re;
  fft_type im;
} fft_data;
/*
#define SADD(a,b) ((a)+(b))

#define SSUB(a,b) ((a)+(b))

#define SMUL(a,b) ((a)*(b))
*/

typedef struct fft_set* fft_object;

fft_object fft_init(int N, int sgn);

struct fft_set{
	int N;
	int sgn;
	int factors[64];
	int lf;
	int lt;
	fft_data twiddle[1];
};

void fft_exec(fft_object obj,fft_data *inp,fft_data *oup);

int divideby(int M,int d);

int dividebyN(int N);

//void arrrev(int M, int* arr);

int factors(int M, int* arr);

void twiddle(fft_data *sig,int N, int radix);

void longvectorN(fft_data *sig,int N, int *array, int M);

void free_fft(fft_object object);

#ifdef __cplusplus
}
#endif




#endif /* HSFFT_H_ */
