#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "real.h"  // Assumes fft_real_object, fft_r2c_exec, fft_c2r_exec, etc., are provided

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

/**
 * @brief Computes the next power of 2 greater than or equal to n.
 * 
 * @param n Input number.
 * @return int The smallest power of 2 >= n.
 */
int next_power_of_two(int n) {
    if (n <= 0) return 1;  // Handle edge case
    return (int)pow(2, ceil(log2(n)));
}

/**
 * @brief Finds the optimal FFT length for convolution.
 * 
 * For linear convolution, pads to the next power of 2 >= min_length.
 * For circular convolution, pads to the next power of 2 >= max(N, L).
 * 
 * @param min_length Minimum required length (e.g., N + L - 1 for linear).
 * @param conv_type Convolution type: "linear" or "circular".
 * @param length1 Length of the first input signal.
 * @param length2 Length of the second input signal.
 * @return int Optimal padded length for FFT.
 */
int find_optimal_fft_length(int min_length, const char *conv_type, int length1, int length2) {
    if (strcmp(conv_type, "linear") == 0) {
        return next_power_of_two(min_length);
    } else if (strcmp(conv_type, "circular") == 0) {
        int max_length = MAX(length1, length2);
        return next_power_of_two(max_length);
    } else {
        fprintf(stderr, "Error: Invalid convolution type\n");
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Performs FFT-based convolution of two real-valued signals.
 * 
 * Supports both linear and circular convolution, with output types "full", "same", or "valid".
 * 
 * @param type Output type: "full" (full convolution), "same" (central portion matching input1), or "valid" (no padding effects).
 * @param conv_type Convolution type: "linear" (standard) or "circular" (periodic).
 * @param input1 First input signal array.
 * @param length1 Length of the first input signal.
 * @param input2 Second input signal array.
 * @param length2 Length of the second input signal.
 * @param output Array to store the convolution result.
 * @return int Length of the output array, or -1 on error.
 * 
 * @note For linear convolution, pads to the next power of 2 >= N + L - 1.
 *       For circular convolution, pads to the next power of 2 >= max(N, L).
 */
int fft_convolve(const char *type, const char *conv_type, fft_type *input1, int length1, 
                 fft_type *input2, int length2, fft_type *output) {
    // Input validation
    if (input1 == NULL || input2 == NULL || output == NULL || length1 <= 0 || length2 <= 0) {
        fprintf(stderr, "Error: Invalid inputs for fft_convolve\n");
        return -1;
    }

    // Determine convolution length based on type
    int conv_length;
    if (strcmp(conv_type, "linear") == 0) {
        conv_length = length1 + length2 - 1;
    } else if (strcmp(conv_type, "circular") == 0) {
        conv_length = MAX(length1, length2);
    } else {
        fprintf(stderr, "Error: Invalid convolution type. Use 'linear' or 'circular'.\n");
        return -1;
    }

    // Find optimal padded length for FFT
    int padded_length = find_optimal_fft_length(conv_length, conv_type, length1, length2);

    // Initialize forward and inverse real FFT objects
    fft_real_object fobj = fft_real_init(padded_length, 1);  // Forward FFT
    fft_real_object iobj = fft_real_init(padded_length, -1); // Inverse FFT
    if (!fobj || !iobj) {
        fprintf(stderr, "Error: Failed to initialize real FFT objects\n");
        if (fobj) free_real_fft(fobj);
        if (iobj) free_real_fft(iobj);
        return -1;
    }

    // Allocate memory for padded inputs and FFT outputs
    fft_type* padded_input1 = (fft_type*)calloc(padded_length, sizeof(fft_type));
    fft_type* padded_input2 = (fft_type*)calloc(padded_length, sizeof(fft_type));
    fft_data* complex_output1 = (fft_data*)malloc(sizeof(fft_data) * padded_length);
    fft_data* complex_output2 = (fft_data*)malloc(sizeof(fft_data) * padded_length);
    fft_data* complex_product = (fft_data*)malloc(sizeof(fft_data) * padded_length);
    fft_type* final_output = (fft_type*)malloc(sizeof(fft_type) * padded_length);

    if (!padded_input1 || !padded_input2 || !complex_output1 || !complex_output2 || !complex_product || !final_output) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free(padded_input1); free(padded_input2); free(complex_output1);
        free(complex_output2); free(complex_product); free(final_output);
        free_real_fft(fobj); free_real_fft(iobj);
        return -1;
    }

    // Copy inputs to padded arrays (zero-padded by calloc)
    memcpy(padded_input1, input1, length1 * sizeof(fft_type));
    memcpy(padded_input2, input2, length2 * sizeof(fft_type));

    // Perform forward real-to-complex FFT on both inputs
    fft_r2c_exec(fobj, padded_input1, complex_output1);
    fft_r2c_exec(fobj, padded_input2, complex_output2);

    // Multiply frequency-domain representations pointwise
    for (int i = 0; i < padded_length; i++) {
        complex_product[i].re = complex_output1[i].re * complex_output2[i].re - complex_output1[i].im * complex_output2[i].im;
        complex_product[i].im = complex_output1[i].re * complex_output2[i].im + complex_output1[i].im * complex_output2[i].re;
    }

    // Perform inverse complex-to-real FFT
    fft_c2r_exec(iobj, complex_product, final_output);

    // Scale the output
    for (int i = 0; i < padded_length; i++) {
        final_output[i] /= padded_length;
    }

    // Determine output length and starting point based on type and conv_type
    int start = 0, padding = 0;
    if (strcmp(conv_type, "linear") == 0) {
        if (strcmp(type, "full") == 0 || type == NULL) {
            start = 0;
            padding = conv_length;
        } else if (strcmp(type, "same") == 0) {
            int larger_base = MAX(length1, length2);
            start = (conv_length - larger_base) / 2;
            padding = larger_base;
        } else if (strcmp(type, "valid") == 0) {
            int smaller_length = MIN(length1, length2);
            start = smaller_length - 1;
            padding = (smaller_length == 0) ? 0 : (MAX(length1, length2) - smaller_length + 1);
        } else {
            fprintf(stderr, "Error: Invalid output type. Use 'full', 'same', or 'valid'.\n");
            padding = -1;
        }
    } else if (strcmp(conv_type, "circular") == 0) {
        // For circular convolution, output is typically the full padded length
        start = 0;
        padding = padded_length;
        // Note: 'same' or 'valid' may not directly apply; adjust based on use case
    }

    if (padding > 0) {
        // Copy the relevant portion to output
        memcpy(output, final_output + start, padding * sizeof(fft_type));
    }

    // Clean up memory
    free(padded_input1); free(padded_input2); free(complex_output1);
    free(complex_output2); free(complex_product); free(final_output);
    free_real_fft(fobj); free_real_fft(iobj);

    return padding;  // Return the length of the output
}