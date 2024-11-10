#include <stdio.h>
#include <algorithm>
#include <pthread.h>
#include <math.h>
#include <immintrin.h>
#include<bitset>
#include "CycleTimer.h"
#include "sqrt_ispc.h"

using namespace ispc;

extern void sqrtSerial(int N, float startGuess, float* values, float* output);
__m256 fabs_avx2(__m256 input) {
    __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
    __m256 result = _mm256_and_ps(input, abs_mask);
    return result;
}

void my_sqrt_ispc(int N, float startGuess, float* values, float* output) {
    __m256 kThreshold = _mm256_set1_ps(0.00001f);
    __m256 one = _mm256_set1_ps(1.0f);
    __m256 three = _mm256_set1_ps(3.0f);
    __m256 half = _mm256_set1_ps(0.5f);

    for (int i = 0; i < N; i += 8) {
        int remaining = N - i;
        int mask_int = (remaining >= 8) ? 0xFF : (1 << remaining) - 1;
        __m256i mask_vec = _mm256_set_epi32(
            (mask_int & 0x80) ? -1 : 0,
            (mask_int & 0x40) ? -1 : 0,
            (mask_int & 0x20) ? -1 : 0,
            (mask_int & 0x10) ? -1 : 0,
            (mask_int & 0x08) ? -1 : 0,
            (mask_int & 0x04) ? -1 : 0,
            (mask_int & 0x02) ? -1 : 0,
            (mask_int & 0x01) ? -1 : 0
        );
        __m256 m_guess = _mm256_set1_ps(startGuess);
        __m256 m_values = _mm256_maskload_ps(&values[i], mask_vec);
        __m256 error = _mm256_mul_ps(m_guess, m_guess);
        error = _mm256_mul_ps(error, m_values);
        error = _mm256_sub_ps(error, one);
        error = _mm256_and_ps(error, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF)));
        __m256 mask_threshold = _mm256_cmp_ps(error, kThreshold, _CMP_GT_OQ);
        unsigned int mask = _mm256_movemask_ps(mask_threshold);

        while (mask) {
            __m256 guess_cubed = _mm256_mul_ps(m_guess, m_guess);
            guess_cubed = _mm256_mul_ps(guess_cubed, m_guess);
            __m256 term = _mm256_mul_ps(m_values, guess_cubed);
            __m256 new_guess = _mm256_sub_ps(_mm256_mul_ps(three, m_guess), term);
            m_guess = _mm256_mul_ps(new_guess, half);
            error = _mm256_mul_ps(m_guess, m_guess);
            error = _mm256_mul_ps(error, m_values);
            error = _mm256_sub_ps(error, one);
            error = _mm256_and_ps(error, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF)));
            mask_threshold = _mm256_cmp_ps(error, kThreshold, _CMP_GT_OQ);
            mask = _mm256_movemask_ps(mask_threshold);
        }
        __m256 result = _mm256_mul_ps(m_values, m_guess);
        _mm256_maskstore_ps(&output[i], mask_vec, result);
    }
}
static void verifyResult(int N, float* result, float* gold) {
    for (int i=0; i<N; i++) {
        if (fabs(result[i] - gold[i]) > 1e-4) {
            printf("Error: [%d] Got %f expected %f\n", i, result[i], gold[i]);
        }
    }
}

int main() {

    const unsigned int N = 20 * 1000 * 1000;
    const float initialGuess = 1.0f;

    float* values = new float[N];
    float* output = new float[N];
    float* gold = new float[N];

    for (unsigned int i=0; i<N; i++)
    {
        // TODO: CS149 students.  Attempt to change the values in the
        // array here to meet the instructions in the handout: we want
        // to you generate best and worse-case speedups
        // starter code populates array with random input values

        // values[i] = .001f + 2.998f * static_cast<float>(rand()) / RAND_MAX; // baseline
        // if (i != N - 1)values[i] = 1.0f;  //  worst case
        // else values[i] = 2.999f;
     
        values[i] = 2.999f; //best case
    }

    // generate a gold version to check results
    for (unsigned int i=0; i<N; i++)
        gold[i] = sqrt(values[i]);

    //
    // And run the serial implementation 3 times, again reporting the
    // minimum time.
    //
    double minSerial = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrtSerial(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minSerial = std::min(minSerial, endTime - startTime);
    }

    printf("[sqrt serial]:\t\t[%.3f] ms\n", minSerial * 1000);

    verifyResult(N, output, gold);

    //
    // Compute the image using the ispc implementation; report the minimum
    // time of three runs.
    //
    double minISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrt_ispc(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minISPC = std::min(minISPC, endTime - startTime);
    }

    printf("[sqrt ispc]:\t\t[%.3f] ms\n", minISPC * 1000);

    verifyResult(N, output, gold);

    // Clear out the buffer
    for (unsigned int i = 0; i < N; ++i)
        output[i] = 0;


    // AVX2 version made by myself
    minISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        my_sqrt_ispc(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minISPC = std::min(minISPC, endTime - startTime);
    }

    printf("[my sqrt avx2]:\t\t[%.3f] ms\n", minISPC * 1000);

    verifyResult(N, output, gold);

    // Clear out the buffer
    for (unsigned int i = 0; i < N; ++i)
        output[i] = 0;
    //
    // Tasking version of the ISPC code
    //
    double minTaskISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrt_ispc_withtasks(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minTaskISPC = std::min(minTaskISPC, endTime - startTime);
    }

    printf("[sqrt task ispc]:\t[%.3f] ms\n", minTaskISPC * 1000);

    verifyResult(N, output, gold);

    printf("\t\t\t\t(%.2fx speedup from ISPC)\n", minSerial/minISPC);
    printf("\t\t\t\t(%.2fx speedup from task ISPC)\n", minSerial/minTaskISPC);

    delete [] values;
    delete [] output;
    delete [] gold;

    return 0;
}
