// http://www.fftw.org/pruned.html

#include <math.h>
#include <complex.h>
#include <fftw3.h>
{
int i, j;
const double TWOPI = 6.2831853071795864769252867665590057683943388;
fftw_complex in[N], out[N], twids[(K - 1) * (N / K - 1)];
fftw_plan plan;

/* plan N/K FFTs of size K */
plan = fftw_plan_many_dft(1, &K, N / K,
    in, NULL, N / K, 1,
    out, NULL, 1, K,
    FFTW_FORWARD, FFTW_ESTIMATE);

/* precompute twiddle factors (since we usually want more than one FFT) */
for (j = 1; j < N / K; ++j)
    for (i = 1; i < K; ++i)
        twids[(j - 1) * (K - 1) + (i - 1)] = cexp((I * FFTW_FORWARD * TWOPI / N) * (i * j));

...initialize in[N] data....

fftw_execute(plan);

/* set *only* first K outputs, in out[], to values for full size-N DFT: */
for (j = 1; j < N / K; ++j) {
    out[0] += out[j * K];
    for (i = 1; i < K; ++i)
        out[i] += out[i + j * K] * twids[(j - 1) * (K - 1) + (i - 1)];
}

fftw_destroy_plan(plan);
}