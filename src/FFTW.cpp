
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fftw3.h>
#include <ctime>
#include <time.h>

using namespace std;

int main()
{
	double thetime;
	clock_t start, end;

	fftwf_init_threads();
	fftwf_plan_with_nthreads(4);

	int N = 1048576;

	float *in;
	fftwf_complex *out;
	fftwf_plan my_plan;

	in = (float*)fftwf_malloc(sizeof(float) * N);
	out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (N/2 + 1));

	std::srand(std::time(0));
	for (int n = 0; n != N; n++)
	{
		in[n] = (float)(((float)std::rand()) / RAND_MAX);
	}

	my_plan = fftwf_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); // FFTW_MEASURE
	fftwf_execute(my_plan);
	fftwf_destroy_plan(my_plan);
	
	for (int k = 0; k < 8; k++)
	{
		cout << out[k][0] << " " << out[k][1] << " \n";
	}
	
	fftwf_free(in);
	fftwf_free(out);

	return 0;
}
