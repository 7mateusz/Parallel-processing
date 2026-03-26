#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

typedef uint64_t u64;
typedef uint32_t u32;

typedef int64_t i64;
typedef int32_t i32;

typedef double f64;
typedef float f32;

static u64 num_steps = 1000000000;

i32 main() {
	i64 i;
	f64 x, pi, sum = 0;
	f64 step = 1.0 / num_steps;

	clock_t start = clock();

	#pragma omp parallel
	{
		f64 x;
		#pragma omp for reduction(+:sum)
			for (i = 0; i < num_steps; i++) {
			x = (i + 0.5) * step;
			sum += 4.0 / (1.0 + x*x);
		}
	}
	pi = step * sum;

	clock_t stop = clock();

	printf("PI = \t\t%.20lf\n", pi);
	printf("math.h: \t%.20lf\n", M_PI);

	printf("\nTime: %f", (f32)(stop - start)/CLOCKS_PER_SEC);
}

