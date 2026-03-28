
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "bignum.h"

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;

typedef int64_t i64;
typedef int32_t i32;
typedef int16_t i16;

u64 calculate_ones(bn_uint* num) {
	u64 result = 0;
	for (u64 i = 0; i < num->nbits; i++) {
		if (bn_uint_get(num, i)) result++;
	}
	return result;
}

u64 bin_len(u64 num) {
	if (num == 0)
		return 0;

	u64 result = 0;
	for (i32 i = 0; i < 64; i++) {
		if (num & 1U << i)
			result = i + 1;
	}
	return result;
}

void print_binary(u64 n, u64 len, const char* text) {
	char buf[len + 1];
	memset(buf, '0', len + 1);
	buf[len] = '\0';
	for (i32 i = 0; i < len; i++) {
		buf[i] = ((n >> (len - i - 1)) & 1U) + '0';

	}
	if (text == NULL) {
		printf("%s\n", buf);
		return;
	}
	printf(text, buf);
}

i32 main(int argc, char* argv[]){

	if (argc < 2) {
		fprintf(stderr, "Usage: %s <graph size> <number of threads>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	u64 nbits = strtoul(argv[1], 0, 0);

	u64 nthreads;
	if (argc == 3) {
		nthreads = strtol(argv[2], 0, 0);
	}
	else {
		nthreads = omp_get_max_threads();
	}

	if (nthreads > omp_get_max_threads()) {
		fprintf(stderr, "Error: you do not possess that many threads sir!\n");
		exit(EXIT_FAILURE);
	}

	omp_set_num_threads((int)nthreads);

	u64 prefix_len = bin_len(nthreads - 1);
	u64 shift = nbits - prefix_len;
	u64 iter = 1UL << shift;

	#pragma omp parallel default(none) shared(shift, nbits)
	{ 
		u64 num = omp_get_thread_num();

		num <<= shift;

		u64 num_max = num | ((1 << shift) - 1);

		char msg_init[100];
		sprintf(msg_init, "Thread no. %d init: %%s\n", omp_get_thread_num());
		print_binary(num, nbits, msg_init);

		#pragma omp barrier
		#pragma omp single
		{
			printf("--------------------------------------------------\n");
		}

		char msg_16ones[100];
		sprintf(msg_16ones, "Thread no. %d 16 ones found: %%s\n", omp_get_thread_num());
		print_binary(num, nbits, msg_16ones);
		while (num < num_max) {
			if (calculate_ones(num) == 16) {
				print_binary(num, nbits, msg_16ones);
			}
			num++;
		}
	}
	return EXIT_SUCCESS;
}
