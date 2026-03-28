
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;

typedef int64_t i64;
typedef int32_t i32;
typedef int16_t i16;

#define PRINT 1

u64 bin_len(u64 num) {
	u64 len = 0;
	while (num != 0) {
		num >>= 1;
		len++;
	}
	return len;
}
u64 isqrt(u64 x) {
	if (x == 0) return 0;

	u64 r = x;

	while (1) {
		u64 nr = (r + x / r) / 2;
		if (nr >= r) break;
		r = nr;
	}

	return r;
}

bool isPrime(u64 num) {
	if (num < 2) return false;
	if (num == 2) return true;
	if (num % 2 == 0) return false;
	for (u64 i = 3; i < isqrt(num) + 2; i += 2) {
		if (num % i == 0) return false;
	}
	return true;
}

bool isBinaryPalindrome(u64 num) {
	u64 len = bin_len(num);
	for (i32 i = 0; i < len / 2; i++) {
		bool left_bit = (num >> (len - 1 - i)) & 1ULL;
		bool right_bit = (num >> i) & 1ULL;
		if (left_bit != right_bit) return false;
	}
	return true;
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
	if (argc < 3) {
		fprintf(stderr, "Usage: %s <number of threads> <bits count>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	u64 nbits = strtoul(argv[2], NULL, 0);

	i32 nthreads = atoi(argv[1]);

	if (nthreads > omp_get_max_threads()) {
		fprintf(stderr, "Error: you do not possess that many threads sir!\n");
		exit(EXIT_FAILURE);
	}

	omp_set_num_threads(nthreads);

	u64 iter = 1UL << nbits;

	u64 palindrome_primes = 0;

	#pragma omp parallel default(none) shared(iter, nbits) reduction(+:palindrome_primes)
	{ 
		#if PRINT
			bool init = false;
			#pragma omp for
			for (u64 i = 0; i < iter; i++) {
				if (!init) {
					if (PRINT) {
						char msg_init[100];
						sprintf(msg_init, "Thread no. %d init: %%s\n", omp_get_thread_num());
						print_binary(i, nbits, msg_init);
					}
					init = !init;
				}
			}

			#pragma omp barrier
			#pragma omp single
			{
				if(PRINT)printf("--------------------------------------------------\n");
			}
		#endif

		#pragma omp for schedule(dynamic, 1024)
		for (u64 i = 0; i < iter; i++) {
			if (isBinaryPalindrome(i) && isPrime(i)) {
				palindrome_primes++;
				#if PRINT
					char msg_prime_palindrome[100];
					sprintf(msg_prime_palindrome, "Thread no. %d found %lu: %%s\n", omp_get_thread_num(), i);
					print_binary(i, nbits, msg_prime_palindrome);
				#endif
			}
		}
	}
	printf("Palindrome primes: %lu\n", palindrome_primes);
	return EXIT_SUCCESS;
}
