#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;

typedef int64_t i64;
typedef int32_t i32;
typedef int16_t i16;

u32 bin_len(u32 num) {
	if (num == 0)
		return 0;

	u32 result = 0;
	for (i32 i = 0; i < 32; i++) {
		if (num & 1U << i)
			result = i + 1;
	}
	return result;
}

void print_binary(u32 n) {
	u32 len = bin_len(n);

	char buf[len + 1];
	buf[len] = '\0';
	for (i32 i = 0; i < len; i++) {
		buf[i] = ((n >> (len - i - 1)) & 1U) + '0';

	}
	printf("%s\n", buf);
}

i32 main(int argc, char* argv[]){

	u32 a[4] = { 1U << 31, 1U << 16, 1U, ~0U };
	for (i32 i = 0; i < 4; i++) {
		printf("%d: %u digs\n", i, bin_len(a[i]));
	}

	#pragma omp parallel default(none)
	{ 
		u32 a = (omp_get_thread_num() << 15);

		while (1) {
			print_binary(a);
			printf("\t tid %d\n", omp_get_thread_num());
		}

	}
	return EXIT_SUCCESS;
}
