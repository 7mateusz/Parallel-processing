#ifndef BIGNUM_H
#define BIGNUM_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef int64_t i64;
typedef int32_t i32;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

typedef struct {
	u64* words;
	u64 cap;
	u64 used;
} bn_uint;

// ---------- lifetime ---------- 
void bn_init(bn_uint* bn);
void bn_init_u64(bn_uint* bn, u64 val);
void bn_init_bits(bn_uint* bn, u64 nbits);
void bn_free(bn_uint* bn);

// ---------- assignment ---------- 
void bn_reserve(bn_uint* bn, u64 nwords);
void bn_copy(bn_uint* dst, const bn_uint* src);
void bn_set_u64(bn_uint* bn, u64 val);
void bn_set_zero(bn_uint* bn);
void bn_swap(bn_uint* a, bn_uint* b);

// ---------- comparison ---------- 
int bn_cmp(const bn_uint* a, const bn_uint* b);
bool bn_eq(const bn_uint* a, const bn_uint* b);
bool bn_iszero(const bn_uint* a);
bool bn_is_odd(const bn_uint* a);
bool bn_is_even(const bn_uint* a);

// ---------- queries ---------- 
u64 bn_bitlen(const bn_uint* bn);
bool bn_getbit(const bn_uint* bn, u64 bit);
void bn_setbit(bn_uint* bn, u64 bit, bool val);

// ---------- arithmetic (bignum × bignum) ---------- 
u32 bn_add(bn_uint* result, const bn_uint* a, const bn_uint* b);
u32 bn_sub(bn_uint* result, const bn_uint* a, const bn_uint* b);
void bn_mul(bn_uint* result, const bn_uint* a, const bn_uint* b);
void bn_div(bn_uint* quot, bn_uint* rem, const bn_uint* a, const bn_uint* b);

// ---------- arithmetic (bignum × u64) ---------- 
u32 bn_add_u64(bn_uint* result, const bn_uint* a, u64 b);
u32 bn_sub_u64(bn_uint* result, const bn_uint* a, u64 b);
void bn_mul_u64(bn_uint* result, const bn_uint* a, u64 b);
u64 bn_div_u64(bn_uint* result, const bn_uint* a, u64 b);
u64 bn_mod_u64(const bn_uint* a, u64 b);

// ---------- bitwise ---------- 
void bn_and(bn_uint* result, const bn_uint* a, const bn_uint* b);
void bn_or (bn_uint* result, const bn_uint* a, const bn_uint* b);
void bn_xor(bn_uint* result, const bn_uint* a, const bn_uint* b);
void bn_not(bn_uint* result, const bn_uint* a);
void bn_shl(bn_uint* result, const bn_uint* a, u64 shift);
void bn_shr(bn_uint* result, const bn_uint* a, u64 shift);

// ---------- modular arithmetic ----------
void bn_mod(bn_uint* result, const bn_uint* a, const bn_uint* m);
void bn_addmod(bn_uint* result, const bn_uint* a, const bn_uint* b, const bn_uint* m);
void bn_submod(bn_uint* result, const bn_uint* a, const bn_uint* b, const bn_uint* m);
void bn_mulmod(bn_uint* result, const bn_uint* a, const bn_uint* b, const bn_uint* m);
void bn_powmod(bn_uint* result, const bn_uint* base, const bn_uint* exp, const bn_uint* m);
bool bn_invmod(bn_uint* result, const bn_uint* a, const bn_uint* m);

// ---------- number theory ----------
void bn_gcd(bn_uint* result, const bn_uint* a, const bn_uint* b);
void bn_lcm(bn_uint* result, const bn_uint* a, const bn_uint* b);

// ---------- string conversion ----------
char* bn_to_hex(const bn_uint* bn);
char* bn_to_dec(const bn_uint* bn);
char* bn_to_bin(const bn_uint* bn);
bool bn_from_hex(bn_uint* bn, const char* str);
bool bn_from_dec(bn_uint* bn, const char* str);

#endif // BIGNUM_H

#define BIGNUM_IMPLEMENTATION
#ifdef BIGNUM_IMPLEMENTATION
// NOLINTBEGIN(misc-definitions-in-headers)

// ---------- lifetime imp ----------

void bn_init(bn_uint* bn) {
	bn->words = NULL;
	bn->cap = 0;
	bn->used = 0;
}

void bn_init_u64(bn_uint* bn, u64 val) {
	bn->words = (u64*)malloc(sizeof(u64));
	if (bn->words == NULL) {
		fprintf(stderr, "BN Error: malloc failed!\n");
		exit(EXIT_FAILURE);
	}
	bn->words[0] = val;

	bn->cap = 1;
	bn->used = (val == 0) ? 0 : 1;
}

void bn_init_bits(bn_uint* bn, u64 nbits) {
	bn->words = (u64*)calloc((nbits + 63) / 64, sizeof(u64));
	if (bn->words == NULL) {
		fprintf(stderr, "BN Error: malloc failed!\n");
		exit(EXIT_FAILURE);
	}

	bn->cap = (nbits + 63) / 64;
	bn->used = 0;
}

void bn_free(bn_uint* bn) {
	free(bn->words);
	bn->words = NULL;
	bn->used = 0;
	bn->cap = 0;
}


// ---------- assignment ---------- 

void bn_reserve(bn_uint* bn, u64 nwords) {
	if (bn->cap >= nwords)
		return;

	bn->words = (u64*)realloc(bn->words, nwords * sizeof(u64));
	if (bn->words == NULL) {
        fprintf(stderr, "BN Error: realloc failed!\n");
        exit(EXIT_FAILURE);
    }
	memset(bn->words + bn->used, 0, (nwords - bn->used) * sizeof(u64));

	bn->cap = nwords;
}

void bn_copy(bn_uint* dst, const bn_uint* src) {
	
}

void bn_set_u64(bn_uint* bn, u64 val);
void bn_set_zero(bn_uint* bn);
void bn_swap(bn_uint* a, bn_uint* b);


// ---------- comparison imp ----------

int bn_cmp(const bn_uint* a, const bn_uint* b) {
	if (a->used > b->used) return 1;
	if (b->used > a->used) return -1;

	for (i64 i = a->used - 1; i >= 0; i--) {
		if (a->words[i] > b->words[i]) return 1;
		if (b->words[i] > a->words[i]) return -1;
	}

	return 0;
}

bool bn_eq(const bn_uint* a, const bn_uint* b) {
	if (a->used != b->used) return false;

	for (i64 i = a->used - 1; i >= 0; i--) {
		if (a->words[i] != b->words[i]) return false;
	}

	return true;
}

bool bn_iszero(const bn_uint* a) {
	if (a->used != 0) return false;
	else return true;
}

bool bn_is_odd(const bn_uint* a) {
	if (bn_iszero(a)) return false;
	if (a->words[0] & 1ULL) return true;
	return false;
}

bool bn_is_even(const bn_uint* a) {
	if (!bn_is_odd(a)) return true;
	return false;
}

// NOLINTEND(misc-definitions-in-headers)
#endif // ! BIGNUM_IMPLEMENTATION
