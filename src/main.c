#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "usfio.h"

#define STEPS_PRINT_RATE (1<<18)
#define NPARALLEL_CHUNKS 12
#define NWORDS 50000000LU
#define FARPAD 1000000LU

#define MSB64 9223372036854775808LU
#define LSB64 1LU

uint64_t asbyte(char *s);

int main(int args, char **argv) {
	/* The goal is to find a number which has the highest steps in the collatz conjecture while not
	 * having more than 99% of its steps be divisions by two.
	 *
	 * Since we can anyways multiply the final result by 2 until there are 99 division operations per
	 * 3x+1 operations, we cannot cheese our way out of finding the biggest non-trivial collatz number */
	uint64_t steps, divsteps, mulsteps;
	steps = divsteps = mulsteps = 0;

	/* Array representing the very big number using machine-words */
	uint64_t *narr = NULL;

	/* To avoid actually shifting by 2 when dividing or as part of multiplying by 3, we use padding to
	 * offset the pointer to the first digit of our number. */
#define SCRATCHPAD (FARPAD * 8 + NWORDS * 8)
	if (posix_memalign((void **) &narr, 64, SCRATCHPAD)) {
		printf("posix_memalign failure, aborting.\n");
		exit(1);
	}

	/* wordoffset points to which machine word the LSB is currently pointing to */
	/* bitoffset masks with the bit which is the LSB */
	uint64_t *wordoffset = narr;
	uint64_t bitmask = 1LU;

	uint64_t *lastword; /* Point to last word and last bit offset */
	uint64_t bitsize; /* Size of number in bits */

	/* Initialize the number to test */
	if (args > 1) {
		char *num, *c;
		uint64_t flen;

		printf("Reading number from file %s\n", argv[1]);
		num = usf_ftos(argv[1], "r", &flen);
		num[--flen] = '\0'; /* Remove \n */

		unsigned int extra = flen % 64;
		bitsize = flen;

		if (bitsize / 64 > NWORDS) {
			printf("Not enough memory (%lu bytes allocated for number representation)!\n", NWORDS * 8);
			exit(1);
		}

		unsigned int i;
		for (i = 0, c = num + flen; c > num + extra; c -= 64, i++)
			wordoffset[i] = asbyte(c - 64);

		/* Last block */
		if (extra)
			wordoffset[i++] = asbyte(c - flen % 64) >> (64 - extra);

		lastword = wordoffset + i;

		free(num);
	} else {
		memset(wordoffset, 255, SCRATCHPAD - FARPAD * 8);
		bitsize = (SCRATCHPAD - FARPAD * 8) * 8;

		lastword = narr + (SCRATCHPAD - FARPAD) / sizeof(uint64_t);
	}

	struct timeval start, end;
	double elapsed;
	gettimeofday(&start, NULL);

	/* Main collatz loop */
	uint64_t *wordptr, oldval, twice;
#define MAX_CVALUE 4052555153018976267LU
#define MAX_KVALUE 6148914691236517204LU
	uint64_t cvalue = 1, kvalue = 0, laststep = 0;
	__uint128_t val128;

	printf("Starting computation...\n");
	while (bitsize != 1) {
		if ((*wordoffset & bitmask) == 0) {
			/* Is even */
			if (bitmask == MSB64) {
				/* Actualize next word only */
				val128 = (__uint128_t) *++wordoffset * cvalue + kvalue;
				*wordoffset = (uint64_t) val128;

				kvalue = (uint64_t) (val128 >> 64);
				bitmask = 1;
			} else bitmask <<= 1;
			steps++; divsteps++;
		} else {
			/* Using deferred multiplication for the 3x+1 step:
			 *
			 * --> Deferred multiplication:
			 * Instead of iterating over the whole number to execute 3x+1, which is extremely slow due to cache
			 * misses and the copious memory accessses, we can defer the evaluation of this step until doing a
			 * single combined pass would overflow the carry (64 bits) between words.
			 *
			 * We can keep lazy-evaluating the multiplication (only evaluating the next word when we get to it)
			 * up until the multiplicator (or kvalue, which is the result of applying 3x+1 to the starting +1
			 * value s times where s is the number of lazy steps) goes above 64 bits.
			 *
			 * The maximum svalue is the biggest one whose kvalue and cvalue are 64 bit.
			 *
			 * Furthermore, since we know that actualizing a lazy evaluation cannot increment the number of bits
			 * by more than 64 (mult by a 64 bit number), we can backshift the entire number if we detect the
			 * lastword as being one smaller than the last allocated word */

			if (cvalue <= MAX_CVALUE && kvalue <= MAX_KVALUE) {
				/* Lazy evaluation */
				steps++; mulsteps++;

				cvalue = (cvalue << 1) + cvalue;

				oldval = *wordoffset;
				twice = *wordoffset << 1;
				*wordoffset = twice + *wordoffset + bitmask;

				kvalue = (kvalue << 1) + kvalue + ((twice < oldval) + (*wordoffset <= twice));

				/* Hail mary for the last word to save itself by jumping to the next */
				if (wordoffset + 1 == lastword && kvalue) goto actualize;
			} else {
actualize:
				/* Actualize to reset cvalue and kvalue then continue as normal */
				for (wordptr = wordoffset + 1; wordptr < lastword; wordptr++) {
					val128 = (__uint128_t) *wordptr * cvalue + kvalue;
					*wordptr = (uint64_t) val128;

					kvalue = (uint64_t) (val128 >> 64);
				}

				if (kvalue) *lastword++ = kvalue; /* Carry out */

				cvalue = 1;
				kvalue = 0;
				bitsize = (lastword - narr - 1) * 64 +
					64 - __builtin_clzll(*(lastword - 1)) - __builtin_ctzll(bitmask);

				if (lastword >= narr + SCRATCHPAD / sizeof(uint64_t) - 2) {
					/* Shift rolling window back to start */
					memmove(narr, wordoffset, (lastword - wordoffset) * sizeof(uint64_t));
					lastword -= wordoffset - narr;
					wordoffset = narr;
				}

				/* Print here to avoid recalculating each time */
				if (steps >= laststep + STEPS_PRINT_RATE) {
					printf("Step %lu has %lu bits. (div/mul %lu %lu)\n", steps, bitsize, divsteps, mulsteps);
					laststep = steps;
				}
			}
		}

		if (wordoffset + 1 == lastword)
			bitsize = 64 - __builtin_clzll(*wordoffset) - __builtin_ctzll(bitmask);
	}

	gettimeofday(&end, NULL);
	elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;

	printf("Finished, took %lu steps and %lf seconds, with step ratios (div/mul) of %lu and %lu.\n", steps, elapsed, divsteps, mulsteps);
}

uint64_t asbyte(char *s) {
	uint64_t v = 0;

	for (unsigned int i = 0; i < 64; i++)
		v += (s[i] == '0' ? 0LU : 1LU) << (63 - i);

	return v;
}
