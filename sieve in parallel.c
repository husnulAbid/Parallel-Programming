/*
  Sieve of Erathostenes
  Assignment 3 - Parallel Programming: Abid and Alex
  Steps to run:
  $ module load GCC
  $ gcc -fopenmp sieve_try4.c -o sieve_try4 -lm
  $ srun sieve_try4 100000 --mem=1G -t 1:00:00
*/

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int main (int argc, char *argv[]) {

  #ifdef _OPENMP
  printf("OpenMp is enabled\n\n");
  #endif

  int i, k, t, N, nr_primes, lastprime=0, chunk_size=100;
  char *prime=NULL;
  const char unmarked = (char)0;
  const char marked = (char)1;
  clock_t start, stop;
  
  if (argc < 2) {
    printf("Usage:  %s N\n", argv[0]);
    exit(-1);
  }
  N = atoi(argv[1]);

  start= clock();   /* Start measuring time */

  /* Allocate marks */
  prime = malloc(N); 
  if (prime == NULL) {
    printf("Could not allocate %d chars of memory\n", N);
    exit(-1);
  }

  /* Mark 0 and 1 as not prime */
  prime[0] = marked;
  prime[1] = marked;

  /* First parallel region */
  #pragma omp parallel private(i) shared(N, prime)
  {

    /* Initialize all numbers to unmarked */
    #pragma omp for nowait
    for (i=2; i<N+1; i++) {
        prime[i] = unmarked;
    }
    
  } /* End of parallel region */

  /* Adapting chunk size for bigger values of N */
  if(N >= 10000000) {
        chunk_size = 10000;
  }

  if(N >= 100000) {
        chunk_size = 1000;
  }

  #pragma omp parallel for schedule(dynamic, chunk_size) // Parallel for with dynamic scheduling
  for (i=2; i <= (int)sqrt(N)+1; i++) {
    if (prime[i]==unmarked) {   // Next unmarked position
      t = i;  // Position i corresponds to the value t
      for (k=t*t; k<=N; k+=t) {
	      prime[k] = marked;  // Mark the multiples of i 
      }
    }
  }
  
  nr_primes = 0;  /* Remember to count 2 as a prime */
  /* Count the marked numbers */
  #pragma omp parallel for lastprivate(lastprime) reduction(+ : nr_primes)
  for (i=2; i<=N; i++) {
    if (prime[i]==unmarked) {
      lastprime = i;
      nr_primes++;
    }
  }
  
  stop = clock();
  printf("Time: %6.5f s\n", (float)(stop-start)/CLOCKS_PER_SEC);

  printf("\n%d primes smaller than or equal to %d\n", nr_primes, N);
  printf("The largest of these primes is %d\n", lastprime);
  printf("\nReady\n");

  exit(0);
}
