#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
 
int main(int argc, char* argv[]) {
 
// {200203:'2.0',200505:'2.5',200805:'3.0',201107:'3.1',201307:'4.0',201511:'4.5',201811:'5.0',202011:'5.1',202111:'5.2'}
#ifdef _OPENMP
  printf("_OPENMP: %d\n",_OPENMP);
  switch (_OPENMP) { 
   case 200203: printf("OMP version 2.0\n"); break;
   case 200505: printf("OMP version 2.5\n"); break;
   case 200805: printf("OMP version 3.0\n"); break;
   case 201107: printf("OMP version 3.1\n"); break;
   case 201307: printf("OMP version 4.0\n"); break;
   case 201511: printf("OMP version 4.5\n"); break;
 //case 201811: printf("OMP version 5.0\n"); break;
 //case 202011: printf("OMP version 5.1\n"); break;
 //case 202111: printf("OMP version 5.2\n"); break;
 //case 202411: printf("OMP version 6.0\n"); break;
   default: 
	 printf("OMP version ???");
  } 
#endif
return EXIT_SUCCESS;
}
