#include <stdio.h>


#include "sito.h"


int main(int argc, char *argv[])
{
  FILE *in, *out;
  char BUFFOR[1024];
  
  long double a[400];
  long double x[21];
  ta A;
  int  N; 
  long LICZ;
 if (argc < 3)
  { printf("Za malo parametrow\n");
    //exit(0);
    return 0;
  }
 else
  {
   /*
   printf("Plik zrodowy       :%s\n",argv[1]);
   printf("Plik przeznaczenia :%s\n",argv[2]); 
   */
   in = fopen(argv[1],"rb");
   if (in != NULL) 
    {
      /* printf("Dostep do pliku zrodlowego OK!\n"); */
      out = fopen(argv[2],"wb"); fclose(out); 
         LICZ = 0;
         while (1)
         {
          if (fgets(BUFFOR,1024,in)!= NULL)
          { LICZ++;
            BMKdecode(BUFFOR,&N,A);
            AToa(N,A,a);
            eigensymmatrix(N, a, 1, N,  x);
            /* printX(N,x); */
            if (isintegral(N, x)) 
            printSout(N, A, argv[2]);
            

            /*  
             if (connected(N,A))
             {
               if (integral(N, A)) { printSout(N, A, argv[2]); }
               dop(N,A);
               if (connected(N,A) && integral(N, A)) { printSout(N, A, argv[2]); }
             }  
            else 
             {  
               dop(N,A);
               if (integral(N, A)) { printSout(N, A, argv[2]); }
             }
             */
          } else break; 
          
         }
         /* printf("\nLicznik %ld\n", LICZ); */
       
      fclose(in);
    } else printf("Brak pliku zrodlowego\n");
   }
  return 0;
}

