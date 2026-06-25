#include <stdio.h>
#include <math.h>
#include "sito.h" 

/*  procedura jest wzorowana na kodzie Pascalowym A.Marciniaka */
int eigensymmatrix( int n,
                     long double* a,
                     int k1, int k2,
                     long double* x
                     )
{
  int i,j,k,k3,k4,L,L1,z;
  long double lambda,eps,g,h,ma,mn,norm,s,t,u,w;
  int cond;
  long double d[20], e[20], e2[20], Lb[20];

  if ((1<=k1) && (k1<=k2) && (k2<=n))
   {
    i = 0;
    for (L=1;L<=n;L++) { i += L; d[L] = a[i]; }
    for (L=n;L>=2;L--)
     {
      i--; j = i; h = a[j]; s = 0;
      for (k=L-2;k>=1;k--) { i--; g = a[i]; s += g*g; }
      i--;
      if (s == 0) { e[L] = h; e2[L] = h*h; a[j] = 0.0; }
       else
        {
          s += h*h; e2[L] = s; g = sqrt(s); if (h>=0.0) g=-g;
          e[L] = g;
          s = 1.0 / (s-h*g);
          a[j] = h - g; h = 0.0; L1 = L - 1; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             k4 = k3; g = 0;
             for (k=1;k<=L1;k++) { g +=a[k4]*a[i+k]; 
                                   if (k<j)  z = 1; else z = k;
                                   k4 += z; }
             k3 += j; g *= s; e[j] = g; h += a[i+j]*g;
           }
          h *= 0.5*s; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             s = a[i+j]; g = e[j]-h*s; e[j] = g;
             for (k=1;k<=j;k++) { a[k3] += -s*e[k]-a[i+k]*g; k3++; }
           }
        }
      h = d[L]; d[L] = a[i+L]; a[i+L] = h;
     }
    h = d[1]; d[1] = a[1]; a[1] = h; e[1] = 0.0; e2[1] = 0.0; s = d[n];
    t = fabs(e[n]); mn = s - t; ma = s + t;
    for (i=n-1;i>=1;i--)
     {
      u = fabs(e[i]); h = t + u; t = u; s = d[i]; u = s - h;
      if (u < mn) mn = u;
      u = s + h;
      if (u > ma) ma = u;
     }
    for (i=1;i<=n;i++) { Lb[i] = mn; x[i] = ma; }
    norm = fabs(mn); s = fabs(ma);
    if (s>norm) norm = s;
    w = ma; lambda = norm;
    for (k=k2;k>=k1;k--)
     {
      eps = 7.28e-17*norm; s = mn; i = k;
      do {cond = 0; g = Lb[i];
         if (s < g) s = g; else { i--; if (i>=k1) cond = 1; }
      } while (cond);
      g = x[k];
      if (w>g) w = g;
      while (w-s>2.91e-16*(fabs(s)+fabs(w))+eps)
       {
         L1 = 0; g = 1.0; t = 0.5*(s+w);
         for (i=1;i<=n;i++)
          {
            if (g!=0)  g = e2[i] / g; else g = fabs(6.87e15*e[i]);
            g = d[i]-t-g;
            if (g<0) L1++;
          }
         if (L1<k1) { s = t; Lb[k1] = s; }
          else
           { if (L1<k)
               {
                 s = t; Lb[L1+1] = s;
                 if (x[L1]>t) x[L1] = t;
               }
              else w = t;
           }
      }
      u = 0.5*(s+w); x[k] = u;
      /* if  (!(( ceil(u) - u  < 10e-5 ) || ( u - floor(u) < 10e-5 ))) { return 0; }; */
    }
  } 
  return 1;
}


/* obliczenia eigenvalues end */


/* konwersja macierzy kwadratowej do postaci koniecznej do obliczenia spektrum 
   dodatkowo jest zerowana glowna przekatna - gdyz w programie mogla byc nadpisana przez 
   procedure sprawdzajaca spojnosc grafu. 

*/
void AToa(int N, ta A, long double * a) {
int poz,i,j;

for (i = 0 ; i < N ; i++) A[i][i] = 0;
a[0] = 0.0;
poz = 1; 
for ( i = 0 ; i < N; i++)
 {
  for ( j = 0; j <= i ; j++)
   { if (A[i][j]) a[poz++]=1.0; else a[poz++]=0.0;}
 }

}

int isintegral(int N, long double * x) {
int i;
 long double u; 
 for (i = 1; i < N; i++) {
  u = x[i];
  if  (!(( ceil(u) - u  < 10e-5 ) || ( u - floor(u) < 10e-5 ))) { return 0; };
 }
 return 1; 
}


int integral(int N, ta A)
{
 long double a[400];
 long double x[21];
 int st; 
 int sum;
 int poz, i, j;

 for (i = 0 ; i < N ; i++) A[i][i] = 0;
 
 a[0] = 0.0;
 poz = 1; 
 for ( i = 0 ; i < N; i++)
 {
  for ( j = 0; j <= i ; j++)
   { if (A[i][j]==1) a[poz++]=1.0; else a[poz++]=0.0;}
 }
 return eigensymmatrix(N, a, 1, N,  x);
}

void spoj2(int N, ta A, int k)
{
  int i;
  for ( i = 0 ; i < N; i++)
   if ((i != k)  && (A[i][k] == 1) && (A[i][i] == 0)) { A[i][i]=1; spoj2( N, A , i); }
} /* spoj2 */

int connected(int N, ta A)
{
 int i, Result;
 A[0][0] = 1;
 for (i = 1; i < N; i++) A[i][i] = 0;
 
 Result = 1; 
 for (i = 0; i < N ; i++)
  {
   if (A[i][i] == 0) { Result = 0; break; }
   spoj2(N,A,i);
  }
  return Result;  
} /* connected */

/** format .g6, for n <= 6 */ 
void BMKdecode(char * BUFFOR, int *N, ta A)
{
  int bit, poz, i, j;
  bit = 32;
  poz = 1;
 
  *N   = BUFFOR[0] - 63;
  for (i = 1; i < *N; i++)
   for (j = 0; j< i; j++)
    {
      if (bit == 0) { bit = 32;  poz++; }
      if ((BUFFOR[poz] - 63) & bit)
                { A[i][j] = A[j][i] = 1; }
               else
                { A[i][j] = A[j][i] = 0; }
      bit = bit >> 1;
    }
}

void dop(int N, ta A)
{ 
 int i,j;
 for (i = 0; i < N - 1; i++)
  for ( j = i + 1; j < N; j++)
   { A[i][j] = A[j][i] = 1 - A[j][i]; }
}

/** 

  Uwaga! Wielokrotne otwieranie pliku do zapisu wprowadza opoznienie
*/
void printSout(int N, ta A, char * filename)
{
 FILE *out;
 int i, j; 
 out = fopen(filename,"ab");
 for (i = 0; i < N - 1; i++)
  { 
    for ( j = i + 1; j < N; j++)
    if (A[i][j]) putc('1',out); else putc('0',out);
  } 
  putc('\n',out);
 fclose(out);
}

/* spektra sa w odwrotnym porzadku np. -1, -1, -1, ..., n-1*/
void printX(int N, long double * x) {
int i; 
  for (i = 1; i <= N; i++ ) printf("%.3Lf, ",x[i]);
printf("\n");
}
