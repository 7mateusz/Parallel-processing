#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <omp.h>

#define NMAX 32
#define DEFAULT_BATCH_SIZE 4096
#define BUFSIZE 1024

volatile sig_atomic_t stop_flag = 0;

void handle_sigint(int sig) {
    stop_flag = 1;
}

int eigensymmatrix(char * BUFOR)
{
  int i,j,k,k3,k4,L,L1,z;
  long double /* lambda, */ eps,g,h,ma,mn,norm,s,t,u,w;
  int cond;
  long double d[NMAX+1], e[NMAX+1], e2[NMAX+1], Lb[NMAX+1];
  long double x[NMAX+1];
  long double a[NMAX*(NMAX-1)/2 + NMAX + 1];
  int n;
  int bit, poz, poz2;
  bit = 32;
  poz = 1;
  poz2 = 1;
  n = BUFOR[0] - 63;
  a[0] = 0.0;
  for (i = 0; i < n; i++)
   for (j = 0; j<=i; j++)
    {
	if (i==j) {a[poz2++] = 0; }
	else {
         if (bit == 0) { bit = 32;  poz++; }
         if ((BUFOR[poz] - 63) & bit)
                { a[poz2++] = 1; }
               else
                { a[poz2++] = 0; }
         bit = bit >> 1;
	}
    }
 
  int k1 = 1;
  int k2 = n;
  // if ((1<=k1) && (k1<=k2) && (k2<=n))
   {
    i = 0;
    for (L=1;L<=n;L++) { i += L; d[L] = a[i]; /* printf("%Lf ",a[i]); */ }
 
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
    w = ma; /* lambda = norm; */ eps = 7.28e-17*norm;
    for (k=k2;k>=k1;k--)
     {
      /* eps = 7.28e-17*norm; */ s = mn; i = k;
      do {cond = 0; g = Lb[i];
         if (s < g) s = g; else { i--; if (i>=k1) cond = 1; }
      } while (cond);
      g = x[k];
      if (w>g) w = g;
      while (w-s>2.91e-16*(fabs(s)+fabs(w))+eps)
       {
         if (floor(w+10e-5)<s-10e-5) return 0;  // przedział nie zawiera liczby całkowitej
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
      } // while
      u = 0.5*(s+w); x[k] = u;
	  if  (!(( ceil(u) - u  < 10e-5 ) || ( u - floor(u) < 10e-5 ))) { return 0; };
    }
  } 
  return 1;
}
long read_config(const char *config_file) {
    FILE *f = fopen(config_file, "r");
    if (!f) return 0;
    long count = 0;
    fscanf(f, "%ld", &count);
    fclose(f);
    return count;
}

void write_config(const char *config_file, long count) {
    sigset_t mask, old_mask;
    sigemptyset(&mask);
    sigaddset(&mask, SIGINT);
    sigaddset(&mask, SIGTERM);
    sigprocmask(SIG_BLOCK, &mask, &old_mask);

    FILE *f = fopen(config_file, "w");
    if (f) {
        fprintf(f, "%ld\n", count);
        fclose(f);
    }

    sigprocmask(SIG_SETMASK, &old_mask, NULL);
}

int main(int argc, char **argv) {
    int n = 0, k = 0, all_graphs = 0, batch_size = DEFAULT_BATCH_SIZE, debug = 0;
    int i = 1;

    if (argc < 3) {
        printf("%s: n k [-a] [-b batch] [-d]\n", argv[0]);
        return 1;
    }

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-a") == 0) all_graphs = 1;
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) batch_size = atoi(argv[++i]);
        else if (strcmp(argv[i], "-d") == 0) debug = 1;
        else if (n == 0) n = atoi(argv[i]);
        else if (k == 0) k = atoi(argv[i]);
    }

    if (n < 1 || k < 0 || k > n * (n - 1) / 2) {
        printf("bad n=%d or k=%d (max edges: %d)\n", n, k, n * (n - 1) / 2);
        return 1;
    }

    if (n > NMAX) {
        printf("n=%d > NMAX=%d\n", n, NMAX);
        return 1;
    }

    struct sigaction sa;
    sa.sa_handler = handle_sigint;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    char config_file[256];
    snprintf(config_file, sizeof(config_file), "config_%d_%d.txt", n, k);
    long skip_count = read_config(config_file);

    char result_file[256];
    snprintf(result_file, sizeof(result_file), "result_%d_%d.txt", n, k);

    char cmd[512];
    if (all_graphs)
        snprintf(cmd, sizeof(cmd), "geng %d %d:%d -g 2>/dev/null", n, k, k);
    else
        snprintf(cmd, sizeof(cmd), "geng -c %d %d:%d -g 2>/dev/null", n, k, k);

    FILE *fp = popen(cmd, "r");
    if (!fp) {
        printf("could not run geng\n");
        return 1;
    }

    FILE *rf = fopen(result_file, "a");
    if (!rf) {
        printf("could not open %s\n", result_file);
        pclose(fp);
        return 1;
    }

    if (debug) printf("searching n=%d, k=%d, %s\n", n, k, all_graphs ? "all" : "connected");

    if (batch_size < 1) {
        printf("batch size must be >= 1\n");
        fclose(rf); pclose(fp); return 1;
    }

    if (skip_count > 0)
        if (debug) printf("resuming, skipping %ld...\n", skip_count);

    char **bufory = malloc(batch_size * sizeof(char *));
    if (!bufory) { 
		printf("malloc fail\n"); 
		fclose(rf);
		pclose(fp);
		return 1;
	}
    for (int i = 0; i < batch_size; i++) {
        bufory[i] = malloc(BUFSIZE);
        if (!bufory[i]) {
            printf("malloc fail\n");
            free(bufory); 
			fclose(rf); 
			pclose(fp); 
			return 1;
        }
    }

    long processed_count = 0;
    long total_found = 0;
    int graphs;

    {
        sigset_t skip_mask, skip_old;
        sigemptyset(&skip_mask);
        sigaddset(&skip_mask, SIGINT);
        sigaddset(&skip_mask, SIGTERM);
        sigprocmask(SIG_BLOCK, &skip_mask, &skip_old);

        while (processed_count < skip_count) {
            if (fgets(bufory[0], BUFSIZE, fp) == NULL) break;
            processed_count++;
        }

        sigprocmask(SIG_SETMASK, &skip_old, NULL);
    }

    do {
        for (graphs = 0; graphs < batch_size; graphs++) {
            if (stop_flag) break;
            if (fgets(bufory[graphs], BUFSIZE, fp) == NULL) {
                if (stop_flag) break;
                break;
            }
        }

        if (graphs == 0) break;

        long batch_found = 0;

        #pragma omp parallel for default(none) shared(graphs, bufory, rf) reduction(+:batch_found) schedule(dynamic) private(i)
        for (int i = 0; i < graphs; i++) {
            if (eigensymmatrix(bufory[i])) {
                batch_found++;
                #pragma omp critical
                {
                    fputs(bufory[i], rf);
                    fflush(rf);
                }
            }
        }

        processed_count += graphs;
        total_found += batch_found;

        if (debug) {
            printf("\rprocessed: %ld, found: %ld", processed_count, total_found);
            fflush(stdout);
        }

        if (stop_flag) {
            printf("\ninterrupted, saving...\n");
            write_config(config_file, processed_count);
            break;
        }
    } while (graphs == batch_size);

    if (!stop_flag) {
        if (debug) printf("\ndone.\n");
        remove(config_file);
    }

    printf("found: %ld\n", total_found);

    for (int i = 0; i < batch_size; i++)
        free(bufory[i]);
    free(bufory);
    fclose(rf);
    pclose(fp);

    return 0;
}
