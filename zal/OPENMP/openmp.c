#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <omp.h>

#define NMAX 32
#define DEFAULT_BATCH_SIZE 65536
#define BUFSIZE 1024
#define DEFAULT_CLASSES 256
#define MAX_CLASSES 1024
#define MAX_HITS 4096
#define LINE_LEN 96

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
    t = fabsl(e[n]); mn = s - t; ma = s + t;
    for (i=n-1;i>=1;i--)
     {
      u = fabsl(e[i]); h = t + u; t = u; s = d[i]; u = s - h;
      if (u < mn) mn = u;
      u = s + h;
      if (u > ma) ma = u;
     }
    for (i=1;i<=n;i++) { Lb[i] = mn; x[i] = ma; }
    norm = fabsl(mn); s = fabsl(ma);
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
      while (w-s>2.91e-16*(fabsl(s)+fabsl(w))+eps)
       {
         if (floor(w+10e-5)<s-10e-5) return 0;  // przedział nie zawiera liczby całkowitej
         L1 = 0; g = 1.0; t = 0.5*(s+w);
         for (i=1;i<=n;i++)
          {
            if (g!=0)  g = e2[i] / g; else g = fabsl(6.87e15*e[i]);
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

// config file: one finished class number per line
void read_done(const char *config_file, int *done) {
    FILE *f = fopen(config_file, "r");
    if (!f) return;
    int r;
    while (fscanf(f, "%d", &r) == 1)
        if (r >= 0 && r < MAX_CLASSES)
            done[r] = 1;
    fclose(f);
}

void mark_done(const char *config_file, int r) {
    FILE *f = fopen(config_file, "a");
    if (f) {
        fprintf(f, "%d\n", r);
        fclose(f);
    }
}

int main(int argc, char **argv) {
    int n = 0, k = 0, all_graphs = 0, batch_size = DEFAULT_BATCH_SIZE;
    int debug = 0, mod = DEFAULT_CLASSES;
    int i = 1;

    if (argc < 3) {
        printf("%s: n k [-a] [-d] [-b batch] [-m classes] [-B threads]\n", argv[0]);
        return 1;
    }

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-a") == 0) all_graphs = 1;
        else if (strcmp(argv[i], "-d") == 0) debug = 1;
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) batch_size = atoi(argv[++i]);
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) mod = atoi(argv[++i]);
        else if (strcmp(argv[i], "-B") == 0 && i + 1 < argc) i++;  // same CLI as cuda, not used here
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

    if (mod < 1 || mod > MAX_CLASSES) {
        printf("bad class count %d (max %d)\n", mod, MAX_CLASSES);
        return 1;
    }

    if (batch_size < 1) {
        printf("batch size must be >= 1\n");
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

    char result_file[256];
    snprintf(result_file, sizeof(result_file), "result_%d_%d.txt", n, k);

    FILE *rf = fopen(result_file, "a");
    if (!rf) {
        printf("could not open %s\n", result_file);
        return 1;
    }

    int done[MAX_CLASSES] = {0};
    read_done(config_file, done);

    char **bufory = malloc(batch_size * sizeof(char *));
    if (!bufory) {
        printf("malloc fail\n");
        fclose(rf);
        return 1;
    }
    for (i = 0; i < batch_size; i++) {
        bufory[i] = malloc(BUFSIZE);
        if (!bufory[i]) {
            printf("malloc fail\n");
            free(bufory);
            fclose(rf);
            return 1;
        }
    }

    static char hits[MAX_HITS][LINE_LEN];

    long total_found = 0;
    int classes_done = 0;

    for (i = 0; i < mod; i++)
        if (done[i]) classes_done++;

    printf("searching n=%d, k=%d, %s, classes %d\n", n, k, all_graphs ? "all" : "connected", mod);
    if (classes_done > 0)
        printf("resuming, %d/%d classes done\n", classes_done, mod);

    for (int r = 0; r < mod && !stop_flag; r++) {
        if (done[r]) continue;

        char cmd[512];
        if (all_graphs)
            snprintf(cmd, sizeof(cmd), "geng %d %d:%d %d/%d -g 2>/dev/null", n, k, k, r, mod);
        else
            snprintf(cmd, sizeof(cmd), "geng -c %d %d:%d %d/%d -g 2>/dev/null", n, k, k, r, mod);

        FILE *fp = popen(cmd, "r");
        if (!fp) {
            printf("could not run geng\n");
            fclose(rf);
            return 1;
        }

        if (debug)
            printf("class %d:\n", r);

        long processed_count = 0;
        long class_found = 0;
        int nhits = 0;
        int graphs;

        do {
            for (graphs = 0; graphs < batch_size; graphs++) {
                if (stop_flag) break;
                if (fgets(bufory[graphs], BUFSIZE, fp) == NULL) break;
            }

            if (graphs == 0) break;

            long batch_found = 0;

            #pragma omp parallel for default(none) reduction(+:batch_found) schedule(dynamic) shared(bufory, graphs, hits, nhits, rf)
            for (int j = 0; j < graphs; j++) {
                if (eigensymmatrix(bufory[j])) {
                    batch_found++;
                    #pragma omp critical
                    {
                        if (nhits < MAX_HITS)
                            strcpy(hits[nhits++], bufory[j]);
                        else
                            fputs(bufory[j], rf);
                    }
                }
            }

            processed_count += graphs;
            class_found += batch_found;

            if (debug) {
                printf("\rprocessed: %ld, found: %ld", processed_count, class_found);
                fflush(stdout);
            }
        } while (graphs == batch_size);

        pclose(fp);

        if (stop_flag) {
            printf("\ninterrupted in class %d, discarding it\n", r);
            break;  // class not finished: hits discarded, not marked, redone on next run
        }

        if (debug)
            printf("\n");

        // class fully processed: save its hits and mark it done
        for (int j = 0; j < nhits; j++)
            fputs(hits[j], rf);
        fflush(rf);

        sigset_t block_mask, old_mask;
        sigemptyset(&block_mask);
        sigaddset(&block_mask, SIGINT);
        sigaddset(&block_mask, SIGTERM);
        sigprocmask(SIG_BLOCK, &block_mask, &old_mask);
        mark_done(config_file, r);
        sigprocmask(SIG_SETMASK, &old_mask, NULL);

        classes_done++;
        total_found += class_found;
    }

    printf("done, classes: %d/%d, found: %ld\n", classes_done, mod, total_found);

    for (i = 0; i < batch_size; i++)
        free(bufory[i]);
    free(bufory);
    fclose(rf);

    return 0;
}
