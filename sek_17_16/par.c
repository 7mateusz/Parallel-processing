#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define N_VERTICES 17
#define K_EDGES 36
#define BATCH_SIZE 4096

int is_integral_graph(char *g6_buffer) {
    int i, j, k, L, L1, z;
    int k3, k4;
    long double g, h, ma, mn, norm, s, t, u, w, eps;
    
    long double d[N_VERTICES + 1], e[N_VERTICES + 1], e2[N_VERTICES + 1];
    long double lb[N_VERTICES + 1], ub[N_VERTICES + 1];
    long double adj_matrix[N_VERTICES * (N_VERTICES - 1) / 2 + N_VERTICES + 1];
    
    int n = g6_buffer[0] - 63; 
    int bit = 32;
    int g6_ptr = 1;
    int mat_ptr = 1;

    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            if (i == j) {
                adj_matrix[mat_ptr++] = 0.0;
            } else {
                if (bit == 0) {
                    bit = 32;
                    g6_ptr++;
                }
                adj_matrix[mat_ptr++] = ((g6_buffer[g6_ptr] - 63) & bit) ? 1.0 : 0.0;
                bit >>= 1;
            }
        }
    }

    i = 0;
    for (L = 1; L <= n; L++) {
        i += L;
        d[L] = adj_matrix[i];
    }

    for (L = n; L >= 2; L--) {
        i--; j = i; h = adj_matrix[j]; s = 0;
        for (k = L - 2; k >= 1; k--) {
            i--; g = adj_matrix[i]; s += g * g;
        }
        i--;
        if (s == 0) {
            e[L] = h; e2[L] = h * h; adj_matrix[j] = 0.0;
        } else {
            s += h * h; e2[L] = s; g = sqrtl(s);
            if (h >= 0.0) g = -g;
            e[L] = g;
            s = 1.0 / (s - h * g);
            adj_matrix[j] = h - g; h = 0.0; L1 = L - 1; k3 = 1;
            for (j = 1; j <= L1; j++) {
                k4 = k3; g = 0;
                for (k = 1; k <= L1; k++) {
                    g += adj_matrix[k4] * adj_matrix[i + k];
                    z = (k < j) ? 1 : k;
                    k4 += z;
                }
                k3 += j; g *= s; e[j] = g; h += adj_matrix[i + j] * g;
            }
            h *= 0.5 * s; k3 = 1;
            for (j = 1; j <= L1; j++) {
                s = adj_matrix[i + j]; g = e[j] - h * s; e[j] = g;
                for (k = 1; k <= j; k++) {
                    adj_matrix[k3] += -s * e[k] - adj_matrix[i + k] * g;
                    k3++;
                }
            }
        }
        h = d[L]; d[L] = adj_matrix[i + L]; adj_matrix[i + L] = h;
    }
    
    h = d[1]; d[1] = adj_matrix[1]; adj_matrix[1] = h;
    e[1] = 0.0; e2[1] = 0.0;
    
    t = fabsl(e[n]); mn = d[n] - t; ma = d[n] + t;
    for (i = n - 1; i >= 1; i--) {
        u = fabsl(e[i]); h = t + u; t = u; s = d[i];
        if (s - h < mn) mn = s - h;
        if (s + h > ma) ma = s + h;
    }
    
    for (i = 1; i <= n; i++) {
        lb[i] = mn; ub[i] = ma;
    }
    
    norm = fabsl(mn) > fabsl(ma) ? fabsl(mn) : fabsl(ma);
    eps = 7.28e-17L * norm;
    w = ma;

    for (k = n; k >= 1; k--) {
        s = mn; i = k;
        int cond;
        do {
            cond = 0;
            if (s < lb[i]) s = lb[i];
            else { i--; if (i >= 1) cond = 1; }
        } while (cond);
        
        if (w > ub[k]) w = ub[k];

        while (w - s > 2.91e-16L * (fabsl(s) + fabsl(w)) + eps) {
            if (floorl(w + 1e-4L) < s - 1e-4L) return 0;

            L1 = 0; g = 1.0; t = 0.5L * (s + w);
            for (i = 1; i <= n; i++) {
                g = (g != 0) ? (e2[i] / g) : fabsl(6.87e15L * e[i]);
                g = d[i] - t - g;
                if (g < 0) L1++;
            }

            if (L1 < 1) { 
                s = t; lb[1] = s; 
            } else if (L1 < k) {
                s = t; lb[L1 + 1] = s;
                if (ub[L1] > t) ub[L1] = t;
            } else {
                w = t;
            }
        }
        
        u = 0.5L * (s + w);
        ub[k] = u;

        if (fabsl(u - roundl(u)) > 1e-5L) return 0;
    }

    return 1;
}

int main() {
    char cmd[256];
    sprintf(cmd, "geng -c %d %d:%d -g", N_VERTICES, K_EDGES, K_EDGES);
    
    FILE *fp = popen(cmd, "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not execute geng pipe.\n");
        return 1;
    }

    char **batch = malloc(BATCH_SIZE * sizeof(char *));
    char *line = NULL;
    size_t line_len = 0;
    long total_integral_found = 0;
    int current_batch_count = 0;

    printf("searching integral graphs (n=%d, k=%d)...\n", N_VERTICES, K_EDGES);

    while (getline(&line, &line_len, fp) != -1) {
        line[strcspn(line, "\r\n")] = 0;
        batch[current_batch_count++] = strdup(line);

        if (current_batch_count == BATCH_SIZE) {
            #pragma omp parallel for reduction(+:total_integral_found) schedule(dynamic)
            for (int i = 0; i < current_batch_count; i++) {
                if (is_integral_graph(batch[i])) {
                    total_integral_found++;
                    #pragma omp critical
                    printf("Found integral graph: %s\n", batch[i]);
                }
                free(batch[i]);
            }
            current_batch_count = 0;
        }
    }

    if (current_batch_count > 0) {
        #pragma omp parallel for reduction(+:total_integral_found) schedule(dynamic)
        for (int i = 0; i < current_batch_count; i++) {
            if (is_integral_graph(batch[i])) {
                total_integral_found++;
                #pragma omp critical
                printf("Found integral graph: %s\n", batch[i]);
            }
            free(batch[i]);
        }
    }

    printf("\nSearch finished.\n");
    printf("Integral graphs found: %ld\n", total_integral_found);

    free(batch);
    free(line);
    pclose(fp);
    return 0;
}
