# Instruction for CUDA agent

## Project: Integral graph generator (CUDA version)

### What this is

This program searches for integral graphs (graphs whose adjacency matrix has all
integer eigenvalues). It uses `geng` (from nauty) to generate graphs in graph6
format, then checks each graph on the GPU using a CUDA kernel implementing the
tridiagonalization + Sturm sequence eigenvalue algorithm.

### Dependencies

```
geng (nauty)      - graph generator, must be in PATH
CUDA toolkit      - nvcc compiler, CUDA runtime
NVIDIA GPU        - any CUDA-capable GPU (sm_60+)
```

Install nauty on Ubuntu/Debian: `sudo apt install nauty`

### Files

```
CUDA/cuda.cu       - main program (kernel + host code)
CUDA/Makefile      - build rules
CUDA/INSTRUCTION.md - this file
OPENMP/openmp.c    - OpenMP reference version (no GPU needed, NMAX=32)
```

### Compile

```bash
cd CUDA
make
```

Or manually:
```bash
nvcc -O3 -arch=sm_60 -o cuda cuda.cu
```

If your GPU architecture is different, change `-arch=sm_60` accordingly:
- GTX 10xx / Quadro P: sm_60 or sm_61
- RTX 20xx: sm_75
- RTX 30xx: sm_80 or sm_86
- RTX 40xx: sm_89

### Run

```bash
./cuda n k [-a] [-b batch_size]
```

- `n` - number of vertices (1..20, limited by NMAX=20)
- `k` - number of edges
- `-a` - all graphs (default: connected only via `geng -c`)
- `-b N` - batch size, how many graphs to process per GPU launch (default: 4096)

Examples:
```bash
./cuda 10 15                     # connected graphs, n=10 k=15
./cuda 10 15 -a                  # all graphs (not just connected)
./cuda 12 20 -b 8192             # custom batch size
```

### Output files

- `result_n_k.txt` - found integral graphs in graph6 format (appended)
- `config_n_k.txt` - resume checkpoint (auto-created on Ctrl+C, auto-deleted on finish)

### How it works

1. Parse args (n, k, flags) — same manual parsing as OPENMP/openmp.c
2. Read `config_n_k.txt` to get skip count from previous interrupted run
3. Open `geng` pipe to generate graphs
4. Read graphs in batches of `batch_size` into a flat host buffer
5. `cudaMemcpy` the whole batch to GPU (one transfer)
6. Launch `eigensymmatrix_kernel` with N threads (one per graph)
   - Each thread checks one graph for integrality
   - Writes 1 or 0 to `d_wyniki` array on GPU
7. `cudaDeviceSynchronize()` and copy results back to host
8. Host iterates results, writes integral graphs to `result_n_k.txt`
9. On Ctrl+C: saves `processed_count` to config, exits cleanly
10. On resume: skips `processed_count` graphs, continues

### Key differences from OpenMP version

| Feature | OPENMP/openmp.c | CUDA/cuda.cu |
|---|---|---|
| NMAX | 32 | 20 |
| Math precision | long double (80-bit) | double (64-bit) |
| Batch transfer | none (in-place check) | cudaMemcpy whole batch |
| Checking | CPU threads via OpenMP | GPU kernel |

### Notes

- NMAX=20 because the kernel uses static local arrays; larger NMAX would
  increase register pressure and local memory usage, degrading GPU performance.
- For `n > 20`, use the OpenMP version with `-m dynamic` (not implemented in
  current openmp.c — needs NMAX bump or dynamic mode).
- The kernel uses `double` precision throughout (not float). The math functions
  are `sqrt`, `fabs`, `floor`, `ceil` (not the float `*f` variants).
- Every CUDA API call is wrapped in `HANDLE_ERROR` macro for academic error
  checking.
- Ctrl+C during skip phase is blocked (signals masked) to avoid inconsistent
  config state.
- Ctrl+C during kernel execution will finish the current batch before stopping
  (GPU kernels are non-interruptible).

### Verification

Compare results with the OpenMP version for small n,k:
```bash
cd OPENMP && make && ./openmp 10 15 && cat result_10_15.txt
cd CUDA && make && ./cuda 10 15 && cat result_10_15.txt
```

Results should match. Known reference: 4 integral connected graphs for n=10, k=15.
