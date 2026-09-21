#include <chrono>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "graph6.hpp"
#include "spectrum.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

static volatile sig_atomic_t stop_requested = 0;
static void handle_sigint(int) { stop_requested = 1; }

static void report(double t, long long processed, long long survivors, const char* tag) {
  long long rate = t > 0 ? (long long)(processed / t) : 0;
  std::cerr << "[sito] " << tag << " t=" << t << "s graphs=" << processed
            << " survivors=" << survivors << " rate=" << rate << "/s\n";
}

int main(int argc, char* argv[]) {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  struct sigaction sa;
  sa.sa_handler = handle_sigint;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = 0;
  sigaction(SIGINT, &sa, nullptr);

  auto t0 = std::chrono::steady_clock::now();
  auto last = t0;
  long long processed = 0, survivors = 0;

#ifdef _OPENMP
  int t = omp_get_max_threads();
  if (argc > 1) t = atoi(argv[1]);
  if (t < 1) t = 1;
  omp_set_num_threads(t);

  const std::ptrdiff_t batch_size = 64 * t;
  std::vector<std::string> batch;
  batch.reserve(batch_size);
  std::string line;

  while (true) {
    batch.clear();
    while ((int)batch.size() < batch_size && std::getline(std::cin, line))
      batch.push_back(std::move(line));
    if (batch.empty()) break;

    std::vector<std::string> found;
    #pragma omp parallel
    {
      std::vector<uint64_t> adj;
      std::vector<std::string> local;
      #pragma omp for schedule(static) nowait
      for (int g = 0; g < (int)batch.size(); g++)
        if (decode_graph6(batch[g], adj) && eigensymmatrix(batch[g].data()))
          local.push_back(std::move(batch[g]));
      #pragma omp critical
      found.insert(found.end(), local.begin(), local.end());
    }
    for (const auto& s : found) std::cout << s << '\n';
    processed += (long long)batch.size();
    survivors += (long long)found.size();
    if (stop_requested) break;

    auto now = std::chrono::steady_clock::now();
    if (now - last >= std::chrono::seconds(10)) {
      last = now;
      report(std::chrono::duration<double>(now - t0).count(), processed, survivors, "progress");
    }
    if ((int)batch.size() < batch_size) break;
  }
#else
  std::string line;
  std::vector<uint64_t> adj;
  while (std::getline(std::cin, line)) {
    processed++;
    if (decode_graph6(line, adj) && eigensymmatrix(line.data())) {
      survivors++;
      std::cout << line << '\n';
    }
    if (stop_requested) break;
  }
#endif

  std::cout.flush();
  report(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(),
         processed, survivors, "done");
  return 0;
}
