#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <signal.h>

#include "graph6.hpp"
#include "spectrum.hpp"

static volatile sig_atomic_t stop_requested = 0;
static void handle_sigint(int) { stop_requested = 1; }

static void report(double t, long long processed, long long survivors, const char* tag) {
  long long rate = t > 0 ? (long long)(processed / t) : 0;
  std::cerr << "[sito] " << tag << " t=" << t << "s graphs=" << processed
            << " survivors=" << survivors << " rate=" << rate << "/s\n";
}

int main() {
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

  std::string line;
  std::vector<uint64_t> adj;
  while (std::getline(std::cin, line)) {
    processed++;
    if (decode_graph6(line, adj) && eigensymmatrix(line.data())) {
      survivors++;
      std::cout << line << '\n';
    }
    if (stop_requested) break;
    if ((processed & 0xFFF) == 0) {
      auto now = std::chrono::steady_clock::now();
      if (now - last >= std::chrono::seconds(10)) {
        last = now;
        report(std::chrono::duration<double>(now - t0).count(), processed, survivors, "progress");
      }
    }
  }

  std::cout.flush();
  report(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(),
         processed, survivors, "done");
  return 0;
}
