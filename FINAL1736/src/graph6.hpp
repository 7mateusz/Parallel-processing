#pragma once

#include <string>
#include <vector>
#include <cstdint>

inline bool decode_graph6(const std::string& line, std::vector<uint64_t>& adj) {
  // getting number of the vertices
  unsigned int n = (unsigned char)line[0] - 63;
  if (n > 62) return false;

  // prepping the matrix
  adj.assign(n, 0);

  int i = 0, j = 1, word_i = 0, dec_word, dec_i = 6;

  while (i < n && j < n) {
    // iterating through next bits of grph6 representation (-63)
    if (dec_i == 6) {
      dec_i = 0;
      word_i++;
      dec_word = (unsigned char)line[word_i] - 63;
    }
    uint64_t bit = (dec_word >> (5 - dec_i)) & 1;
    adj[i] |= bit << (n - j - 1); // top tri
    adj[j] |= bit << (n - i - 1); // bot tri

    i++;
    // go to next column if we step on the diagonal
    if (i == j) {
      j++; i = 0;
    }

    dec_i++;
  }
  
  return true;
}
