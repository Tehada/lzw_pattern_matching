#include <cassert>
#include <fstream>
#include <iostream>
#include <istream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "codewords_trie.hpp"
#include "suffix_tree.hpp"

void fill_lps(const std::string& pattern, std::vector<size_t>* lps) {
  (*lps)[0] = 0;
  size_t len = 0;

  for (size_t i = 0; i + 1 < pattern.size(); ++i) {
    while (pattern[i + 1] != pattern[len] && len > 0) {
      len = (*lps)[len - 1];
    }

    if (pattern[i + 1] == pattern[len]) {
      ++len;
    }
    (*lps)[i + 1] = len;
  }
}

void naive_pattern_matching(const std::string& snippet_before,
                            const std::vector<std::string>& snippets,
                            const std::string& snippet_after) {}

struct NaiveAlgo {
  void step(const CodewordsTrie::Node& node) {
    std::cout << "   -------> naive algo step: " << node << "\n";
  }
};

void search(std::istream& is, SuffixTree* st) {
  NaiveAlgo algo;
  CodewordsTrie trie(st);
  CodeType k;

  while (is.read(reinterpret_cast<char*>(&k), sizeof(CodeType))) {
    std::cout << "read codeword: " << k << "\n";
    trie.step(k);
    const CodewordsTrie::Node& node = trie._nodes.at(k);
    if (node.is_snippet) {
      algo.step(node);
    }
  }
}

int main(int argc, char* argv[]) {
  std::string pattern = "abacabadabdl";
  assert(pattern.find("$") == std::string::npos);
  std::string prefix_snippet = "abacab";
  std::string snippet = "cabad";

  std::string reversed_pattern(pattern.rbegin(), pattern.rend());

  std::vector<size_t> lps(pattern.size());
  fill_lps(pattern, &lps);

  std::vector<size_t> reversed_lps(reversed_pattern.size());
  fill_lps(reversed_pattern, &reversed_lps);

  auto st = SuffixTree(pattern, true);
  auto reversed_st = SuffixTree(reversed_pattern);

  std::ifstream input_file(argv[1], std::ios_base::binary);
  search(input_file, &st);
  return 0;
}
