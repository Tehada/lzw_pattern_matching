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

void naive_pattern_matching(const std::string& snippet_before,
                            const std::vector<std::string>& snippets,
                            const std::string& snippet_after) {}

struct NaiveAlgo {
  void step(const CodewordsTrie::Node& node) {
    std::cout << "naive algo step\n";
  }
};

void search(std::istream& is, const UsefulStructs& useful_structs) {
  NaiveAlgo algo;
  CodewordsTrie trie(useful_structs);
  CodeType k;

  while (is.read(reinterpret_cast<char*>(&k), sizeof(CodeType))) {
    std::cout << "================================================\n";
    std::cout << "read codeword: " << k << "\n";
    trie.step(k);
    const CodewordsTrie::Node& node = trie._nodes.at(k);
    if (node.is_snippet) {
      algo.step(node);
    }
  }
}

int main(int argc, char* argv[]) {
  // std::string pattern = "abcd";
  std::string pattern = "aab";
  assert(pattern.find("$") == std::string::npos);
  std::string prefix_snippet = "abacab";
  std::string snippet = "cabad";

  std::string reversed_pattern(pattern.rbegin(), pattern.rend());
  auto useful_structs = UsefulStructs(pattern, reversed_pattern);
  std::ifstream input_file(argv[1], std::ios_base::binary);
  search(input_file, useful_structs);
  return 0;
}
