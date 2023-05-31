#include <spdlog/spdlog.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <istream>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "codewords_trie.hpp"
#include "lemmas.hpp"
#include "suffix_tree.hpp"

void naive_pattern_matching(const std::string& snippet_before,
                            const std::vector<std::string>& snippets,
                            const std::string& snippet_after) {}

struct CodewordInfo {
  CodeType k;
  size_t counter;
  size_t start_pos;
};

struct NaiveAlgo {
  NaiveAlgo(const CodewordsTrie& trie) : _trie(trie) {}

  std::vector<size_t> step2(const std::vector<CodewordInfo>& codewords) {
    std::vector<size_t> matches;
    if (codewords.empty()) {
      return matches;
    }

    std::cout << "start_pos: " << codewords[0].start_pos << "\n";

    const CodewordsTrie::Node& snippet_node = _trie._nodes[codewords[0].k];
    std::cout << "s_0: " << snippet_node.full_string << "\n";

    std::cout << snippet_node << "\n";

    size_t l = snippet_node.len_of_suffix_intersection_with_prefix_of_pattern;
    // size_t l =
    // snippet_node.len_of_prefix_intersection_with_suffix_of_pattern;

    size_t k = 1;

    // TODO: here add only from middle.
    size_t len_of_leftover_snippets = 0;
    for (auto elem : codewords) {
      len_of_leftover_snippets += _trie._nodes.at(elem.k).str_codeword_size;
    }

    size_t curr_offset = 0;

    const std::string& pattern = _trie._useful_structs.pattern;
    while (k < codewords.size() &&
           l + len_of_leftover_snippets >= pattern.size()) {
      std::cout << "\n\n\n";
      std::cout << "-----> l: " << l << ", k: " << k
                << ", curr_offset: " << curr_offset << "\n";
      std::string p1l = pattern.substr(0, l);
      std::cout << "p1l: `" << p1l << "`\n";
      const CodewordsTrie::Node& prev_codeword_node =
          _trie._nodes[codewords[k].k - 1];
      const CodewordsTrie::Node& codeword_node = _trie._nodes[codewords[k].k];
      std::cout << "s_" << k << ": `" << codeword_node.full_string << "`\n";
      std::cout << codeword_node << "\n";

      if (!codeword_node.is_snippet) {
        size_t start_pos = codewords[0].start_pos;
        for (auto elem : codeword_node.pattern_occurences) {
          matches.push_back(elem + start_pos + curr_offset);
        }
      }

      if (l == 0) {
        std::cout << "l is zero, skipping...\n";
        l = codeword_node.len_of_suffix_intersection_with_prefix_of_pattern;
        ++k;
        curr_offset += prev_codeword_node.str_codeword_size;
        continue;
      }
      std::cout << "calling lemma28\n";
      std::vector<size_t> res = lemma28_simplified(
          pattern.substr(0, l), codeword_node, _trie._useful_structs);

      if (!res.empty()) {
        size_t start_pos = codewords[0].start_pos;
        for (auto elem : res) {
          matches.push_back(elem + start_pos + curr_offset);
        }
        l = codeword_node.len_of_suffix_intersection_with_prefix_of_pattern;
        ++k;
        curr_offset += prev_codeword_node.str_codeword_size;
        continue;
      }

      if (k + 1 == codewords.size()) {
        break;
      }

      assert(l < _trie._useful_structs.pattern.size());

      std::cout << "checking if p[1...l]s_k is a prefix of p...\n";

      std::string substr_after_p1l = pattern.substr(l, std::string::npos);
      std::cout << "substr_after_p1l: `" << substr_after_p1l << "`\n";
      std::string s = _trie._useful_structs.st.get_longest_common_prefix(
          codeword_node.full_string, substr_after_p1l);

      if (s.size() == codeword_node.full_string.size()) {
        std::cout << "p[1...l]s_k is a prefix of p!\n";
        std::cout << "p[1...l]: `" << p1l << "`, s_k: `"
                  << codeword_node.full_string << "`, pattern: `" << pattern
                  << "`\n";
        l += codeword_node.full_string.size();
        ++k;
        curr_offset += prev_codeword_node.str_codeword_size;
        continue;
      }

      std::cout << "calling lemma29\n";
      size_t pos_in_concat =
          lemma29(pattern.substr(0, l), codeword_node, _trie._useful_structs);

      std::cout << "calculated pos_in_concat: " << pos_in_concat << "\n";
      // if (pos_in_concat == -1) {
      //   l = _trie._useful_structs.half_snippets_info.at(l);
      //   continue;
      // }

      assert(pos_in_concat <= l + codeword_node.str_codeword_size);
      if (pos_in_concat == l + codeword_node.str_codeword_size) {
        l = 0;
      } else {
        l = l + codeword_node.str_codeword_size - pos_in_concat;
      }
      std::cout << "set l: " << l << "\n";
      ++k;
      curr_offset += prev_codeword_node.str_codeword_size;
    }

    std::cout << "\nread " << k << " codewords\n\n=============\n\n";
    return matches;
  }

  // void step(CodeType codeword) {
  //   std::cout << "naive algo step\n";
  //   std::cout << codeword << "\n";

  //   if (!initialized) {
  //     longest_prefix_of_p_ending_s_1 =
  //         _trie._nodes.at(codeword)
  //             .len_of_suffix_intersection_with_prefix_of_pattern;
  //     initialized = true;
  //     return;
  //   }

  //   bool res = lemma28_simplified();

  //   if (prefix_of_p) {
  //     l += codeword.size;
  //     ++k;
  //     return;
  //   }

  //   b = lemma29();

  //   if (b == -1) {
  //     l = 0;
  //     return;
  //   }

  //   l = b + codeword.size;
  //   ++k;
  // }

  void step3(const std::vector<CodewordInfo>& elems) {
    for (auto elem : elems) {
      std::string s = _trie._nodes.at(elem.k).full_string;
      std::cout << "<" << elem.counter << ", " << s << ", " << elem.start_pos
                << "> ";
    }
    std::cout << "\n";
    return;
  }

  // bool initialized = false;
  // int k = 2;
  size_t longest_prefix_of_p_ending_s_1;
  const CodewordsTrie& _trie;
  std::vector<CodeType> codewords;
};

void search(std::istream& is, const UsefulStructs& useful_structs) {
  CodewordsTrie trie(useful_structs);
  NaiveAlgo algo(trie);
  CodeType k;

  bool is_in_sequence_state = false;

  std::string decompressed_string;

  std::vector<CodewordInfo> shit_sequence;

  size_t start_pos = 0;

  size_t counter = 0;

  CodewordInfo cw_info;
  std::optional<CodewordInfo> prev_cw_info;

  std::vector<size_t> matches;

  while (is.read(reinterpret_cast<char*>(&k), sizeof(CodeType))) {
    // std::cout << "================================================\n";
    // std::cout << "read codeword: " << k << "\n";
    trie.step(k);
    const CodewordsTrie::Node& node = trie._nodes.at(k);
    // std::cout << "read codeword str: " << trie._nodes.at(k).full_string
    //           << ", counter: " << counter << "\n";
    cw_info = {k, counter, start_pos};
    ++counter;
    /*
      here we should consider these cases for is_snippet value:
      ...010...
      10...
      ...01010...
    */
    if (node.is_snippet) {
      // std::cout << "here1\n";
      is_in_sequence_state = true;
      if (prev_cw_info.has_value()) {
        // std::cout << "here1.1\n";
        // algo.step(prev_codeword.value());
        shit_sequence.push_back(prev_cw_info.value());
        prev_cw_info.reset();
      }
      // algo.step(k);
      shit_sequence.push_back(cw_info);
    } else if (is_in_sequence_state) {
      // std::cout << "here2\n";
      // algo.step(k);
      shit_sequence.push_back(cw_info);
      is_in_sequence_state = false;
      algo.step3(shit_sequence);

      std::vector<size_t> res = algo.step2(shit_sequence);
      for (auto elem : res) {
        matches.push_back(elem);
      }

      shit_sequence.clear();
      prev_cw_info = cw_info;
    } else {
      // std::cout << "here3\n";
      prev_cw_info = cw_info;
    }

    decompressed_string += node.full_string;
    // if (node.is_snippet) {
    //   snippets_sequence.push_back(k);
    // }
    start_pos += node.str_codeword_size;
  }
  // algo.step2(snippets_sequence);

  std::cout << "reporting matches:\n";
  for (auto elem : matches) {
    std::cout << "match: " << elem << "\n";
  }

  std::cout << "\nreference search:\n";
  size_t pos = decompressed_string.find(useful_structs.pattern, 0);
  while (pos != std::string::npos) {
    std::cout << "found occurence of the pattern: " << pos << "\n";
    pos = decompressed_string.find(useful_structs.pattern, pos + 1);
  }
}

int main(int argc, char* argv[]) {
  spdlog::info("Welcome to spdlog!");
  std::string pattern = "abcd";
  // std::string pattern = "aabasfaabasfasf";
  // std::string pattern = "abababa";
  assert(pattern.find("$") == std::string::npos);
  std::string prefix_snippet = "abacab";
  std::string snippet = "cabad";

  std::string reversed_pattern(pattern.rbegin(), pattern.rend());
  auto useful_structs = UsefulStructs(pattern, reversed_pattern);
  std::ifstream input_file(argv[1], std::ios_base::binary);
  search(input_file, useful_structs);

  // int res = lemma29("a", "a", useful_structs);
  // std::cout << "res: " << res << "\n";
  return 0;
}
