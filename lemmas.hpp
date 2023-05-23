#pragma once

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "codewords_trie.hpp"
#include "suffix_tree.hpp"

std::vector<size_t> lemma28_simplified(const std::string& prefix_snippet,
                                       const CodewordsTrie::Node& snippet_node,
                                       const UsefulStructs& useful_structs) {
  const std::string& pattern = useful_structs.pattern;
  const std::string& snippet = snippet_node.full_string;
  std::cout << "pattern: " << pattern << "\n";
  std::cout << "prefix_snippet: " << prefix_snippet << "\n";
  std::cout << "snippet: " << snippet_node << "\n";
  std::string concat = prefix_snippet + snippet;
  std::cout << "concat: " << concat << "\n";

  assert(!prefix_snippet.empty());
  assert(!snippet.empty());

  if (pattern.size() > concat.size()) {
    return {};
  }

  // this part is not described in lemma 2.8 by gawry, but it seems necessary.
  // Basically lemma 2.8 is only concerned about long borders of both snippets.
  // It seems that there is no word on what we should do, if both snippets don't
  // have long border. But it is possible to get a match, if there are no long
  // borders in both snippets still. Particularly, match can happen either from
  // the beginning of the concatenation, or from the end -- only these two cases
  // seem to be possible here. Both of these cases could be checked in O(1) so
  // we will perform the check.

  // first check match from the beginning. From pattern we subtract
  // prefix_snippet and compare the resulting part with snippet.
  //
  // | prefix_snippet|| snippet|
  // |        pattern       |
  //                  |xxxxx|
  //                     ^
  //                     |
  // the part that we want to check.
  assert(pattern.size() > prefix_snippet.size());
  std::string substr_to_check = pattern.substr(prefix_snippet.size());
  std::cout << "[simple check] comparing `" << substr_to_check
            << "` with left edge of snippet `" << snippet << "`\n";
  std::string lcpp =
      useful_structs.st.get_longest_common_prefix(substr_to_check, snippet);
  if (lcpp.size() == substr_to_check.size()) {
    std::cout << "We have found match from the left edge!\n";
    return {0};
  }

  // TODO: here we can actually find a last match, when checking right edge, but
  // this could be not the only match and, hence, we will return not the first
  // occurence of pattern.

  // perform a similar check from the end of concat.
  //
  // | prefix_snippet|| snippet|
  //           |        pattern       |
  //           |xxxxx|
  //              ^
  //              |
  // the part that we want to check.
  //
  // to check for equality of 2 substrings here we will use query on reversed
  // suffix tree.
  assert(pattern.size() > snippet.size());
  substr_to_check = pattern.substr(0, pattern.size() - snippet.size());
  std::cout << "[simple check] comparing `" << substr_to_check
            << "` with right edge of prefix_snippet `" << prefix_snippet
            << "`\n";
  std::string reversed_prefix_snippet(prefix_snippet.rbegin(),
                                      prefix_snippet.rend());
  std::string reversed_substr_to_check(substr_to_check.rbegin(),
                                       substr_to_check.rend());
  lcpp = useful_structs.reversed_st.get_longest_common_prefix(
      reversed_substr_to_check, reversed_prefix_snippet);
  if (lcpp.size() == substr_to_check.size()) {
    std::cout << "We have found match from the right edge!\n";
    return {concat.size() - pattern.size()};
  }
}

int lemma28(const std::string& pattern, const std::string& prefix_snippet,
            const std::string& suffix_snippet, const std::vector<size_t>& lps,
            const std::vector<size_t>& reversed_lps, const SuffixTree& st,
            const SuffixTree& reversed_st) {
  std::cout << "pattern: " << pattern << "\n";
  std::cout << "prefix_snippet: " << prefix_snippet << "\n";
  std::cout << "suffix_snippet: " << suffix_snippet << "\n";
  std::string concat = prefix_snippet + suffix_snippet;
  std::cout << "concat: " << concat << "\n";

  assert(!prefix_snippet.empty());
  assert(!suffix_snippet.empty());
  assert(prefix_snippet.size() >= suffix_snippet.size());

  if (pattern.size() > concat.size()) {
    return -1;
  }

  // this part is not described in lemma 2.8 by gawry, but it seems necessary.
  // Basically lemma 2.8 is only concerned about long borders of both snippets.
  // It seems that there is no word on what we should do, if both snippets don't
  // have long border. But it is possible to get a match, if there are no long
  // borders in both snippets still. Particularly, match can happen either from
  // the beginning of the concatenation, or from the end -- only these two cases
  // seem to be possible here. Both of these cases could be checked in O(1) so
  // we will perform the check.

  // first check match from the beginning. From pattern we subtract
  // prefix_snippet and compare the resulting part with suffix_snippet.
  //
  // | prefix_snippet|| suffix_snippet|
  // |        pattern       |
  //                  |xxxxx|
  //                     ^
  //                     |
  // the part that we want to check.
  assert(pattern.size() > prefix_snippet.size());
  std::string substr_to_check = pattern.substr(prefix_snippet.size());
  std::cout << "[simple check] comparing `" << substr_to_check
            << "` with left edge of suffix_snippet `" << suffix_snippet
            << "`\n";
  std::string lcpp =
      st.get_longest_common_prefix(substr_to_check, suffix_snippet);
  if (lcpp.size() == substr_to_check.size()) {
    std::cout << "We have found match from the left edge!\n";
    return 0;
  }

  // TODO: here we can actually find a last match, when checking right edge, but
  // this could be not the only match and, hence, we will return not the first
  // occurence of pattern.

  // perform a similar check from the end of concat.
  //
  // | prefix_snippet|| suffix_snippet|
  //           |        pattern       |
  //           |xxxxx|
  //              ^
  //              |
  // the part that we want to check.
  //
  // to check for equality of 2 substrings here we will use query on reversed
  // suffix tree.
  assert(pattern.size() > suffix_snippet.size());
  substr_to_check = pattern.substr(0, pattern.size() - suffix_snippet.size());
  std::cout << "[simple check] comparing `" << substr_to_check
            << "` with right edge of prefix_snippet `" << prefix_snippet
            << "`\n";
  std::string reversed_prefix_snippet(prefix_snippet.rbegin(),
                                      prefix_snippet.rend());
  std::string reversed_substr_to_check(substr_to_check.rbegin(),
                                       substr_to_check.rend());
  lcpp = reversed_st.get_longest_common_prefix(reversed_substr_to_check,
                                               reversed_prefix_snippet);
  if (lcpp.size() == substr_to_check.size()) {
    std::cout << "We have found match from the right edge!\n";
    return concat.size() - pattern.size();
  }

  return -1;

  std::cout << "lps:\n";
  for (auto elem : lps) {
    std::cout << elem << " ";
  }
  std::cout << "\n";

  std::cout << "reversed_lps:\n";
  for (auto elem : reversed_lps) {
    std::cout << elem << " ";
  }
  std::cout << "\n";

  size_t x = lps.at(prefix_snippet.size() - 1);
  size_t y = reversed_lps.at(suffix_snippet.size() - 1);

  std::cout << "borders, x: " << x << ", y: " << y << "\n";

  if (x + y < pattern.size()) {
    std::cout << "borders x and y are too short: x + y < m, leaving...\n";
    return -1;
  }

  size_t d = prefix_snippet.size() - lps.at(prefix_snippet.size() - 1);
  std::cout << "period d for prefix_snippet: " << d << ", " << prefix_snippet
            << "\n";

  std::cout << "calculating p[1...k] by finding longest common prefix of: `"
            << pattern.substr(d) << "` and `" << pattern << "`\n";

  std::string substr_of_pattern = pattern.substr(d);
  std::string str = st.get_longest_common_prefix(substr_of_pattern, pattern);

  size_t k = pattern.size() - substr_of_pattern.size() + str.size();
  std::cout << "k: " << k << ", pattern: " << pattern
            << ", p1k: " << pattern.substr(0, k) << "\n";
  std::string p1k = pattern.substr(0, k);
  std::cout << "longest common prefix for `" << substr_of_pattern << "` and `"
            << pattern << "` is: `" << str << "`\n";

  std::cout << "p[1...k]: " << p1k << "\n";

  int shift = std::min(prefix_snippet.size(), concat.size() - pattern.size());

  std::cout << "shift: " << shift << ", d: " << d
            << ", shift / d: " << shift / d
            << ", final_shift / d * d: " << shift / d * d << "\n";

  int final_shift = shift / d * d;

  std::cout << "looking for leftmost mismatch between p1k: `" << p1k
            << "` and p[j...m]: `" << suffix_snippet << "` in concat: `"
            << concat << "`\n";

  std::cout << concat << "\n";
  std::cout << std::string(final_shift, '_') + p1k << "\n";

  // find leftmost mismatch between shifted p1k and pjm
  std::optional<bool> leftmost_mismatch_was_found;
  if (p1k.size() + final_shift <= prefix_snippet.size()) {
    // means that we don't even have any intersection between shifted p[1...k]
    // and p[j...m], hence no mismatch.
    leftmost_mismatch_was_found = false;
  } else {
    assert(prefix_snippet.size() >= final_shift);
    size_t r = prefix_snippet.size() - final_shift + 1;
    assert(r >= 1);
    std::cout << "r: " << r << "\n";
    std::string prm = pattern.substr(r - 1);
    std::cout << "prm: `" << prm << "`\n";
    std::cout << "pjm: `" << suffix_snippet << "`\n";
    std::string str1 = st.get_longest_common_prefix(prm, suffix_snippet);
    std::cout << "longest common prefix for `" << prm << "` and `"
              << suffix_snippet << "` is `" << str1 << "`\n";
    if (str1.size() == prm.size()) {
      leftmost_mismatch_was_found = false;
    } else {
      throw std::runtime_error("Not implemented!");
    }
  }

  if (leftmost_mismatch_was_found.has_value()) {
    if (!leftmost_mismatch_was_found.value()) {
      if (p1k.size() == pattern.size()) {
        // TODO: here we can find not the first occurence, adjust shift on
        // period to return the first occurence.
        std::cout << "We have found occurence of pattern in concat!\n";
        std::cout << "Exact occurence is for shift " << final_shift
                  << " inside concat `" << concat << "`\n";
        std::cout << concat << "\n";
        std::cout << std::string(final_shift, '_') + pattern << "\n";
        return -1;

      } else {
        std::cout << "pattern `" << pattern << "` doesn't occur inside concat `"
                  << concat << "`!\n";
      }
    } else {
      throw std::runtime_error("Not implemented!");
    }
  } else {
    throw std::runtime_error("Not implemented!");
  }

  return -1;
}

size_t lemma29(const std::string& prefix_snippet,
               const CodewordsTrie::Node& snippet_node,
               const UsefulStructs& useful_structs) {
  const std::string& pattern = useful_structs.pattern;
  const std::string& snippet = snippet_node.full_string;

  std::cout << "pattern: " << pattern << "\n";
  std::cout << "prefix_snippet: " << prefix_snippet << "\n";
  std::cout << "snippet: " << snippet << "\n";
  std::string concat = prefix_snippet + snippet;
  std::cout << "concat: " << concat << "\n";

  size_t d =
      prefix_snippet.size() - useful_structs.lps.at(prefix_snippet.size() - 1);
  std::cout << "period d for prefix_snippet: " << d << ", " << prefix_snippet
            << "\n";

  std::cout << "calculating p[1...ii] by finding longest common prefix of: `"
            << pattern.substr(d) << "` and `" << pattern << "`\n";

  std::string substr_of_pattern = pattern.substr(d);
  std::string str =
      useful_structs.st.get_longest_common_prefix(substr_of_pattern, pattern);

  size_t ii = pattern.size() - substr_of_pattern.size() + str.size();
  std::cout << "ii: " << ii << ", pattern: " << pattern
            << ", p1ii: " << pattern.substr(0, ii) << "\n";
  std::string p1ii = pattern.substr(0, ii);
  std::cout << "longest common prefix for `" << substr_of_pattern << "` and `"
            << pattern << "` is: `" << str << "` which_is_a_prefix_of_pn";

  std::cout << "p[1...ii]: " << p1ii << "\n";

  assert(prefix_snippet.size() + snippet.size() >= pattern.size());
  int shift = prefix_snippet.size() + snippet.size() - pattern.size();
  shift = (shift + d) / d * d;
  std::cout << "smallest possible shift: " << shift << "\n";
  std::cout << prefix_snippet << " -- prefix_snippet\n";
  std::cout << concat << " -- concat\n";
  std::cout << std::string(shift, '_') + p1ii << " -- shifted p1ii\n";
  std::cout << std::string(shift, '_') + pattern << " -- shifted pattern\n";

  // find leftmost mismatch between shifted p1k and snippet
  std::optional<bool> leftmost_mismatch_was_found;
  if (p1ii.size() + shift <= prefix_snippet.size()) {
    throw std::runtime_error("Not implemented!");
    // means that we don't even have any intersection between shifted p[1...k]
    // and p[j...m], hence no mismatch.
    leftmost_mismatch_was_found = false;
  } else {
    assert(prefix_snippet.size() >= shift);
    size_t r = prefix_snippet.size() - shift + 1;
    assert(r >= 1);
    std::cout << "r: " << r << "\n";
    std::string prm = pattern.substr(r - 1);
    std::cout << "prm: `" << prm << "`\n";
    std::cout << "pjk: `" << snippet << "`\n";
    std::string str1 =
        useful_structs.st.get_longest_common_prefix(prm, snippet);
    if (str1.size() == snippet.size()) {
      if (ii == pattern.size() || ii + 1 + shift > concat.size()) {
        return shift;
      }
      std::cout << "mismatch was not found. shift is: " << shift << "\n";
      std::cout << concat << "\n";
      std::cout << std::string(shift, '_') + pattern << "\n";
      return shift;
    } else {
      return 0;
    }
  }

  return 0;
}

// TODO: write tests for this method.
std::vector<size_t> lemma41(const std::string& pattern,
                            const std::vector<size_t>& lps) {
  std::vector<size_t> res;
  res.push_back(1);

  size_t curr_longest_prefix = 1;
  for (size_t i = 2; i <= pattern.size(); ++i) {
    size_t left = i / 2 - 1, right = i - 1;
    // std::cout << "left: " << left << ", right: " << right << "\n";
    std::string substr = pattern.substr(left, right - left + 1);
    // std::cout << "substr: " << substr << "\n";
    // First check is for the case, when we have out of bounds on the left. The
    // second check when on the right we got a mismatch.
    size_t substr_size = right - left + 1;

    if (right % 2 == 1 && curr_longest_prefix == substr_size) {
      curr_longest_prefix = lps.at(curr_longest_prefix - 1);
    }

    while (pattern.at(curr_longest_prefix) != pattern.at(right) &&
           curr_longest_prefix > 0) {
      curr_longest_prefix = lps.at(curr_longest_prefix - 1);
    }

    if (pattern.at(curr_longest_prefix) == pattern.at(right)) {
      ++curr_longest_prefix;
    }

    res.push_back(curr_longest_prefix);
    // std::cout << "added: " << pattern.substr(0, curr_longest_prefix) << "\n";
    // std::cout << "curr_longest_prefix: " << curr_longest_prefix << "\n";
    // std::cout << "\n\n";
  }
  return res;
}

int main2(int argc, char* argv[]) {
  // // test1
  // std::string pattern = "abcabfa";
  // std::string prefix_snippet = pattern.substr(0, 5);
  // std::string suffix_snippet = pattern.substr(2);

  // test 1.5
  // std::string pattern = "abacaba";

  // // test2
  // std::string pattern = "abcabfa";
  // std::string prefix_snippet = pattern.substr(0, 5);
  // std::string suffix_snippet = pattern.substr(3);

  // // test3 -- passing, but we are not returning the first occurence.
  // std::string pattern = "ababab";
  // std::string prefix_snippet = pattern.substr(0, 5);
  // std::string suffix_snippet = pattern.substr(1);

  // // test4
  // std::string pattern = "abcabcabcab";
  // std::string prefix_snippet = "abcabcabc";
  // std::string suffix_snippet = "abcab";

  // // test5, lemma29
  // std::string pattern = "abacabadab";
  // std::string prefix_snippet = "abacab";
  // std::string snippet = "acabad";

  // test6, lemma29
  // std::string pattern = "abacabadab";
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

  // int pos = lemma28(pattern, prefix_snippet, suffix_snippet, lps,
  // reversed_lps,
  //                   st, reversed_st);
  // std::cout << "position of occurence: " << pos << "\n";

  // std::string::size_type pos_etalon;
  // std::string concat = prefix_snippet + suffix_snippet;
  // pos_etalon = concat.find(pattern);

  // int pos_etalon_int;
  // if (pos_etalon == std::string::npos) {
  //   pos_etalon_int = -1;
  // } else {
  //   pos_etalon_int = pos_etalon;
  // }

  // std::cout << "etalon position of occurence: " << pos_etalon_int << "\n";

  // int pos = lemma29(pattern, prefix_snippet, snippet, lps, reversed_lps, st,
  //                   reversed_st);

  std::ifstream input_file(argv[1], std::ios_base::binary);
  // search(input_file, &st);
  return 0;
}