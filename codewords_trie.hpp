#pragma once

#include "suffix_tree.hpp"

using CodeType = std::uint16_t;

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

std::vector<size_t> lemma41(const std::string& pattern,
                            const std::vector<size_t>& lps);

struct UsefulStructs {
  std::string pattern, reversed_pattern;
  std::vector<size_t> lps, reversed_lps;
  SuffixTree st, reversed_st;
  std::vector<size_t> half_snippets_info;

  UsefulStructs(const std::string& pattern, const std::string& reversed_pattern)
      : pattern(pattern), st(pattern, false), reversed_st(reversed_pattern) {
    lps.resize(pattern.size());
    fill_lps(pattern, &(lps));

    // std::cout << "lps:\n";
    // for (auto elem : lps) {
    //   std::cout << elem << " ";
    // }
    // std::cout << "\n====\n";

    reversed_lps.resize(pattern.size());
    fill_lps(reversed_pattern, &(reversed_lps));

    half_snippets_info = lemma41(pattern, lps);
  }
};

struct CodewordsTrie {
  struct Node {
    Node() {}

    Node(char first_letter) : first_letter(first_letter) {}
    bool is_snippet = false;
    SuffTreeNodeInfo node_in_suffix_tree;
    size_t len_of_prefix_intersection_with_suffix_of_pattern = 0;
    size_t len_of_suffix_intersection_with_prefix_of_pattern = 0;
    std::vector<size_t> pattern_occurences;

    bool is_prefix = false;
    bool is_suffix = false;

    size_t str_codeword_size = 0;

    CodeType codeword;
    char letter;

    // Временно буду хранить тут всю строку, чтобы легче было дебажить. В
    // продовой реализации это поле надо будет убрать:
    std::string full_string;

    // // храним здесь только первую букву, чтобы память не росла:
    char first_letter;
    std::unordered_map<char, size_t> children;
  };

  CodewordsTrie(const UsefulStructs& useful_structs)
      : _useful_structs(useful_structs) {
    const long int minc = std::numeric_limits<char>::min();
    const long int maxc = std::numeric_limits<char>::max();
    for (long int c = minc; c <= maxc; ++c) {
      _nodes.emplace_back(c);
      _nodes.back().str_codeword_size = 1;
      CodeType new_codeword = _nodes.size() - 1;
      _nodes.back().codeword = new_codeword;
      // std::cout << "added codeword: " << new_codeword << "\n";
      _nodes.back().full_string += static_cast<char>(c);
      _nodes.back().letter += static_cast<char>(c);
      // std::cout << "finding next_node for char: `" << static_cast<char>(c)
      //           << "`\n";
      SuffTreeNodeInfo next_node = _useful_structs.st.get_next_node(
          _useful_structs.st.get_root_node(), c);
      _nodes.back().node_in_suffix_tree = next_node;
      // std::cout << next_node << "\n";
      if (next_node.node_id != -1) {
        if (next_node.edge == nullptr &&
            _useful_structs.st.is_marked_node(next_node.node_id)) {
          _nodes.back().len_of_prefix_intersection_with_suffix_of_pattern =
              _nodes.back().str_codeword_size;
          _nodes.back().is_suffix = true;
        }
        _nodes.back().is_snippet = true;
        // std::cout << "set is_snippet for codeword: <" << c << ", `"
        //           << _nodes.back().full_string << "`>\n";
      }
      if (_useful_structs.pattern.at(0) == static_cast<char>(c)) {
        _nodes.back().len_of_suffix_intersection_with_prefix_of_pattern = 1;
        _nodes.back().is_prefix = true;
      };
    }
  }

  void step(CodeType k);

  int _curr_node_id = -1;
  std::vector<Node> _nodes;
  const UsefulStructs& _useful_structs;
};

std::ostream& operator<<(std::ostream& os, const CodewordsTrie::Node& node) {
  os << "CodewordsTrie::Node: <codeword: " << node.codeword
     << ", is_snippet: " << node.is_snippet
     << ", node_in_suffix_tree: " << node.node_in_suffix_tree
     << ", full_string: " << node.full_string
     << ", first_letter: " << node.first_letter << ", prefix len: "
     << node.len_of_prefix_intersection_with_suffix_of_pattern
     << ", suffix len: "
     << node.len_of_suffix_intersection_with_prefix_of_pattern
     << ", is_prefix: " << node.is_prefix << ", is_suffix: " << node.is_suffix
     << ", pattern_occurences: <";
  for (auto elem : node.pattern_occurences) {
    os << elem << ", ";
  }
  os << ">>";
  return os;
}

void CodewordsTrie::step(CodeType k) {
  if (_curr_node_id == -1) {
    _curr_node_id = k;
    return;
  }

  char new_letter;
  if (k == _nodes.size()) {
    new_letter = _nodes.at(_curr_node_id).first_letter;
  } else {
    new_letter = _nodes.at(k).first_letter;
  }

  if (!_nodes[_curr_node_id].children.contains(new_letter)) {
    _nodes.emplace_back(_nodes[_curr_node_id].first_letter);
    _nodes.back().str_codeword_size =
        _nodes[_curr_node_id].str_codeword_size + 1;
    CodeType new_codeword = _nodes.size() - 1;
    _nodes.back().codeword = new_codeword;
    _nodes.back().letter = new_letter;
    _nodes.back().full_string = _nodes[_curr_node_id].full_string + new_letter;
    _nodes[_curr_node_id].children.emplace(new_letter, new_codeword);
    // The idea of constructing of the field
    // len_of_prefix_intersection_with_suffix_of_pattern is to take the current
    // codeword and to make a step in the suffix tree with a new letter. If we
    // have made a successful step, it means that the new codeword has this
    // prefix-suffix intersection on 1 letter longer than it was in the current
    // codeword. We want to consider only proper suffixes of the pattern (the
    // length of the suffix should be smaller than the length of the pattern).
    // Hence, if we already have this prefix-suffix intersection equal to
    // `pattern.size() - 1` we don't want to make steps in the suffix tree any
    // more.

    // // remove
    // assert(_nodes[_curr_node_id]
    //                .len_of_prefix_intersection_with_suffix_of_pattern +
    //            1 <=
    //        _useful_structs.pattern.size());

    // By default we will initialize the size of the intersection with the value
    // of the current codeword. It will have at least this value. It is possible
    // that the real value for intersection will be +1 from the current and we
    // are checking below for this case.

    // // remove
    // _nodes.back().len_of_prefix_intersection_with_suffix_of_pattern =
    //     _nodes[_curr_node_id].len_of_prefix_intersection_with_suffix_of_pattern;

    // Do we want to know the suff tree node for every codeword? Once we have a
    // -1 node for some codeword, all its children will also have -1 node in
    // suff tree. Right?
    SuffTreeNodeInfo next_node = _useful_structs.st.get_next_node(
        _nodes[_curr_node_id].node_in_suffix_tree, new_letter);
    _nodes.back().node_in_suffix_tree = next_node;
    // Below we check, whether the prefix-suffix intersection of the current
    // codeword is of the same length as the pattern. If so, it means that we
    // have found a full occurence of the pattern in the beginning of this
    // codeword, but we are interested only in proper suffix of the pattern, so
    // we will not increase the length of the instersection here.
    if (next_node.node_id != -1 && next_node.edge == nullptr &&
        _useful_structs.st.is_marked_node(next_node.node_id) &&
        _nodes.back().str_codeword_size < _useful_structs.pattern.size()) {
      // If we are here, it means that new codeword corresponds to some
      // proper suffix of the pattern.
      _nodes.back().len_of_prefix_intersection_with_suffix_of_pattern =
          _nodes.back().str_codeword_size;
      _nodes.back().is_suffix = true;
      if (_nodes[_curr_node_id].is_snippet) {
        _nodes.back().is_snippet = true;
      }

    } else {
      _nodes.back().len_of_prefix_intersection_with_suffix_of_pattern =
          _nodes[_curr_node_id]
              .len_of_prefix_intersection_with_suffix_of_pattern;
    }

    if (_nodes[_curr_node_id].is_prefix) {
      if (_nodes.back().str_codeword_size < _useful_structs.pattern.size()) {
        assert(_nodes.back().str_codeword_size > 0);
        size_t i = _nodes.back().str_codeword_size - 1;
        if (_useful_structs.pattern.at(i) == new_letter) {
          _nodes.back().is_prefix = true;
        }
      }
    }

    size_t p_i =
        _nodes[_curr_node_id].len_of_suffix_intersection_with_prefix_of_pattern;

    if (p_i == _useful_structs.pattern.size()) {
      p_i -= 1;
      p_i = _useful_structs.lps.at(p_i);
    }

    while (_useful_structs.pattern.at(p_i) != new_letter && p_i > 0) {
      p_i = _useful_structs.lps[p_i - 1];
    }
    if (_useful_structs.pattern.at(p_i) == new_letter) {
      _nodes.back().len_of_suffix_intersection_with_prefix_of_pattern = p_i + 1;
    }

    assert(_nodes.back().len_of_suffix_intersection_with_prefix_of_pattern <=
           _useful_structs.pattern.size());

    // Check for an occurence of the pattern inside this snippet. If we have
    // found an occurence, then remember about it and also adjust suffix-prefix
    // intersection using prefix-function table to spot further occurences
    // correctly.
    if (_nodes.back().len_of_suffix_intersection_with_prefix_of_pattern ==
        _useful_structs.pattern.size()) {
      assert(_nodes.back().str_codeword_size >= _useful_structs.pattern.size());
      _nodes.back().pattern_occurences.push_back(
          _nodes.back().str_codeword_size - _useful_structs.pattern.size());
      assert(_useful_structs.pattern.size() > 0);
      _nodes.back().len_of_suffix_intersection_with_prefix_of_pattern =
          _useful_structs.lps.at(_useful_structs.pattern.size() - 1);
    }

    // if (_nodes.back().full_string == "abcd") {
    //   std::cout << "added codeword: " << _nodes.back() << "\n\n";

    //   std::cout << "\n";
    //   std::cout << "curr: " << _nodes.at(_curr_node_id) << "\n";
    //   std::cout << "read: " << _nodes.at(k) << "\n";
    //   std::cout << "\n";
    // }

    _curr_node_id = k;
  }
}