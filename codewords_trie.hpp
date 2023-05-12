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

struct UsefulStructs {
  std::string pattern, reversed_pattern;
  std::vector<size_t> lps, reversed_lps;
  SuffixTree st, reversed_st;

  UsefulStructs(const std::string& pattern, const std::string& reversed_pattern)
      : pattern(pattern), st(pattern, true), reversed_st(reversed_pattern) {
    lps.resize(pattern.size());
    fill_lps(pattern, &(lps));

    std::cout << "lps:\n";
    for (auto elem : lps) {
      std::cout << elem << " ";
    }
    std::cout << "\n====\n";

    reversed_lps.resize(pattern.size());
    fill_lps(reversed_pattern, &(reversed_lps));
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
    // bool is_pattern_inside;
    size_t p_i = 0;

    size_t str_codeword_size = 0;

    CodeType codeword;
    char symbol_which_brought_here;

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
      std::cout << "added codeword: " << new_codeword << "\n";
      _nodes.back().full_string += static_cast<char>(c);
      std::cout << "finding next_node for char: `" << static_cast<char>(c)
                << "`\n";
      SuffTreeNodeInfo next_node = _useful_structs.st.get_next_node(
          _useful_structs.st.get_root_node(), c);
      _nodes.back().node_in_suffix_tree = next_node;
      std::cout << next_node << "\n";
      if (next_node.node_id != -1) {
        if (next_node.edge == nullptr &&
            _useful_structs.st.is_marked_node(next_node.node_id)) {
          _nodes.back().len_of_prefix_intersection_with_suffix_of_pattern =
              _nodes.back().str_codeword_size;
        }
        _nodes.back().is_snippet = true;
        std::cout << "set is_snippet for codeword: <" << c << ", `"
                  << _nodes.back().full_string << "`>\n";
      }
      if (_useful_structs.pattern.at(0) == static_cast<char>(c)) {
        _nodes.back().len_of_suffix_intersection_with_prefix_of_pattern = 1;
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
     << node.len_of_suffix_intersection_with_prefix_of_pattern << ">";
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
    _nodes.back().full_string = _nodes[_curr_node_id].full_string + new_letter;
    _nodes[_curr_node_id].children.emplace(new_letter, new_codeword);
    SuffTreeNodeInfo next_node = _useful_structs.st.get_next_node(
        _nodes[_curr_node_id].node_in_suffix_tree, new_letter);
    _nodes.back().node_in_suffix_tree = next_node;
    if (next_node.node_id != -1) {
      if (next_node.edge == nullptr &&
          _useful_structs.st.is_marked_node(next_node.node_id)) {
        _nodes.back().len_of_prefix_intersection_with_suffix_of_pattern =
            _nodes.back().str_codeword_size;
      }

      if (_nodes[_curr_node_id].is_snippet) {
        _nodes.back().is_snippet = true;
      }
    } else {
      _nodes.back().len_of_prefix_intersection_with_suffix_of_pattern =
          _nodes[_curr_node_id]
              .len_of_prefix_intersection_with_suffix_of_pattern;
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

    std::cout << "added codeword: " << _nodes.back() << "\n\n";
  }

  std::cout << "\n";
  std::cout << "curr: " << _nodes.at(_curr_node_id) << "\n";
  std::cout << "read: " << _nodes.at(k) << "\n";
  std::cout << "\n";

  _curr_node_id = k;
}