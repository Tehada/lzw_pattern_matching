#pragma once

#include "suffix_tree.hpp"

using CodeType = std::uint16_t;

struct CodewordsTrie {
  struct Node {
    Node() {}

    Node(char first_letter) : first_letter(first_letter) {}
    bool is_snippet = false;
    SuffTreeNodeInfo node_in_suffix_tree;
    // bool is_pattern_inside;
    // size_t len_of_prefix_intersection;
    // size_t len_of_suffix_intersection;

    // CodeType k;
    char symbol_which_brought_here;

    // Временно буду хранить тут всю строку, чтобы легче было дебажить. В
    // продовой реализации это поле надо будет убрать:
    std::string full_string;

    // // храним здесь только первую букву, чтобы память не росла:
    char first_letter;
    std::unordered_map<char, size_t> childs;
  };

  CodewordsTrie(SuffixTree* st) : _st(st) {
    const long int minc = std::numeric_limits<char>::min();
    const long int maxc = std::numeric_limits<char>::max();
    for (long int c = minc; c <= maxc; ++c) {
      _nodes.emplace_back(c);
      CodeType new_codeword = _nodes.size() - 1;
      std::cout << "added codeword: " << new_codeword << "\n";
      _nodes.back().full_string += static_cast<char>(c);
      std::cout << "finding next_node for char: `" << static_cast<char>(c)
                << "`\n";
      SuffTreeNodeInfo next_node = _st->get_next_node(_st->get_root_node(), c);
      std::cout << next_node << "\n";
      if (next_node.node_id != -1) {
        _nodes.back().is_snippet = true;
        _nodes.back().node_in_suffix_tree = next_node;
        std::cout << "set is_snippet for codeword: <" << c << ", `"
                  << _nodes.back().full_string << "`>\n";
      }
    }
  }

  void step(CodeType k) {
    if (_curr_node_id == -1) {
      _curr_node_id = k;
      return;
    }
    char first_letter = _nodes.at(k).first_letter;
    std::cout << "here: <first_letter: `" << first_letter
              << "`, full_string: " << _nodes[_curr_node_id].full_string
              << ", curr_node_id: " << _curr_node_id << ", k: " << k << ">\n";
    if (!_nodes[_curr_node_id].childs.contains(first_letter)) {
      _nodes.emplace_back(_nodes[_curr_node_id].first_letter);
      CodeType new_codeword = _nodes.size() - 1;
      _nodes.back().full_string =
          _nodes[_curr_node_id].full_string + first_letter;
      std::cout << "added codeword: " << new_codeword << " for `"
                << _nodes.back().full_string << "`, first_letter: `"
                << _nodes.back().first_letter << "`\n";
      _nodes[_curr_node_id].childs.emplace(first_letter, _nodes.size() - 1);
      if (_nodes[_curr_node_id].is_snippet) {
        std::cout << "here2, checking is_snippet for: `"
                  << _nodes.back().full_string << "`\n";
        SuffTreeNodeInfo next_node = _st->get_next_node(
            _nodes[_curr_node_id].node_in_suffix_tree, first_letter);
        if (next_node.node_id != -1) {
          _nodes.back().is_snippet = true;
          _nodes.back().node_in_suffix_tree = next_node;
          std::cout << "set is_snippet for codeword: <" << new_codeword << ", `"
                    << _nodes.back().full_string << "`>\n";
        }
      }
    }
    _curr_node_id = k;
  }

  int _curr_node_id = -1;
  std::vector<Node> _nodes;
  SuffixTree* _st;
};

std::ostream& operator<<(std::ostream& os, const CodewordsTrie::Node& node) {
  os << "CodewordsTrie::Node: <is_snippet: " << node.is_snippet
     << ", node_in_suffix_tree: " << node.node_in_suffix_tree
     << ", full_string: " << node.full_string
     << ", first_letter: " << node.first_letter << ">";
  return os;
}
