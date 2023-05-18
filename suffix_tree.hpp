#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "rmq/pm_rmq.hpp"

class SuffixTree;

struct Node {
  int suffix_node = -1;
};

struct Edge {
  Edge() {}
  Edge(int first_char_index, int last_char_index, int source_node_index,
       int dest_node_index)
      : first_char_index(first_char_index),
        last_char_index(last_char_index),
        source_node_index(source_node_index),
        dest_node_index(dest_node_index) {}

  int first_char_index, last_char_index, source_node_index, dest_node_index;

  size_t size() const {
    // if (leads_to_leaf) {
    //   assert(last_char_index >= first_char_index);
    //   return last_char_index - first_char_index;
    // }

    // return last_char_index - first_char_index + 1;
    return last_char_index - first_char_index;
  }

  bool leads_to_leaf = false;
};

struct Suffix {
  Suffix(int source_node_index, int first_char_index, int last_char_index)
      : source_node_index(source_node_index),
        first_char_index(first_char_index),
        last_char_index(last_char_index) {}

  int source_node_index, first_char_index, last_char_index;

  size_t size() const { return last_char_index - first_char_index; }

  bool is_explicit() const { return first_char_index > last_char_index; }

  bool is_implicit() const { return last_char_index >= first_char_index; }
};

struct Key {
  int first;
  char second;
};

struct KeyHash {
  std::size_t operator()(const Key& k) const {
    return (std::hash<int>()(k.first) << 8) | std::hash<char>()(k.second);
  }
};

struct KeyEqual {
  bool operator()(const Key& lhs, const Key& rhs) const {
    return lhs.first == rhs.first && lhs.second == rhs.second;
  }
};

struct SuffTreeNodeInfo {
  int node_id = -1;
  const Edge* edge = nullptr;
  size_t offset = 0;
  const SuffixTree* sf_ptr = nullptr;
};

std::ostream& operator<<(std::ostream& os, const SuffTreeNodeInfo& info);

struct SuffixTree {
  std::string _string;
  Suffix _active;
  int _N;
  std::unordered_map<Key, Edge, KeyHash, KeyEqual> _edges;
  std::vector<Node> _nodes;
  // necessary for dfs for euler path:
  std::vector<std::unordered_set<int>> _nodes_childs;

  std::vector<int> _euler;
  std::vector<ssize_t> _level;
  std::unique_ptr<pm_rmq<std::vector<ssize_t>::const_iterator>> _rmq;

  std::vector<int> _repr;
  std::vector<std::unordered_set<Edge*>> _edges_from_node;
  std::vector<std::pair<int, int>> _nodes_info;
  bool _print_debug_info;

  std::unordered_set<int> _leaf_nodes;

  std::unordered_set<int> _marked_nodes;

  SuffixTree(const std::string& str, bool print_debug_info = false)
      : _active(0, 0, -1), _nodes(1), _print_debug_info(print_debug_info) {
    _string = str + "$";
    _N = _string.size() - 1;

    std::unordered_set<int> printed_dest_node;
    for (int i = 0; i < _string.size(); i++) {
      _add_prefix(i);
      for (auto& [key, edge] : _edges) {
        if (!printed_dest_node.contains(edge.dest_node_index)) {
          printed_dest_node.insert(edge.dest_node_index);
          if (print_debug_info) {
            std::cout << "first seen: " << edge.source_node_index << ", "
                      << edge.dest_node_index << ", " << edge.first_char_index
                      << ", " << edge.last_char_index << "\n";
          }
        }
      }
    }
    for (auto& [key, edge] : _edges) {
      if (edge.source_node_index >= _nodes_childs.size()) {
        _nodes_childs.resize(edge.source_node_index + 1);
      }
      _nodes_childs[edge.source_node_index].insert(edge.dest_node_index);

      if (edge.source_node_index >= _edges_from_node.size()) {
        _edges_from_node.resize(edge.source_node_index + 1);
      }
      _edges_from_node[edge.source_node_index].insert(&edge);

      if (_leaf_nodes.contains(edge.dest_node_index)) {
        edge.leads_to_leaf = true;
        if (edge.size() == 0) {
          _marked_nodes.insert(edge.source_node_index);
        } else {
          _marked_nodes.insert(edge.dest_node_index);
        }
      }

      if (print_debug_info) {
        std::cout << "---> " << edge.source_node_index << " <--> "
                  << edge.dest_node_index << ", str: "
                  << _string.substr(
                         edge.first_char_index,
                         edge.last_char_index - edge.first_char_index + 1)
                  << ", suffix_node: " << _nodes[key.first].suffix_node << " | "
                  << edge.first_char_index << ", " << edge.last_char_index
                  << "\n";
      }
    }

    if (print_debug_info) {
      for (auto node : _marked_nodes) {
        std::cout << "marked node: " << node << "\n";
      }

      for (auto node : _leaf_nodes) {
        std::cout << "leaf node: " << node << "\n";
      }
    }

    if (print_debug_info) {
      for (auto& edges : _edges_from_node) {
        for (Edge* edge : edges) {
          std::cout << "---> " << edge->source_node_index << ", "
                    << edge->dest_node_index << "\n";
        }
      }
    }

    _build_euler_path(0, 0, std::back_inserter(_euler),
                      std::back_inserter(_level));

    _build_str_info_for_nodes(0, 0, -1);

    // std::cout << "euler:\n";
    // for (auto elem : _euler) {
    //   std::cout << elem << " ";
    // }
    // std::cout << "\n";

    // std::cout << "level:\n";
    // for (auto elem : _level) {
    //   std::cout << elem << " ";
    // }
    // std::cout << "\n";

    _rmq.reset(new pm_rmq<std::vector<ssize_t>::const_iterator>(_level.begin(),
                                                                _level.end()));
  }

  template <typename OutputIterator1, typename OutputIterator2>
  void _build_euler_path(int node, ssize_t level, OutputIterator1 eulerit,
                         OutputIterator2 levelit) {
    *eulerit++ = node;
    *levelit++ = level;
    if (node >= _nodes_childs.size()) {
      return;
    }
    // std::cout << "-- node: " << node << "\n-- childs: ";
    // for (auto elem : _nodes_childs[node]) {
    //   std::cout << elem << " ";
    // }
    // std::cout << "\n";

    if (node >= _repr.size()) {
      _repr.resize(node + 1);
    }
    _repr[node] = _level.end() - 1 - _level.begin();
    std::for_each(_nodes_childs[node].begin(), _nodes_childs[node].end(),
                  [this, node, level, &eulerit, &levelit](int node_child) {
                    _build_euler_path(node_child, level + 1, eulerit, levelit);
                    *eulerit++ = node;
                    *levelit++ = level;
                  });
  }

  void _build_str_info_for_nodes(int node, int substr_len,
                                 int substr_last_char_index) {
    if (_print_debug_info) {
      std::cout << "dfs: " << node << ", " << substr_len << ", "
                << substr_last_char_index << "\n";
    }

    if (node >= _nodes_info.size()) {
      _nodes_info.resize(node + 1);
    }
    _nodes_info[node] = {substr_len, substr_last_char_index};
    if (_print_debug_info) {
      std::cout << "set node_info -- node: " << node
                << ", substr_len: " << substr_len
                << ", last_char_index: " << substr_last_char_index << "\n";
    }

    if (node >= _edges_from_node.size()) {
      return;
    }

    std::for_each(_edges_from_node[node].begin(), _edges_from_node[node].end(),
                  [&](Edge* edge) {
                    int new_substr_len = substr_len + edge->last_char_index -
                                         edge->first_char_index + 1;
                    _build_str_info_for_nodes(edge->dest_node_index,
                                              new_substr_len,
                                              edge->last_char_index);
                  });
  }

  void _add_prefix(int last_char_index) {
    int last_parent_node = -1;
    int parent_node;
    while (true) {
      parent_node = _active.source_node_index;
      if (_active.is_explicit()) {
        if (_edges.contains(
                {_active.source_node_index, _string[last_char_index]})) {
          break;
        }
      } else {
        Edge e = _edges.at(
            {_active.source_node_index, _string[_active.first_char_index]});
        if (_string[e.first_char_index + _active.size() + 1] ==
            _string[last_char_index]) {
          break;
        }
        parent_node = _split_edge(e, _active);
      }

      _nodes.push_back(Node());
      size_t dest_node_index = _nodes.size() - 1;
      Edge e = Edge(last_char_index, _N, parent_node, dest_node_index);
      _insert_edge(e);
      _leaf_nodes.insert(dest_node_index);

      if (last_parent_node > 0) {
        _nodes[last_parent_node].suffix_node = parent_node;
      }
      last_parent_node = parent_node;

      if (_active.source_node_index == 0) {
        _active.first_char_index += 1;
      } else {
        _active.source_node_index =
            _nodes[_active.source_node_index].suffix_node;
      }
      _canonize_suffix(_active);
    }
    if (last_parent_node > 0) {
      _nodes[last_parent_node].suffix_node = parent_node;
    }
    _active.last_char_index += 1;
    _canonize_suffix(_active);
  }

  void _insert_edge(const Edge& edge) {
    _edges.insert(
        {{edge.source_node_index, _string[edge.first_char_index]}, edge});
  }

  void _remove_edge(const Edge& edge) {
    _edges.erase({edge.source_node_index, _string[edge.first_char_index]});
  }

  int _split_edge(Edge& edge, const Suffix& suffix) {
    _nodes.push_back(Node());
    Edge e = Edge(edge.first_char_index, edge.first_char_index + suffix.size(),
                  suffix.source_node_index, _nodes.size() - 1);
    _remove_edge(edge);
    _insert_edge(e);
    _nodes[e.dest_node_index].suffix_node = suffix.source_node_index;
    edge.first_char_index += suffix.size() + 1;
    edge.source_node_index = e.dest_node_index;
    _insert_edge(edge);
    return e.dest_node_index;
  }

  void _canonize_suffix(Suffix& suffix) {
    if (!suffix.is_explicit()) {
      Edge e = _edges.at(
          {suffix.source_node_index, _string[suffix.first_char_index]});
      if (e.size() <= suffix.size()) {
        suffix.first_char_index += e.size() + 1;
        suffix.source_node_index = e.dest_node_index;
        _canonize_suffix(suffix);
      }
    }
  }

  int find_substring(const std::string& substring) {
    if (substring.empty()) {
      return -1;
    }
    int curr_node = 0;
    int i = 0;
    Edge edge;
    int ln;
    while (i < substring.size()) {
      if (!_edges.contains({curr_node, substring[i]})) {
        return -1;
      }
      edge = _edges.at({curr_node, substring[i]});
      ln = std::min(edge.size() + 1, substring.size() - i);
      for (size_t ii = 0; ii < ln; ++ii) {
        if (substring[i + ii] != _string[edge.first_char_index + ii]) {
          return -1;
        }
      }
      i += edge.size() + 1;
      curr_node = edge.dest_node_index;
    }
    return edge.first_char_index - substring.size() + ln;
  }

  bool has_substring(const std::string& substring) {
    return find_substring(substring) != -1;
  }

  int query_lca(int u, int v) const {
    auto ui = _repr[u];
    auto vi = _repr[v];

    // The RMQ interface uses an exclusive upper bound so we need to go
    // one past that to include the node represented by the upper bound
    // here.
    auto idx =
        (ui <= vi ? _rmq->query(_level.begin() + ui, _level.begin() + vi + 1)
                  : _rmq->query(_level.begin() + vi, _level.begin() + ui + 1));

    // std::cout << "idx for euler: " << idx << "\n";

    // Future work: return the node itself, rather than the node's id?
    return _euler[idx];
  }

  // here we assume that substring is substr of _string of suff tree.
  // ideally this should be checked strictly.
  SuffTreeNodeInfo get_node(const std::string& substring) const {
    SuffTreeNodeInfo res;
    res.sf_ptr = this;
    if (substring.empty()) {
      return res;
    }
    int curr_node = 0;
    int i = 0;
    const Edge* edge;
    size_t ln;
    while (i < substring.size()) {
      if (!_edges.contains({curr_node, substring[i]})) {
        return {.node_id = curr_node};
      }
      edge = &(_edges.at({curr_node, substring[i]}));
      ln = std::min(edge->size() + 1, substring.size() - i);
      for (size_t ii = 0; ii < ln; ++ii) {
        if (substring[i + ii] != _string[edge->first_char_index + ii]) {
          throw std::runtime_error("Not implemented!");
          return {.node_id = curr_node, .edge = edge, .offset = ii};
        }
      }
      if (ln < edge->size() + 1) {
        return {.node_id = curr_node, .edge = edge, .offset = ln};
      }
      i += edge->size() + 1;
      curr_node = edge->dest_node_index;
    }
    return {.node_id = curr_node};
  }

  bool has_symbol_on_next_node(const SuffTreeNodeInfo& info, char c) {
    if (info.edge == nullptr) {
      return _edges.contains({info.node_id, c});
    }
    return _string.at(info.edge->first_char_index + info.offset + 1) == c;
  }

  SuffTreeNodeInfo get_next_node(const SuffTreeNodeInfo& info, char c) const {
    SuffTreeNodeInfo res = {.sf_ptr = this};
    if (c == '$') {
      return res;
    }

    if (info.node_id == -1) {
      return res;
    }

    // Проверяем находимся ли мы в explicit ноде сейчас:
    if (info.edge == nullptr) {
      std::cout << "-------> here3\n";
      // Мы в explicit ноде.
      // Можем ли мы шагнуть из неё по символу `c`:
      if (_edges.contains({info.node_id, c})) {
        const Edge* edge = &(_edges.at({info.node_id, c}));
        std::cout << "-------> here4\n";
        // Проверим, сколько у нас символов на этом ребре, чтобы понять, будем
        // ли мы возвращать explicit или implicit ноду:
        if (edge->size() == 0 || (edge->size() == 1 && edge->leads_to_leaf)) {
          std::cout << "-----> here5\n";
          res.node_id = edge->dest_node_index;
          return res;
        } else {
          std::cout << "-----> here6\n";
          res.node_id = info.node_id;
          res.edge = edge;
          res.offset = 0;
          return res;
        }
      } else {
        // Если не можем шагнуть, то подстрока не найдена, тогда вернём ноду с
        // node_id == -1 (по дефолту там уже это значение стоит у
        // SuffTreeNodeInfo):
        return res;
      }
    } else {
      std::cout << "-------> here7\n";
      // Мы в implicit ноде.
      // Проверим совпадение символов:
      if (_string.at(info.edge->first_char_index + info.offset + 1) == c) {
        std::cout << "-------> here8\n";
        // Теперь надо понять, будет ли следющая нода explicit или implicit.
        if (info.offset + 1 == info.edge->size() ||
            (info.offset + 2 == info.edge->size() &&
             info.edge->leads_to_leaf)) {
          std::cout << "-------> here9\n";
          res.node_id = info.edge->dest_node_index;
          return res;
        } else {
          std::cout << "-------> here10\n";
          res.node_id = info.node_id;
          res.edge = info.edge;
          res.offset = info.offset + 1;
          return res;
        }
      } else {
        return res;
      }
    }
  }

  std::string get_common_prefix_str(int lca_node,
                                    const std::string& substring) const {
    if (substring.empty()) {
      return "";
    }
    int curr_node = 0;
    int i = 0;
    Edge edge;
    int ln;
    while (i < substring.size()) {
      if (curr_node == lca_node) {
        break;
      }

      if (!_edges.contains({curr_node, substring[i]})) {
        return "";
      }
      edge = _edges.at({curr_node, substring[i]});
      ln = std::min(edge.size() + 1, substring.size() - i);
      for (size_t ii = 0; ii < ln; ++ii) {
        if (substring[i + ii] != _string[edge.first_char_index + ii]) {
          return "";
        }
      }
      i += edge.size() + 1;
      curr_node = edge.dest_node_index;
      if (curr_node == lca_node) {
        break;
      }
    }
    return substring.substr(0, i);
  }

  std::string _get_string_by_explicit_node_id(int node_id) const {
    auto [substr_len, last_char_index] = _nodes_info[node_id];
    size_t index = last_char_index - (substr_len - 1);
    return _string.substr(index, substr_len);
  }

  // TODO: tests!
  std::string get_longest_common_prefix(const std::string& s1,
                                        const std::string& s2) const {
    auto node_info1 = get_node(s1);
    auto node_info2 = get_node(s2);

    std::cout << "s1: " << s1 << ", node_info1: {" << node_info1.node_id << ", "
              << node_info1.offset;
    if (node_info1.edge == nullptr) {
      assert(node_info1.offset == 0);
      std::cout << "}\n";
    } else {
      assert(node_info1.offset != 0);
      std::cout << ", " << node_info1.edge->dest_node_index << "}\n";
    }

    std::cout << "s2: " << s2 << ", node_info2: {" << node_info2.node_id << ", "
              << node_info2.offset;
    if (node_info2.edge == nullptr) {
      assert(node_info2.offset == 0);
      std::cout << "}\n";
    } else {
      assert(node_info2.offset != 0);
      std::cout << ", " << node_info2.edge->dest_node_index << "}\n";
    }

    std::string res;
    int lca_node_id;

    if (node_info1.node_id == node_info2.node_id) {
      if (node_info1.edge == node_info2.edge) {
        if (s1.size() < s2.size()) {
          res = s1;
          goto end;
        } else {
          res = s2;
          goto end;
        }
      } else {
        res = s1.substr(0, s1.size() - node_info1.offset);
        goto end;
      }
    }

    lca_node_id = query_lca(node_info1.node_id, node_info2.node_id);
    std::cout << "lca: " << lca_node_id << "\n";

    if (node_info1.node_id == lca_node_id && node_info1.offset != 0) {
      char first_symbol_after_lca_for_s1 =
          _string[node_info1.edge->first_char_index];
      assert(s1.size() >= node_info1.offset);
      int i = s1.size() - node_info1.offset;
      char first_symbol_after_lca_for_s2 = s2.at(i);

      if (first_symbol_after_lca_for_s1 == first_symbol_after_lca_for_s2) {
        res = s1;
        goto end;
      }
      res = s1.substr(0, s1.size() - node_info1.offset);
      goto end;
    }

    if (node_info2.node_id == lca_node_id && node_info2.offset != 0) {
      char first_symbol_after_lca_for_s2 =
          _string[node_info2.edge->first_char_index];
      assert(s2.size() >= node_info2.offset);
      int i = s2.size() - node_info2.offset;
      char first_symbol_after_lca_for_s1 = s1.at(i);

      if (first_symbol_after_lca_for_s2 == first_symbol_after_lca_for_s1) {
        res = s2;
        goto end;
      }
      res = s2.substr(0, s2.size() - node_info2.offset);
      goto end;
    }

    res = _get_string_by_explicit_node_id(lca_node_id);

  end:
    std::cout << "found longest common prefix for s1: `" << s1 << "` and s2: `"
              << s2 << "` -- `" << res << "`\n";
    return res;
  }

  SuffTreeNodeInfo get_root_node() const { return {.node_id = 0}; }

  bool is_marked_node(int node_id) const {
    return _marked_nodes.contains(node_id);
  }
};

std::ostream& operator<<(std::ostream& os, const SuffTreeNodeInfo& info) {
  os << "SuffTreeNodeInfo: <node_id: " << info.node_id
     << ", offset: " << info.offset;
  if (info.edge != nullptr) {
    os << ", edge size: " << info.edge->size() << ", str: "
       << info.sf_ptr->_string.substr(
              info.edge->first_char_index,
              info.edge->last_char_index - info.edge->first_char_index + 1);
  }
  os << ">";
  return os;
}
