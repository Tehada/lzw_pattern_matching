#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "pm_rmq.hpp"

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

  size_t size() const { return last_char_index - first_char_index; }
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

template <typename rmq_impl>
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
  std::unique_ptr<rmq_impl> _rmq;

  std::vector<int> _repr;

  SuffixTree(const std::string& str)
      : _string(str), _active(0, 0, -1), _N(str.size() - 1), _nodes(1) {
    for (int i = 0; i < _string.size(); i++) {
      _add_prefix(i);
    }
    for (auto& [key, edge] : _edges) {
      if (key.first >= _nodes_childs.size()) {
        _nodes_childs.resize(key.first + 1);
      }
      _nodes_childs[key.first].insert(edge.dest_node_index);
      std::cout << "---> " << key.first << " <--> " << edge.dest_node_index
                << "\n";
    }

    _build_euler_path(0, 0, std::back_inserter(_euler),
                      std::back_inserter(_level));

    std::cout << "euler:\n";
    for (auto elem : _euler) {
      std::cout << elem << " ";
    }
    std::cout << "\n";

    std::cout << "level:\n";
    for (auto elem : _level) {
      std::cout << elem << " ";
    }
    std::cout << "\n";

    _rmq.reset(new rmq_impl(_level.begin(), _level.end()));
  }

  template <typename OutputIterator1, typename OutputIterator2>
  void _build_euler_path(int node, ssize_t level, OutputIterator1 eulerit,
                         OutputIterator2 levelit) {
    *eulerit++ = node;
    *levelit++ = level;
    if (node >= _nodes_childs.size()) {
      return;
    }
    std::cout << "-- node: " << node << "\n-- childs: ";
    for (auto elem : _nodes_childs[node]) {
      std::cout << elem << " ";
    }
    std::cout << "\n";

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
      Edge e = Edge(last_char_index, _N, parent_node, _nodes.size() - 1);
      _insert_edge(e);

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

    std::cout << "idx for euler: " << idx << "\n";

    // Future work: return the node itself, rather than the node's id?
    return _euler[idx];
  }
};

int main() {
  auto st =
      SuffixTree<pm_rmq<std::vector<ssize_t>::const_iterator>>("abacaba$");

  int lca = st.query_lca(4, 6);
  std::cout << "lca: " << lca << "\n";

  // std::cout << st.has_substring("aba") << "\n";
  // std::cout << st.has_substring("abb") << "\n";
  // std::cout << st.has_substring("abaca") << "\n";
  // std::cout << st.has_substring("aca") << "\n";
  // std::cout << st.has_substring("aba$") << "\n";
  // std::string s1 = "aba", s2 = "aca";
}
