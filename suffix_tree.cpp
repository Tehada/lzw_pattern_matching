#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
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

  SuffixTree(const std::string& str) : _active(0, 0, -1), _nodes(1) {
    _string = str + "$";
    _N = _string.size() - 1;
    for (int i = 0; i < _string.size(); i++) {
      _add_prefix(i);
    }
    for (auto& [key, edge] : _edges) {
      if (key.first >= _nodes_childs.size()) {
        _nodes_childs.resize(key.first + 1);
      }
      _nodes_childs[key.first].insert(edge.dest_node_index);
      std::cout << "---> " << key.first << " <--> " << edge.dest_node_index
                << ", suffix_node: " << _nodes[key.first].suffix_node << "\n";
    }

    _build_euler_path(0, 0, std::back_inserter(_euler),
                      std::back_inserter(_level));

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

    // std::cout << "idx for euler: " << idx << "\n";

    // Future work: return the node itself, rather than the node's id?
    return _euler[idx];
  }

  int get_node(const std::string& substring) const {
    if (substring.empty()) {
      return -1;
    }
    int curr_node = 0;
    int i = 0;
    Edge edge;
    int ln;
    while (i < substring.size()) {
      if (!_edges.contains({curr_node, substring[i]})) {
        return curr_node;
      }
      edge = _edges.at({curr_node, substring[i]});
      ln = std::min(edge.size() + 1, substring.size() - i);
      for (size_t ii = 0; ii < ln; ++ii) {
        if (substring[i + ii] != _string[edge.first_char_index + ii]) {
          return curr_node;
        }
      }
      i += edge.size() + 1;
      curr_node = edge.dest_node_index;
    }
    return curr_node;
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

  std::string get_longest_common_prefix(const std::string& s1,
                                        const std::string& s2) const {
    int node1 = get_node(s1);
    int node2 = get_node(s2);
    int lca_node = query_lca(node1, node2);

    std::cout << "lca: " << lca_node << "\n";

    // TODO: probably I should just precalculate info about corresponding string
    // for every node, instead of doing this traversal.
    std::string str = get_common_prefix_str(lca_node, s1);
    std::cout << "found longest common prefix: `" << str << "`\n";
    return str;
  }
};

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

bool lemma28(const std::string& pattern, const std::string& prefix_snippet,
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
    return false;
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
    return true;
  }

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
    return true;
  }

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
    return false;
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
      std::cout << "Not implemented3\n";
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
        return true;

      } else {
        std::cout << "pattern `" << pattern << "` doesn't occur inside concat `"
                  << concat << "`!\n";
      }
    } else {
      std::cout << "Not implemented1\n";
    }
  } else {
    std::cout << "Not implemented2\n";
  }

  return true;
}

int main() {
  // // test1
  // std::string pattern = "abcabfa";
  // std::string prefix_snippet = pattern.substr(0, 5);
  // std::string suffix_snippet = pattern.substr(2);

  // test 1.5
  // std::string pattern = "abacaba";

  // test2
  std::string pattern = "abcabfa";
  std::string prefix_snippet = pattern.substr(0, 5);
  std::string suffix_snippet = pattern.substr(3);

  // // test3 -- passing, but we are not returning the first occurence.
  // std::string pattern = "ababab";
  // std::string prefix_snippet = pattern.substr(0, 5);
  // std::string suffix_snippet = pattern.substr(1);

  std::string reversed_pattern(pattern.rbegin(), pattern.rend());

  std::vector<size_t> lps(pattern.size());
  fill_lps(pattern, &lps);

  std::vector<size_t> reversed_lps(reversed_pattern.size());
  fill_lps(reversed_pattern, &reversed_lps);

  auto st = SuffixTree(pattern);
  auto reversed_st = SuffixTree(reversed_pattern);

  bool res = lemma28(pattern, prefix_snippet, suffix_snippet, lps, reversed_lps,
                     st, reversed_st);
  std::cout << res << "\n";

  // std::cout << st.has_substring("aba") << "\n";
  // std::cout << st.has_substring("abb") << "\n";
  // std::cout << st.has_substring("abaca") << "\n";
  // std::cout << st.has_substring("aca") << "\n";
  // std::cout << st.has_substring("aba$") << "\n";
  // std::string s1 = "aba", s2 = "aca";
}
