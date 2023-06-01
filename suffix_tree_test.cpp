#include "suffix_tree.hpp"

#include <gtest/gtest.h>

#include <string>

TEST(SuffixTree, BasicChecks) {
  std::string pattern = "abacaba";
  SuffixTree st(pattern, true);

  EXPECT_TRUE(st.has_substring("baca"));
  EXPECT_FALSE(st.has_substring("bacad"));
  EXPECT_TRUE(st.has_substring("abacaba"));
  // EXPECT_FALSE(st.has_substring("abacaba$"));
}
