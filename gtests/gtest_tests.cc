#include <gtest/gtest.h>

TEST(GTest, working)
{
  unsigned int i = 1;
  ++i;
  EXPECT_EQ((unsigned int)2, i);
}
