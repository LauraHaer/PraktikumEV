#include <gtest/gtest.h>

TEST(AsIlike, whatever)
{
  unsigned int i = 1;
  ++i;
  EXPECT_EQ((unsigned int)2, i);
}


//FIXME!!!
TEST(AsIlike, DISABLED_failure)
{
  unsigned int i = 1;
  ++i;
  EXPECT_EQ((unsigned int)3, i);
}
