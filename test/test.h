#ifndef LP_MP_TEST_H
#define LP_MP_TEST_H

#include <stdexcept>
#include <string>

inline void test(const bool& pred)
{
  if (!pred)
    throw std::runtime_error("Test failed.");
}

#endif // LP_MP_TEST_H
