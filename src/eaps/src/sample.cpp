#include <iostream>

#include "eaps/config.hpp"

int main()
{
#if defined(EAPS_OS_LINUX)
    std::cout << "LINUX" << std::endl;
#endif
}