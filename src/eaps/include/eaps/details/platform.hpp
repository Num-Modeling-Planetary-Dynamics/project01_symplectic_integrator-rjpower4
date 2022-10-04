/// @file eaps/details/platform.hpp
/// @brief Pre-processor macro definitions for exporting platform identification symbols
/// @author Rolfe Power

#ifndef EAPS_CONFIG_PLATFORM_HPP_
#define EAPS_CONFIG_PLATFORM_HPP_

#if defined(__clang__)
#    define EAPS_COMPILER_CLANG
#elif defined(__GNUC__)
#    define EAPS_COMPILER_GNU
#elif defined(_MSC_VER_)
#    define EAPS_COMPILER_MSVC
#else
#    error "Unknown compiler"
#endif

#if defined(__APPLE__)
#    define EAPS_OS_MACOS
#elif defined(__linux__)
#    define EAPS_OS_LINUX
#elif defined(WIN32) || defined(_WIN32) || defined(_WIN64)
#    define EAPS_OS_WINDOWS
#else
#    error "Unknown operating system"
#endif

#if defined(EAPS_OS_MACOS) || defined(EAPS_OS_LINUX)
#    define EAPS_POSIX
#endif

#endif // EAPS_CONFIG_PLATFORM_HPP_