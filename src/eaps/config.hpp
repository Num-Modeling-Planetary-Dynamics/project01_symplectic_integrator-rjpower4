/// @file eaps/config.hpp
/// @brief library-wide pre-processor configurations
/// @author Rolfe Power

#ifndef EAPS_CONFIG_HPP_
#define EAPS_CONFIG_HPP_

#include <cstdint>

// ---------------------------------------------------------------------------------------
// Compiler Identification
// ---------------------------------------------------------------------------------------
#if defined(__clang__)
#    define EAPS_CLANG
#elif defined(__GNUC__)
#    define EAPS_GCC
#elif defined(_MSC_VER)
#    define EAPS_MSVC
#else
#    error Unknown compiler
#endif

// ---------------------------------------------------------------------------------------
// Operating System Identification
// ---------------------------------------------------------------------------------------
#if defined(__APPLE__)
#    define EAPS_OSX
#elif defined(__linux__)
#    define EAPS_LINUX
#elif defined(WIN32) || defined(_WIN32) || defined(_WIN64)
#    define EAPS_WINDOWS
#else
#    error Unknown platform
#endif

#if defined(EAPS_LINUX) || defined(EAPS_OSX)
#    define EAPS_POSIX
#endif

// ---------------------------------------------------------------------------------------
// Windows Symbol Export
// ---------------------------------------------------------------------------------------
#if defined(EAPS_WINDOWS) && defined(EAPS_BUILT_AS_SHARED)
#    if defined(EAPS_EXPORT)
#        define EAPS_API __declspec(dllexport)
#    else
#        define EAPS_API __declspec(dllimport)
#    endif
#else
#    define EAPS_API
#endif

// ---------------------------------------------------------------------------------------
// Type Aliases
// ---------------------------------------------------------------------------------------
namespace eaps
{

using NaifId = std::int32_t;

} // namespace eaps

#endif // EAPS_CONFIG_HPP_