/// @file eaps/config.hpp
/// @brief configuration amalgamation header; included by most other headers
/// @author Rolfe Power

#ifndef EAPS_CONFIG_HPP_
#define EAPS_CONFIG_HPP_

// ---------------------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------------------

#include "eaps/details/platform.hpp"

#include <cstdint>

// ---------------------------------------------------------------------------------------
// Windows Symbol Export
// ---------------------------------------------------------------------------------------
#if defined(EAPS_OS_WINDOWS) && defined(EAPS_BUILT_AS_SHARED)
#    if defined(EAPS_EXPORT)
#        define EAPS_API __declspec(dllexport)
#    else
#        define EAPS_API __declspec(dllimport)
#    endif
#endif

#if !defined(EAPS_API)
#    define EAPS_API
#endif

// ---------------------------------------------------------------------------------------
// Library-wide type aliases
// ---------------------------------------------------------------------------------------
using NaifId = std::int32_t;


#endif // EAPS_CONFIG_HPP_