add_library(
    eaps 
    body/celestial_body.cpp
    twobody/conic.cpp
    twobody/kepler.cpp
    twobody/keplerian_elements.cpp
)
add_library(EAPS::eaps ALIAS eaps)
target_link_libraries(eaps PUBLIC eaps-headers)
target_compile_features(eaps PRIVATE cxx_std_17)