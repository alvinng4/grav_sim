cmake_minimum_required(VERSION 3.14)

add_library(pcg STATIC
    src/pcg_basic.c
)

add_library(pcg::pcg ALIAS pcg)

target_include_directories(pcg
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties(pcg PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS ON
    POSITION_INDEPENDENT_CODE ON
)
