cmake_minimum_required(VERSION 3.14)

project(grav_sim
    VERSION 1.0.1
    LANGUAGES C
    DESCRIPTION "N-body gravity simulation library with C and Python API"
)

### Dependencies & Modules ###
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# PROJECT_IS_TOP_LEVEL was introduced after CMake 3.14, so use a compatible
# equivalent.
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(GRAV_SIM_IS_TOP_LEVEL ON)
else()
    set(GRAV_SIM_IS_TOP_LEVEL OFF)
endif()

# Non-debug builds used optimization.
# Do not modify the parent project's build type when used as a subdirectory.
if(GRAV_SIM_IS_TOP_LEVEL
   AND NOT CMAKE_CONFIGURATION_TYPES
   AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the build type" FORCE
    )
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel
    )
endif()

### Options ###
option(USE_OPENMP "Enable OpenMP support" OFF)
option(USE_HDF5 "Enable HDF5 support" OFF)
option(USE_FFTW3 "Enable FFTW3 support" OFF)
option(ENABLE_SANITIZERS
    "Enable address and undefined-behavior sanitizers"
    OFF
)
option(GRAV_SIM_INSTALL
    "Generate installation and CMake package rules"
    ${GRAV_SIM_IS_TOP_LEVEL}
)

# Backward compatibility with the old build interface.
option(DEBUG
    "Legacy debug mode; prefer CMAKE_BUILD_TYPE=Debug and ENABLE_SANITIZERS=ON"
    OFF
)

if(DEBUG)
    message(DEPRECATION
        "The DEBUG option is deprecated. Prefer "
        "-DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZERS=ON"
    )
    set(ENABLE_SANITIZERS ON)
endif()

### Sources ###
set(GRAV_SIM_SOURCES
    src/acceleration.c
    src/acceleration_barnes_hut.c
    src/cosmology.c
    src/error.c
    src/grav_sim.c
    src/integrator.c
    src/integrator_rk_embedded.c
    src/integrator_ias15.c
    src/integrator_whfast.c
    src/linear_octree.c
    src/math_functions.c
    src/output.c
    src/progress_bar.c
    src/python_interface.c
    src/settings.c
    src/system.c
    src/utils.c
)

set(GRAV_SIM_FFTW3_SOURCES
    src/acceleration_PM.c
    src/integrator_cosmology_leapfrog.c
)

### Library ###
add_library(grav_sim SHARED ${GRAV_SIM_SOURCES})
add_library(grav_sim::grav_sim ALIAS grav_sim)

### Standard ###
set_target_properties(grav_sim PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS ON
    POSITION_INDEPENDENT_CODE ON
)

### Includes ###
target_include_directories(grav_sim
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src/internal"
)

### Warnings ###
if(MSVC)
    target_compile_options(grav_sim PRIVATE /W4)
elseif(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(grav_sim PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wpointer-arith
        -Wcast-qual
        -Wstrict-prototypes
        -Wmissing-prototypes
    )
endif()

### Legacy Debug Behavior ###
if(DEBUG AND NOT MSVC)
    target_compile_options(grav_sim PRIVATE -O0 -g)
endif()

### Sanitizers ###
if(ENABLE_SANITIZERS)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(grav_sim PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
        )
        target_link_options(grav_sim PRIVATE
            -fsanitize=address,undefined
        )
    else()
        message(WARNING
            "ENABLE_SANITIZERS is ON, but address/undefined sanitizers "
            "are not configured for compiler '${CMAKE_C_COMPILER_ID}'"
        )
    endif()
endif()

### Version Definition ###
set(GRAV_SIM_VERSION_INFO "${PROJECT_VERSION}")

target_compile_definitions(grav_sim PRIVATE
    VERSION_INFO="${GRAV_SIM_VERSION_INFO}"
)

### Bundled PCG ###
if(NOT TARGET pcg)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/src/libs/pcg"
        "${CMAKE_CURRENT_BINARY_DIR}/_deps/pcg"
    )
endif()

if(NOT TARGET pcg)
    message(FATAL_ERROR
        "src/libs/pcg did not define the expected 'pcg' target"
    )
endif()

target_link_libraries(grav_sim PUBLIC pcg)

### OpenMP ###
if(USE_OPENMP)
    find_package(OpenMP QUIET COMPONENTS C)

    if(OpenMP_C_FOUND)
        target_link_libraries(grav_sim PRIVATE OpenMP::OpenMP_C)
        target_compile_definitions(grav_sim PRIVATE USE_OPENMP)
        message(STATUS "Building grav_sim with OpenMP support")
    else()
        message(WARNING
            "OpenMP support was requested, but OpenMP was not found"
        )
    endif()
endif()

### HDF5 ###
if(USE_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS C)

    target_compile_definitions(grav_sim PRIVATE USE_HDF5)

    # Imported HDF5 targets are preferred when available, but the variable
    # fallback keeps this compatible with older CMake/HDF5 packages.
    if(TARGET HDF5::HDF5)
        target_link_libraries(grav_sim PRIVATE HDF5::HDF5)
    elseif(TARGET hdf5::hdf5)
        target_link_libraries(grav_sim PRIVATE hdf5::hdf5)
    else()
        target_include_directories(grav_sim PRIVATE ${HDF5_INCLUDE_DIRS})
        target_link_libraries(grav_sim PRIVATE ${HDF5_LIBRARIES})

        if(HDF5_DEFINITIONS)
            target_compile_options(grav_sim PRIVATE ${HDF5_DEFINITIONS})
        endif()
    endif()

    message(STATUS "Building grav_sim with HDF5 support")
endif()

### FFTW3 ###
if(USE_FFTW3)
    find_package(FFTW3 QUIET)

    if(FFTW3_FOUND)
        target_sources(grav_sim PRIVATE ${GRAV_SIM_FFTW3_SOURCES})
        target_compile_definitions(grav_sim PRIVATE USE_FFTW3)

        # Support common imported-target names as well as the variables
        # exposed by the existing custom FindFFTW3.cmake module.
        if(TARGET FFTW3::fftw3)
            target_link_libraries(grav_sim PRIVATE FFTW3::fftw3)
        elseif(TARGET FFTW3::FFTW3)
            target_link_libraries(grav_sim PRIVATE FFTW3::FFTW3)
        else()
            target_include_directories(grav_sim PRIVATE
                ${FFTW3_INCLUDE_DIRS}
            )
            target_link_libraries(grav_sim PRIVATE
                ${FFTW3_LIBRARIES}
            )
        endif()

        message(STATUS "Building grav_sim with FFTW3 support")
    else()
        message(WARNING
            "FFTW3 support was requested, but FFTW3 was not found"
        )
    endif()
endif()

### Platform Libraries ###
if(NOT MSVC)
    target_link_libraries(grav_sim PRIVATE m)
endif()

if(UNIX AND NOT APPLE)
    target_link_libraries(grav_sim PRIVATE rt)
endif()

### Installation ###
if(GRAV_SIM_INSTALL)
    install(TARGETS
            grav_sim
            pcg
        EXPORT grav_simTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # grav_sim public headers
    install(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )

    # PCG is exposed by grav_sim's public API, so its header must be installed.
    install(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/libs/pcg/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )

    ### Package Config ###
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/grav_simConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}"
        COMPATIBILITY SameMajorVersion
    )

    file(WRITE
        "${CMAKE_CURRENT_BINARY_DIR}/grav_simConfig.cmake.in"
        "@PACKAGE_INIT@\n\n"
        "include(\"\${CMAKE_CURRENT_LIST_DIR}/grav_simTargets.cmake\")\n\n"
        "check_required_components(grav_sim)\n"
    )

    configure_package_config_file(
        "${CMAKE_CURRENT_BINARY_DIR}/grav_simConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/grav_simConfig.cmake"
        INSTALL_DESTINATION
            "${CMAKE_INSTALL_LIBDIR}/cmake/grav_sim"
    )

    install(EXPORT grav_simTargets
        FILE grav_simTargets.cmake
        NAMESPACE grav_sim::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/grav_sim"
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/grav_simConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/grav_simConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/grav_sim"
    )
endif()

### Testing ###
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    enable_testing()
    add_subdirectory(tests)
endif()
