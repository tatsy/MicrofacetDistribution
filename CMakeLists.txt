cmake_minimum_required(VERSION 2.8.0)
project(Microfacet)

if (NOT MSVC)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11")
endif()

if (NOT APPLE)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()


set(BUILD_TARGET "microfacet")
add_definitions("-D_CRT_SECURE_NO_WARNINGS")
file(GLOB SOURCE_FILES "*.cpp" "*.h")
add_executable(${BUILD_TARGET} ${SOURCE_FILES})
source_group("Source Files" FILES ${SOURCE_FILES})

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
  set_property(TARGET ${BUILD_TARGET} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
endif()
