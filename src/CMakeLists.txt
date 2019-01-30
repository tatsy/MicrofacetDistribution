set(BUILD_TARGET "microfacet")
add_definitions("-D_CRT_SECURE_NO_WARNINGS")

file(GLOB SOURCE_FILES "*.cpp" "*.h")
add_executable(${BUILD_TARGET} ${SOURCE_FILES})
set_target_properties(${BUILD_TARGET} PROPERTIES DEBUG_POSTFIX "-debug")

source_group("Source Files" FILES ${SOURCE_FILES})

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
  set_property(TARGET ${BUILD_TARGET} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
endif()
