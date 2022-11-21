cmake_minimum_required(VERSION 3.15)
project(hash_cdbg CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
find_package(LibSDSL REQUIRED)

message("LIBSDSL_FOUND: " ${LIBSDSL_FOUND})
message("LIBSDSL_INCLUDE_DIRS: " ${LIBSDSL_INCLUDE_DIRS})
message("LIBSDSL_LIBRARIES: " ${LIBSDSL_LIBRARIES})
message("LIBSDSL_DEFINITIONS: " ${LIBSDSL_DEFINITIONS})

add_library(${PROJECT_NAME} STATIC)

add_subdirectory(src)
target_include_directories(${PROJECT_NAME}
  PUBLIC
    ${LIBSDSL_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBSDSL_LIBRARIES} z pthread)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
target_compile_options(${PROJECT_NAME}
  PUBLIC
    -Wall -Wpedantic -Wextra -msse4.2 -O3 -funroll-loops -fomit-frame-pointer -ffast-math
    $<$<CONFIG:Debug>:-g>
)
target_compile_definitions(${PROJECT_NAME}
  INTERFACE
    $<$<CONFIG:Debug>:_DEBUG>
)

install(TARGETS ${PROJECT_NAME}
  EXPORT
    hash_cdbg-export
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  PUBLIC_HEADER
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}
)
