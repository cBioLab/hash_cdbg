cmake_policy(SET CMP0091 NEW)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
list(FILTER CPP_FILES EXCLUDE REGEX ".*boss_tests.cpp$")
file(GLOB HPP_FILES ${PROJECT_SOURCE_DIR}/include/*.h*)
target_sources(${PROJECT_NAME}
  PUBLIC
    ${CPP_FILES}
)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HPP_FILES}")
