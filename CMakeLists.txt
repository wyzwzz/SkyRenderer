cmake_minimum_required(VERSION 3.17)

project(SkyRenderer LANGUAGES CXX C)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(BUILD_ENABLE_OPENGL ON)
add_subdirectory(ext/CGUtils)

file(
        GLOB_RECURSE SRCS
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.hpp"
)

add_executable(SkyRenderer ${SRCS})

target_link_libraries(SkyRenderer PRIVATE CGUtils)

target_include_directories(SkyRenderer PRIVATE ${PROJECT_SOURCE_DIR}/ext/cy)

target_compile_features(SkyRenderer PRIVATE cxx_std_20)