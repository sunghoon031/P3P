cmake_minimum_required(VERSION 3.16)
project(p3p_lee_demo LANGUAGES CXX)

add_executable(p3p_lee_demo
  src/main.cpp
  src/p3p.h
)

target_include_directories(p3p_lee_demo PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_features(p3p_lee_demo PRIVATE cxx_std_17)

if(MSVC)
  target_compile_options(p3p_lee_demo PRIVATE /W4 /permissive-)
else()
  target_compile_options(p3p_lee_demo PRIVATE -Wall -Wextra -Wpedantic)
endif()

