cmake_minimum_required(VERSION 3.5.1)
# Assembler files (.s) in the source list are ignored by CMake unless we
# “enable” the assembler by telling CMake in the project definition that we’re
# using assembly files.
project(LevMarTest)

add_subdirectory(modules)
