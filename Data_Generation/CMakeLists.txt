cmake_minimum_required(VERSION 3.7)
project(MoranDrive)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        Parameters.cpp
        Parameters.hpp
        Population.cpp
        Population.hpp
        Simulation.cpp
        Simulation.hpp StructPop.cpp StructPop.h StructPopSim.cpp StructPopSim.h main.cpp)

add_executable(MoranDrive ${SOURCE_FILES})