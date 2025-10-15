cmake_minimum_required(VERSION 3.28)
project(OPTYMALIZATION)

set(CMAKE_CXX_STANDARD 20)

add_executable(OPTYMALIZATION main.cpp
        opt_alg.cpp
        matrix.cpp
        ode_solver.cpp
        solution.cpp
        user_funs.cpp)
