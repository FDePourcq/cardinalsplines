cmake_minimum_required(VERSION 3.13)
project(splines)

set(CMAKE_CXX_STANDARD 11)

add_executable(cardinalsplines main.cpp)

add_custom_target(gnuplot
        COMMAND ./cardinalsplines 20 2> test.data && gnuplot -e \" plot 'test.data' using 2:3:1 w lines palette \; pause -100;\"
        DEPENDS cardinalsplines)