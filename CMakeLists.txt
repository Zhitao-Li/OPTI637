cmake_minimum_required( VERSION 2.8 )

project( OPTI637 )

set(CMAKE_CXX_FLAGS "-std=c++1y")

set(Source src/main.cpp
           src/basicFunc.h
           src/CTFunc.h
           src/MRFunc.h
)

add_executable(ForwardModel ${Source})

target_link_libraries(ForwardModel -lpthread)

