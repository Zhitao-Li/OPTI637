cmake_minimum_required( VERSION 2.8 )

project( OPTI637 )

set(CMAKE_CXX_FLAGS "-std=c++0x")

set(Source src/main.cpp)

add_executable(ForwardModel ${Source})

target_link_libraries(ForwardModel -lpthread)

