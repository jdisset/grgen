cmake_minimum_required(VERSION 2.8)

project(grgen_unit CXX)

set(CMAKE_CXX_FLAGS "-O3 -g -std=c++11 -Wall -Wextra -Wundef -Wshadow -Winit-self -Wconversion -Wformat=2 -pedantic ")
file(GLOB SRC
	"*.hpp"
	"*.cpp"
	"../*.h"
	"../*.hpp"
	"../*.cpp"
	)

add_executable(unit_test ${SRC})
