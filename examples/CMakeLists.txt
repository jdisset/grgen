cmake_minimum_required (VERSION 2.6)
project (ViewerExample)
set(CMAKE_CXX_FLAGS "-O3 -g -std=c++14 -Wall -Wextra -Wundef -Wshadow -Winit-self -Wconversion -Wformat=2 -pedantic ")
file(GLOB SRC
	"*.hpp"
	"*.cpp"
	"../*.h"
	"../*.hpp"
	"../*.cpp"
	)

add_executable(grn-example ${SRC})
