cmake_minimum_required(VERSION 3.8)
project(ICDE_kcore)

set(Boost_NO_BOOST_CMAKE ON)
set(CMAKE_CXX_STANDARD 14)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

#set(BOOST_ROOT ../boost_1_75_0)
find_package(Boost REQUIRED  unit_test_framework filesystem system date_time program_options)
#find_package(Threads)


include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

set(CMAKE_CXX_STANDARD 14)

add_executable(ICDE_kcore main.cpp src/graph.cpp src/config.cpp src/mylib.cpp src/myalgo.cpp src/anchoredkcore.cpp)

target_link_libraries(ICDE_kcore pthread boost_serialization boost_filesystem boost_system)
