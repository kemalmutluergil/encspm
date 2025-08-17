cmake_minimum_required(VERSION 3.10)

project(EncSPMV)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES "*.cpp")

add_executable(main ${SOURCES})

find_package(SEAL 4.1 REQUIRED) 
target_link_libraries(main SEAL::seal)

include_directories(include)


