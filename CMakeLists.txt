cmake_minimum_required(VERSION 3.9)
project(huawei)

set(CMAKE_CXX_STANDARD 11)

add_executable(huawei
        ARMAModel.h ARModel.h MAModel.h ARIMAModel.h  main.cpp)
