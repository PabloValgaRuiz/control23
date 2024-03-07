cmake_minimum_required(VERSION 2.8)


project(control23_tests)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Bstatic -lpthread")


set(CMAKE_CXX_STANDARD 17)


add_executable(main
    main.cpp
    MC_Dist.cpp
    MC_DistDiaNocheF.cpp)

include_directories(main PUBLIC
    ./
    )
