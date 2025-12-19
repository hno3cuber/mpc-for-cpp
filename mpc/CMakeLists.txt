cmake_minimum_required(VERSION 3.10)


project(MPC_Project)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(Eigen3 3.3 REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)


set(SOURCE_FILES
    src/Main.cpp
    src/MPCCal.cpp
    src/GetRefPath.cpp
)


add_executable(mpc_run ${SOURCE_FILES})


target_include_directories(mpc_run PUBLIC
    ${CMAKE_SOURCE_DIR}/include  
    Eigen3::Eigen
    OsqpEigen::OsqpEigen
    ${Python3_INCLUDE_DIRS} 
)


target_link_libraries(mpc_run PUBLIC Eigen3::Eigen OsqpEigen::OsqpEigen ${Python3_LIBRARIES} )


file(COPY ${CMAKE_SOURCE_DIR}/src/ref.csv DESTINATION ${CMAKE_BINARY_DIR})