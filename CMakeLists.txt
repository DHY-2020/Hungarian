cmake_minimum_required(VERSION 2.8)
project(Hungarian)

find_package(Eigen3 3.1.0 REQUIRED)


include_directories(
${PROJECT_SOURCE_DIR}/include
${EIGEN3_INCLUDE_DIR}
)

add_executable(Hungarian src/testmain.cpp  src/hungarian_myself.cpp)
target_link_libraries(Hungarian ${EIGEN3_LIBS})