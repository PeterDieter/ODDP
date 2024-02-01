cmake_minimum_required(VERSION 3.5)

project(onlineAssignment)

set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
# List all source files
set(SRC_FILES
    src/main.cpp
    src/Environment.cpp
    src/WQAssign.cpp
    src/Data.cpp
)

# List all header files
set(HDR_FILES
    src/Environment.h
    src/Data.h
    src/Matrix.h
    src/WQAssign.h
    src/xorshift128.h
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# Set the path to Gurobi (modify this according to your installation)
find_package(GUROBI REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})

if(CXX)
  set(CMAKE_CXX_STANDARD 11)
  target_link_libraries(${onlineAssignment} optimized ${GUROBI_CXX_LIBRARY}
                        debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()


# Create an executable target
add_executable(onlineAssignment ${SRC_FILES} ${HDR_FILES})
set_property(TARGET onlineAssignment PROPERTY CXX_STANDARD 14)

target_link_libraries(onlineAssignment ${GUROBI_LIBRARIES})




