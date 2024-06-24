cmake_minimum_required(VERSION 3.0.0)
 
project(HeatExchanger)
 
set(CMAKE_BUILD_TYPE Debug)
 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++20 -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++20") -lnlopt
 
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VTK COMPONENTS 
  IOXML
)
if (NOT VTK_FOUND)
  message(FATAL_ERROR "WriteVTU: Unable to find the VTK build folder.")
endif()

set(NLOPT_SOURCE_DIR "/usr/local/lib/")  # Adjust path if needed
set(NLOPT_LIBRARY_DIR "/usr/local/include")  # Adjust path accordingly

find_library(NLOPT_LIBRARY NAMES nlopt PATHS ${NLOPT_LIBRARY_DIR})
#add_subdirectory(${NLOPT_SOURCE_DIR} EXCLUDE_FROM_ALL)
 
#find_package(NLOPT)

#include_directories("${PROJECT_SOURCE_DIR}/include/")
 
set (source_dir "${PROJECT_SOURCE_DIR}/src/")
#file (GLOB source_files "${source_dir}/*.cpp")

set(source_files
    "${source_dir}/main.cpp"
    "${source_dir}/MeshHelperFunctions.cpp"
    "${source_dir}/FowardAnalysisFunctions.cpp"
    "${source_dir}/ProjectionHelperFunctions.cpp"
    "${source_dir}/FeaHelperFunctions.cpp"
    "${source_dir}/ThermalFeaFunction.cpp"
    "${source_dir}/ElasticFeaFunction.cpp"
    "${source_dir}/WriteVtuFunction.cpp"
    # ... other source files ...
)

# Path to Eigen3 library
#set(EIGEN3_INCLUDE_DIR "C:/Users/wbello3/toolbox/eigen-3.4.0/")
set(EIGEN3_INCLUDE_DIR "/mnt/c/Users/wbello3/toolbox/eigen-3.4.0/")

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")
#add_executable(HeatExchanger MACOSX_BUNDLE WriteVTU.cpp )



find_package(OpenMP REQUIRED)
add_executable (HeatExchanger ${source_files})

target_include_directories(HeatExchanger PRIVATE ${EIGEN3_INCLUDE_DIR})
target_link_libraries(HeatExchanger PRIVATE ${VTK_LIBRARIES} ${NLOPT_LIBRARY} nlopt PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(HeatExchanger PUBLIC "${OpenMP_CXX_FLAGS}")
message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")


#target_link_libraries (minifem PUBLIC Eigen3::Eigen)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS HeatExchanger
  MODULES ${VTK_LIBRARIES}
)