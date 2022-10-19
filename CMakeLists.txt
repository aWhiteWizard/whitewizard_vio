#whitewizard created at Sep 1st
cmake_minimum_required(VERSION 2.8)
project(whitewizard_vio)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -g")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(Pangolin REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem program_options system) 

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${CERES_INCLUDE_DIRS}
    ${Pangolin_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)
#camera model lib
add_library(camera_model SHARED
    src/camera_models/chessboard/Chessboard.cc
    src/camera_models/calib/CameraCalibration.cc
    src/camera_models/camera/Camera.cc
    src/camera_models/camera/CameraFactory.cc
    src/camera_models/camera/CostFunctionFactory.cc
    src/camera_models/camera/PinholeCamera.cc
    src/camera_models/camera/CataCamera.cc
    src/camera_models/camera/EquidistantCamera.cc
    src/camera_models/camera/ScaramuzzaCamera.cc
    src/camera_models/sparse_graph/Transform.cc
    src/camera_models/gpl/gpl.cc
    src/camera_models/gpl/EigenQuaternionParameterization.cc
)
target_link_libraries(camera_model ${Boost_LIBRARIES} ${OpenCV_LIBS} ${CERES_LIBRARIES})

#vio lib
add_library(vio_lib SHARED
    src/System.cpp
    src/parameter.cpp
    src/estimator.cpp
    src/feature_manager.cpp
    src/feature_track.cpp
    src/utility.cpp

    src/initial/solve_5pts.cpp
    src/initial/initial_aligment.cpp
    src/initial/initial_sfm.cpp
    src/initial/initial_ex_rotation.cpp

    src/backend/vertex.cpp
    src/backend/edge.cpp
    src/backend/solve.cpp
    src/backend/vertex_pose.cpp
    src/backend/edge_reprojection.cpp
    src/backend/edge_imu.cpp
    src/backend/edge_prior.cpp
    src/backend/loss_function.cpp
    src/backend/imu_dataset.cpp
)
target_link_libraries( vio_lib
    ${OpenCV_LIBS}  
    ${CERES_LIBRARIES}
    ${Pangolin_LIBRARIES}
    camera_model
)

#运行EUROC数据集
add_executable(test_euroc app/test_euroc.cpp)
target_link_libraries(test_euroc
    vio_lib
    ${Pangolin_LIBRARIES}
    -lpthread
)