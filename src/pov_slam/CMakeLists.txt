cmake_minimum_required(VERSION 2.8)
project(pov_slam)

add_definitions(-w)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})


set(CMAKE_CXX_STANDARD 14)

include(CMakePrintHelpers)


find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

find_package(OpenCV 3.4.19 REQUIRED)

find_package(Boost REQUIRED COMPONENTS
    system
    filesystem
)

find_package(Eigen3 3.1.0 REQUIRED)

find_package(PkgConfig)


SET(G2O_USE_OPENMP OFF CACHE BOOL "Build g2o with OpenMP support (EXPERIMENTAL)")
IF(OPENMP_FOUND AND G2O_USE_OPENMP)
  SET (G2O_OPENMP 1)
  SET(g2o_C_FLAGS "${g2o_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(g2o_CXX_FLAGS "${g2o_CXX_FLAGS} -DEIGEN_DONT_PARALLELIZE ${OpenMP_CXX_FLAGS}")
  MESSAGE(STATUS "Compiling with OpenMP support")
ENDIF(OPENMP_FOUND AND G2O_USE_OPENMP)

# activate warnings !!!
SET(g2o_C_FLAGS "${g2o_C_FLAGS} -Wall -W")
SET(g2o_CXX_FLAGS "${g2o_CXX_FLAGS} -Wall -W")

# specifying compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${g2o_CXX_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${g2o_C_FLAGS}")


SET(G2O_CXX_COMPILER "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}")
configure_file(Thirdparty/g2o/g2o/config.h.in Thirdparty/g2o/config.h)

cs_add_library(libg2o
	#types
	Thirdparty/g2o/g2o/types/types_sba.h
	Thirdparty/g2o/g2o/types/types_six_dof_expmap.h
	Thirdparty/g2o/g2o/types/types_sba.cpp
	Thirdparty/g2o/g2o/types/types_six_dof_expmap.cpp
	Thirdparty/g2o/g2o/types/types_seven_dof_expmap.cpp
	Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h
	Thirdparty/g2o/g2o/types/se3quat.h
	Thirdparty/g2o/g2o/types/se3_ops.h
	Thirdparty/g2o/g2o/types/se3_ops.hpp
	
	#added types
	Thirdparty/g2o/g2o/types/dquat2mat.h
	Thirdparty/g2o/g2o/types/dquat2mat.cpp
	Thirdparty/g2o/g2o/types/isometry3d_gradients.h
	Thirdparty/g2o/g2o/types/isometry3d_gradients.cpp
	Thirdparty/g2o/g2o/types/isometry3d_mappings.h
	Thirdparty/g2o/g2o/types/isometry3d_mappings.cpp
	Thirdparty/g2o/g2o/types/vertex_se3.h
	Thirdparty/g2o/g2o/types/vertex_se3.cpp
	Thirdparty/g2o/g2o/types/edge_se3.h
	Thirdparty/g2o/g2o/types/edge_se3.cpp
	
	Thirdparty/g2o/g2o/types/parameter_se3_offset.h
	Thirdparty/g2o/g2o/types/parameter_se3_offset.cpp
	Thirdparty/g2o/g2o/types/vertex_pointxyz.h
	Thirdparty/g2o/g2o/types/vertex_pointxyz.cpp
	Thirdparty/g2o/g2o/types/edge_xyz_prior.h
	Thirdparty/g2o/g2o/types/edge_xyz_prior.cpp
	Thirdparty/g2o/g2o/types/edge_se3_pointxyz.h
	Thirdparty/g2o/g2o/types/edge_se3_pointxyz.cpp
	Thirdparty/g2o/g2o/types/edge_se3_prior.h
	Thirdparty/g2o/g2o/types/edge_se3_prior.cpp
	
	Thirdparty/g2o/g2o/types/edge_se3_mixture.h
	Thirdparty/g2o/g2o/types/edge_se3_mixture.cpp
	Thirdparty/g2o/g2o/types/edge_se3_pointxyz_noninformative.h
	Thirdparty/g2o/g2o/types/edge_se3_pointxyz_noninformative.cpp
	Thirdparty/g2o/g2o/types/edge_se3_pointxyz_noninformative_mixture_vem.h
	Thirdparty/g2o/g2o/types/edge_se3_pointxyz_noninformative_mixture_vem.cpp
	
	Thirdparty/g2o/g2o/types/types_slam.h
	
	#core
	Thirdparty/g2o/g2o/core/base_edge.h
	Thirdparty/g2o/g2o/core/base_binary_edge.h
	Thirdparty/g2o/g2o/core/hyper_graph_action.cpp
	Thirdparty/g2o/g2o/core/base_binary_edge.hpp
	Thirdparty/g2o/g2o/core/hyper_graph_action.h
	Thirdparty/g2o/g2o/core/base_multi_edge.h           
	Thirdparty/g2o/g2o/core/hyper_graph.cpp
	Thirdparty/g2o/g2o/core/base_multi_edge.hpp         
	Thirdparty/g2o/g2o/core/hyper_graph.h
	Thirdparty/g2o/g2o/core/base_unary_edge.h          
	Thirdparty/g2o/g2o/core/linear_solver.h
	Thirdparty/g2o/g2o/core/base_unary_edge.hpp         
	Thirdparty/g2o/g2o/core/marginal_covariance_cholesky.cpp
	Thirdparty/g2o/g2o/core/base_vertex.h               
	Thirdparty/g2o/g2o/core/marginal_covariance_cholesky.h
	Thirdparty/g2o/g2o/core/base_vertex.hpp             
	Thirdparty/g2o/g2o/core/matrix_structure.cpp
	Thirdparty/g2o/g2o/core/batch_stats.cpp             
	Thirdparty/g2o/g2o/core/matrix_structure.h
	Thirdparty/g2o/g2o/core/batch_stats.h               
	Thirdparty/g2o/g2o/core/openmp_mutex.h
	Thirdparty/g2o/g2o/core/block_solver.h              
	Thirdparty/g2o/g2o/core/block_solver.hpp            
	Thirdparty/g2o/g2o/core/parameter.cpp               
	Thirdparty/g2o/g2o/core/parameter.h                 
	Thirdparty/g2o/g2o/core/cache.cpp                   
	Thirdparty/g2o/g2o/core/cache.h
	Thirdparty/g2o/g2o/core/optimizable_graph.cpp       
	Thirdparty/g2o/g2o/core/optimizable_graph.h         
	Thirdparty/g2o/g2o/core/solver.cpp                  
	Thirdparty/g2o/g2o/core/solver.h
	Thirdparty/g2o/g2o/core/creators.h                 
	Thirdparty/g2o/g2o/core/optimization_algorithm_factory.cpp
	Thirdparty/g2o/g2o/core/estimate_propagator.cpp     
	Thirdparty/g2o/g2o/core/optimization_algorithm_factory.h
	Thirdparty/g2o/g2o/core/estimate_propagator.h       
	Thirdparty/g2o/g2o/core/factory.cpp                 
	Thirdparty/g2o/g2o/core/optimization_algorithm_property.h
	Thirdparty/g2o/g2o/core/factory.h                   
	Thirdparty/g2o/g2o/core/sparse_block_matrix.h
	Thirdparty/g2o/g2o/core/sparse_optimizer.cpp  
	Thirdparty/g2o/g2o/core/sparse_block_matrix.hpp
	Thirdparty/g2o/g2o/core/sparse_optimizer.h
	Thirdparty/g2o/g2o/core/hyper_dijkstra.cpp 
	Thirdparty/g2o/g2o/core/hyper_dijkstra.h
	Thirdparty/g2o/g2o/core/parameter_container.cpp     
	Thirdparty/g2o/g2o/core/parameter_container.h
	Thirdparty/g2o/g2o/core/optimization_algorithm.cpp 
	Thirdparty/g2o/g2o/core/optimization_algorithm.h
	Thirdparty/g2o/g2o/core/optimization_algorithm_with_hessian.cpp 
	Thirdparty/g2o/g2o/core/optimization_algorithm_with_hessian.h
	Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.cpp 
	Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h
	Thirdparty/g2o/g2o/core/optimization_algorithm_gauss_newton.cpp 
	Thirdparty/g2o/g2o/core/optimization_algorithm_gauss_newton.h
	Thirdparty/g2o/g2o/core/jacobian_workspace.cpp 
	Thirdparty/g2o/g2o/core/jacobian_workspace.h
	Thirdparty/g2o/g2o/core/robust_kernel.cpp 
	Thirdparty/g2o/g2o/core/robust_kernel.h
	Thirdparty/g2o/g2o/core/robust_kernel_factory.cpp
	Thirdparty/g2o/g2o/core/robust_kernel_factory.h
	Thirdparty/g2o/g2o/core/robust_kernel_impl.cpp 
	Thirdparty/g2o/g2o/core/robust_kernel_impl.h
	#stuff
	Thirdparty/g2o/g2o/stuff/string_tools.h
	Thirdparty/g2o/g2o/stuff/color_macros.h 
	Thirdparty/g2o/g2o/stuff/macros.h
	Thirdparty/g2o/g2o/stuff/timeutil.cpp
	Thirdparty/g2o/g2o/stuff/misc.h
	Thirdparty/g2o/g2o/stuff/timeutil.h
	Thirdparty/g2o/g2o/stuff/os_specific.c    
	Thirdparty/g2o/g2o/stuff/os_specific.h
	Thirdparty/g2o/g2o/stuff/string_tools.cpp
	Thirdparty/g2o/g2o/stuff/property.cpp       
	Thirdparty/g2o/g2o/stuff/property.h       
)

include_directories(
	include
	Thirdparty/g2o
	${catkin_INCLUDE_DIRS}
	${OpenCV_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
)


set(POV_SLAM_SRC
	src/types/kalman_filter.cpp
	src/types/Pose2.cpp
)

cs_add_library(povslam_lib
    ${POV_SLAM_SRC}
)
target_link_libraries(povslam_lib
    ${catkin_LIBRARIES}
    ${OpenCV_LIBS}
)

cs_add_executable(pov_slam_demo demo/demo_dynamic.cpp)
target_link_libraries(pov_slam_demo
	povslam_lib
	libg2o
	GLEW
	boost_serialization
	crypto
)

cs_install()
cs_export()


















