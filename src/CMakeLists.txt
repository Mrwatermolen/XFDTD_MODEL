add_executable(fdtd_model_msh main.cpp)

target_link_libraries(fdtd_model_msh PRIVATE xfdtd::xfdtd_core)
target_include_directories(fdtd_model_msh PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(model_scatter model_scatter.cpp)

target_link_libraries(model_scatter PRIVATE xfdtd::xfdtd_core)
target_include_directories(model_scatter PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
