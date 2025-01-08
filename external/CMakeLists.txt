cmake_minimum_required(VERSION 3.15)
# Tracy
FetchContent_Declare(
 tracy
 GIT_REPOSITORY https://github.com/wolfpld/tracy.git
 GIT_TAG 5d542dc09f3d9378d005092a4ad446bd405f819a
 GIT_SHALLOW TRUE
 GIT_PROGRESS TRUE
 )
 FetchContent_MakeAvailable(tracy)

 #GLM
add_library(glm INTERFACE)
target_include_directories(glm INTERFACE glm)

# STB Image
add_library(stb_image STATIC ${CMAKE_CURRENT_SOURCE_DIR}/stb_image/stb_image.cpp)
target_include_directories(stb_image PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/stb_image)
