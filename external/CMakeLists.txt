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

 # glm
FetchContent_Declare(
  glm
  GIT_REPOSITORY https://github.com/g-truc/glm.git
  GIT_TAG        0.9.9.8 # Or a specific branch/tag you need
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
)

# STB Image
add_library(stb_image STATIC stb_image/stb_image.cpp)
target_include_directories(stb_image 
    PUBLIC stb_image
    INTERFACE stb_image
)