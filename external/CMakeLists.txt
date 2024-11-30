# STB Image
add_library(stb_image STATIC stb_image/stb_image.cpp)
target_include_directories(stb_image 
    PUBLIC stb_image
    INTERFACE stb_image
)