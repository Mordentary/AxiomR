
# Collect all source files
file(GLOB_RECURSE SOURCES 
    "*.cpp"
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} "windows_bitmap.cpp" "vec.cpp")

# Add include directories
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${SDL2_INCLUDE_DIRS}
)

# Link stb_image 
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        stb_image
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/external/stb
)

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)