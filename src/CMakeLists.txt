
# Collect all source files
file(GLOB_RECURSE SOURCES 
    "*.cpp"
)



# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} "windows_bitmap.cpp" "math.cpp" "framebuffer.cpp" "texture.cpp" "window.cpp")


# Add include directories
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${glm_SOURCE_DIR}
)



add_custom_command(TARGET AxiomR POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
        "$<TARGET_FILE_DIR:AxiomR>/assets"  # Create 'assets' directory if it doesn't exist
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/assets" "$<TARGET_FILE_DIR:AxiomR>/assets"
    COMMENT "Copying assets to output directory"
)

# Link stb_image 
target_link_libraries(${PROJECT_NAME}
   PUBLIC 
        glm      
        stb_image
        TracyClient
)
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/external/stb
)

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)