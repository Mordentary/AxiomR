
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

set(DISABLE_PROFILING ON CACHE BOOL "DISABLE TRACY PROFILE?")
target_compile_definitions(${PROJECT_NAME} PRIVATE
    $<$<OR:$<CONFIG:Release>>:DISABLE_PROFILING=FALSE>
)
if (DISABLE_PROFILING)
  get_target_property(TRACY_INTERFACE TracyClient INTERFACE_COMPILE_DEFINITIONS)
  list(REMOVE_ITEM TRACY_INTERFACE "TRACY_ENABLE")
  set_target_properties(TracyClient PROPERTIES INTERFACE_COMPILE_DEFINITIONS "${TRACY_INTERFACE}")
endif ()

add_custom_command(TARGET AxiomR POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
        "$<TARGET_FILE_DIR:AxiomR>/assets"  # Create 'assets' directory if it doesn't exist
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/assets" "$<TARGET_FILE_DIR:AxiomR>/assets"
    COMMENT "Copying assets to output directory"
)

# Link stb_image 
target_link_libraries(${PROJECT_NAME}
  PRIVATE 
        glm      
        stb_image
        TracyClient 
)
