set(MODULES_DIR "${CMAKE_CURRENT_LIST_DIR}")

# Don't let levmar build it's demo program
set(BUILD_DEMO OFF CACHE BOOL "Build demo program?" FORCE)
add_subdirectory(levmar)

add_subdirectory(tests)
