set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

set(PROJECT_SOURCES source/model.cpp source/solver.cpp source/plotter.cpp)
add_library(solver SHARED ${PROJECT_SOURCES})
target_include_directories(solver PUBLIC include/)
target_link_libraries(solver Threads::Threads)

add_executable(main main.cpp)
target_include_directories(main PUBLIC "${PROJECT_FOLDER}/project/include/")
target_link_libraries(main solver)
