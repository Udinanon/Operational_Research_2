cmake_minimum_required(VERSION 3.0.0)
project(tsp VERSION 0.1.0)

#include(CTest)
#enable_testing()

include_directories("/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include/ilcplex")
include_directories("opt/ibm/ILOG/CPLEX_Studio2211/cplex/include/ilcplex")
include_directories("/home/lorenzo/concorde")
link_directories("/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic")
link_directories("/opt/ibm/ILOG/CPLEX_Studio2211/concert/lib/x86-64_linux/static_pic")
link_directories("/home/lorenzo/concorde")

add_executable(tsp tsp_cplex.c tsp_meta.c tsp_generative.c tsp.c utils.c io_utils.c chrono.c main.c)

#set(CPACK_PROJECT_NAME ${PROJECT_NAME})
#set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(CMAKE_CXX_FLAG "-Wall -O3")
add_compile_options(-Wall -O3)
target_link_libraries(${PROJECT_NAME} m -lcplex -lpthread -ldl -lconcorde)
#include(CPack)
