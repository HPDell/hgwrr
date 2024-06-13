cmake_minimum_required(VERSION 3.12.0)
project(hgwrr)

set(HGWRR_VERSION_PATCH_ADD 1)
math(EXPR HGWRR_VERSION_PATCH "${CMAKE_PROJECT_VERSION_PATCH}+${HGWRR_VERSION_PATCH_ADD}" OUTPUT_FORMAT DECIMAL)
set(PROJECT_VERSION_R "${CMAKE_PROJECT_VERSION_MAJOR}.${CMAKE_PROJECT_VERSION_MINOR}-${HGWRR_VERSION_PATCH}")

set(PROJECT_RBUILD_DIR ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
make_directory(${PROJECT_RBUILD_DIR})

find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIRS} ${RCPP_ARMADILLO_INCLUDE_DIR} ${RCPP_INCLUDE_DIR})

configure_file(DESCRIPTION.in ${PROJECT_RBUILD_DIR}/DESCRIPTION)
file(COPY cleanup DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY configure DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY configure.ac DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/R DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/src DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/man DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/data DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/tests DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/include/hlmgwr.h DESTINATION ${PROJECT_RBUILD_DIR}/src)
file(COPY ${CMAKE_SOURCE_DIR}/src/hlmgwr.cpp DESTINATION ${PROJECT_RBUILD_DIR}/src)

include_directories(../include)
add_definitions(-DRCPP_USE_GLOBAL_ROSTREAM)

add_custom_target(hgwrr_roxygen
    VERBATIM
    WORKING_DIRECTORY ${PROJECT_RBUILD_DIR}
    COMMAND_EXPAND_LISTS
    COMMAND ${R_EXECUTABLE} -e "roxygen2::roxygenize(load_code = 'source')"
)

add_custom_target(hgwrr_install_temp
    VERBATIM
    WORKING_DIRECTORY ${PROJECT_RBUILD_DIR}/..
    COMMAND_EXPAND_LISTS
    COMMAND ${CMAKE_COMMAND} -E
        make_directory "${PROJECT_NAME}.library"
    COMMAND ${R_EXECUTABLE} CMD INSTALL --preclean --library=${PROJECT_NAME}.library ${PROJECT_NAME}
)
add_dependencies(hgwrr_install_temp hgwrr_roxygen)

add_custom_target(hgwrr_rpkg
    VERBATIM
    WORKING_DIRECTORY ${PROJECT_RBUILD_DIR}/..
    COMMAND_EXPAND_LISTS
    COMMAND ${R_EXECUTABLE} CMD build ${PROJECT_NAME}
)
add_dependencies(hgwrr_rpkg hgwrr_roxygen)

add_library(hgwrr_rcpp_export SHARED src/RcppExports.cpp)
target_link_libraries(hgwrr_rcpp_export hgwrbml)
set_target_properties(hgwrr_rcpp_export PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_test(
    NAME R_Check_CRAN
    COMMAND ${R_EXECUTABLE} CMD check ${PROJECT_NAME}_${PROJECT_VERSION_R}.tar.gz --as-cran --no-manual
    WORKING_DIRECTORY ${PROJECT_RBUILD_DIR}/..
)

add_test(
    NAME R_TestThat
    COMMAND ${R_EXECUTABLE} -e "testthat::test_local()"
    WORKING_DIRECTORY ${PROJECT_RBUILD_DIR}
)
