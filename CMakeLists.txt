cmake_minimum_required(VERSION 3.12.0)
project(hgwrr VERSION 0.1.0)

set(PROJECT_VERSION_R "0.1-0")

set(PROJECT_RBUILD_DIR ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
make_directory(${PROJECT_RBUILD_DIR})

find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIRS} ${RCPP_ARMADILLO_INCLUDE_DIR} ${RCPP_INCLUDE_DIR})

configure_file(DESCRIPTION.in ${PROJECT_RBUILD_DIR}/DESCRIPTION)
configure_file(NAMESPACE.in ${PROJECT_RBUILD_DIR}/NAMESPACE)
file(COPY ${CMAKE_SOURCE_DIR}/R/R DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/src DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/man DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/R/data DESTINATION ${PROJECT_RBUILD_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/include/hlmgwr.h DESTINATION ${PROJECT_RBUILD_DIR}/src)
file(COPY ${CMAKE_SOURCE_DIR}/src/hlmgwr.cpp DESTINATION ${PROJECT_RBUILD_DIR}/src)

include_directories(../include)

add_custom_target(hgwrr_rbuild
    VERBATIM
    WORKING_DIRECTORY ${PROJECT_RBUILD_DIR}/..
    COMMAND_EXPAND_LISTS
    COMMAND ${CMAKE_COMMAND} -E
        make_directory "${PROJECT_NAME}.library"
    COMMAND ${R_EXECUTABLE} CMD INSTALL --preclean --clean --library="${PROJECT_NAME}.library" ${PROJECT_NAME}
    COMMAND ${R_EXECUTABLE} -e "devtools::document('${PROJECT_NAME}')"
    COMMAND ${R_EXECUTABLE} CMD build ${PROJECT_NAME}
)

add_library(hgwrr_rcpp_export SHARED src/RcppExports.cpp)
target_link_libraries(hgwrr_rcpp_export hgwrbml)
set_target_properties(hgwrr_rcpp_export PROPERTIES EXCLUDE_FROM_ALL TRUE)

set(HGWRR_TEST_SCRIPT_DIR ${hgwrr_BINARY_DIR}/bin)
file(COPY test/test.R DESTINATION ${HGWRR_TEST_SCRIPT_DIR})

add_test(
    NAME Test_R_hgwr
    COMMAND ${R_EXECUTABLE} < test.R TEST_DATA_DIR=${TEST_DATA_DIR}
    WORKING_DIRECTORY ${HGWRR_TEST_SCRIPT_DIR}
)
