cmake_minimum_required(VERSION 3.12.0)
project(hgwrr VERSION 0.1.0)

set(PROJECT_VERSION_R "0.1-0")

set(PROJECT_RBUILD_DIR ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
make_directory(${PROJECT_RBUILD_DIR})

configure_file(DESCRIPTION.in ${PROJECT_RBUILD_DIR}/DESCRIPTION)
configure_file(NAMESPACE.in ${PROJECT_RBUILD_DIR}/NAMESPACE)

include_directories(../include)

add_custom_target(hgwrr_rbuild
    VERBATIM
    COMMAND_EXPAND_LISTS
    COMMAND ${CMAKE_COMMAND} -E
        make_directory
        ${PROJECT_RBUILD_DIR}
    COMMAND ${CMAKE_COMMAND} -E
        copy_directory ${CMAKE_SOURCE_DIR}/R/R ${PROJECT_RBUILD_DIR}/R
    COMMAND ${CMAKE_COMMAND} -E
        copy_directory ${CMAKE_SOURCE_DIR}/R/src ${PROJECT_RBUILD_DIR}/src
    COMMAND ${CMAKE_COMMAND} -E
        copy_directory ${CMAKE_SOURCE_DIR}/include ${PROJECT_RBUILD_DIR}/src/libhgwr
    COMMAND ${CMAKE_COMMAND} -E
        copy_directory ${CMAKE_SOURCE_DIR}/src ${PROJECT_RBUILD_DIR}/src/libhgwr
    WORKING_DIRECTORY ${PROJECT_RBUILD_DIR}
)
