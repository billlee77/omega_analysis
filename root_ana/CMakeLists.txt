cmake_minimum_required(VERSION 2.8)
project(DEMPgen)
set (DEMPgen_VERSION_MAJOR 0)
set (DEMPgen_VERSION_MINOR 1)

add_definitions (-DROOT_MATH_LOG)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#include_directories(include)
#include_directories("src/eic_evgen")
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
#include_directories(include/json)

#file(GLOB SOURCES "src/*" "src/eic_evgen/eic.cc" "src/eic_evgen/eic_pim.cc" "src/eic_evgen/tssa_sig_Para.cc" "src/eic_evgen/reaction_routine.cc" "src/eic_evgen/legacy_function.cc" "src/eic_evgen/legacy_function.cc" "src/eic_evgen/process_routine/*")
# file(GLOB SOURCES "src/*.cpp")
# file(GLOB SOURCES "src/*.cpp")

file(GLOB HEEP_SOURCES "root_ana_heep.cpp",  "analysis.cpp", "read_setting.cpp", "general_utility.cpp", "cut.cpp", "print_data.cpp", "analysis_heep.cpp")

file(GLOB HEEP_SINGLE_SOURCES "root_ana_heep_single.cpp",  "analysis_heep.cpp",  "analysis.cpp", "read_setting.cpp", "general_utility.cpp", "cut.cpp", "print_data.cpp", "analysis_heep_single.cpp")

file(GLOB PIONPLUS_SOURCES "root_ana_pl.cpp",  "analysis.cpp", "read_setting.cpp", "general_utility.cpp", "cut.cpp", "print_data.cpp" )

file(GLOB OMEGA_SOURCES "root_ana_omega.cpp",  "analysis.cpp", "analysis_heep.cpp", "read_setting.cpp", "general_utility.cpp", "cut.cpp", "print_data.cpp", "analysis_omega.cpp")

file(GLOB OMEGA_ITT_SOURCES "root_ana_omega_itt.cpp",  "analysis.cpp", "analysis_heep.cpp", "read_setting.cpp", "general_utility.cpp", "cut.cpp", "print_data.cpp", "analysis_omega.cpp" )

message("Root Include Dirs: " ${ROOT_INCLUDE_DIRS})
message("Root Libraries: " ${ROOT_LIBRARIES})
message("Mathmore: " ${ROOT_mathmore_FOUND})

add_executable (root_ana_heep ${HEEP_SOURCES})
target_link_libraries(root_ana_heep ${ROOT_LIBRARIES})

add_executable (root_ana_heep_single ${HEEP_SINGLE_SOURCES})
target_link_libraries(root_ana_heep_single ${ROOT_LIBRARIES})

add_executable (root_ana_pl ${PIONPLUS_SOURCES})
target_link_libraries(root_ana_pl ${ROOT_LIBRARIES})

add_executable (root_ana_omega ${OMEGA_SOURCES})
target_link_libraries(root_ana_omega ${ROOT_LIBRARIES})

add_executable (root_ana_omega_itt ${OMEGA_ITT_SOURCES})
target_link_libraries(root_ana_omega_itt ${ROOT_LIBRARIES})






