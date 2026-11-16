# One-command build on any platform (Windows, macOS, Linux):
#     cmake -P build.cmake
# Configures and builds SNaQ++ in ./build and copies the program to ./bin.
cmake_minimum_required(VERSION 3.16)
set(src "${CMAKE_CURRENT_LIST_DIR}")
set(bld "${src}/build")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${src}" -B "${bld}" -DCMAKE_BUILD_TYPE=Release
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "CMake configuration failed")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${bld}" --config Release --parallel
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "Build failed")
endif()
file(MAKE_DIRECTORY "${src}/bin")
foreach(cand "${bld}/snaq" "${bld}/snaq.exe" "${bld}/Release/snaq.exe")
  if(EXISTS "${cand}")
    file(COPY "${cand}" DESTINATION "${src}/bin")
    get_filename_component(name "${cand}" NAME)
    message(STATUS "Done: the snaq program is at bin/${name}")
  endif()
endforeach()
