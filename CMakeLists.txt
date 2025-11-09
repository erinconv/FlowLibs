cmake_minimum_required (VERSION 3.9)

project(FlowLibs)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(GDAL)

add_subdirectory(submodules/Barnes2019-DepressionHierarchy EXCLUDE_FROM_ALL)
# add_subdirectory(submodules/Barnes2019-DepressionHierarchy/submodules/richdem EXCLUDE_FROM_ALL)

include_directories(submodules)

add_library(FlowLibs INTERFACE)
target_include_directories(FlowLibs INTERFACE include/)
target_compile_features(FlowLibs INTERFACE cxx_auto_type cxx_std_17)

if(USE_GDAL AND GDAL_FOUND)
  
  #Executable for main custom carving algorithm
  add_executable(main_carve.exe src/main_carve.cpp)
  target_link_libraries(main_carve.exe PRIVATE dephier richdem)
  target_include_directories(main_carve.exe PRIVATE src)
  target_compile_features(main_carve.exe PUBLIC cxx_auto_type cxx_std_17)
  
  # Standalone watershed delineation executable
  add_executable(watershed_delineation.exe src/main_watershed.cpp)
  target_link_libraries(watershed_delineation.exe PRIVATE dephier richdem)
  target_include_directories(watershed_delineation.exe PRIVATE src)
  target_compile_features(watershed_delineation.exe PUBLIC cxx_auto_type cxx_std_17)

  install(TARGETS main_carve.exe DESTINATION bin)
  install(TARGETS watershed_delineation.exe DESTINATION bin)
else()
  message(WARNING "GDAL is required to build the DepHier example programs; it is not needed to use the DepHier library in another project.")
endif()
