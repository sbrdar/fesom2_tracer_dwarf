cmake_minimum_required(VERSION 3.5)

project(fesom Fortran)

if(DEFINED ENV{FESOM_PLATFORM_STRATEGY})
   set(FESOM_PLATFORM_STRATEGY $ENV{FESOM_PLATFORM_STRATEGY} CACHE STRING "switch to platform specific compile settings, this is usually determined via the env.sh script")
else()
   set(FESOM_PLATFORM_STRATEGY "notset" CACHE STRING "switch to platform specific compile settings, this is usually determined via the env.sh script")
endif()

option(ENABLE_OPENACC "compile with OpenACC support" ON)
option(DISABLE_OPENACC_ATOMICS "disable kernels using atomic statement for reproducible results" OFF)
option(ENABLE_NVTX "enable use of nvtx tools" OFF)

option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)

# Precision options
option(BUILD_SINGLE_PRECISION "Build single precision version" OFF)
option(BUILD_DOUBLE_PRECISION "Build double precision version (default)" ON)
option(BUILD_BOTH_PRECISIONS "Build both single and double precision versions" OFF)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
#if(${USE_ICEPACK})
#   file(GLOB sources_Fortran ${src_home}/*.F90
#                             ${src_home}/icepack_drivers/*.F90
#                             ${src_home}/icepack_drivers/Icepack/columnphysics/*.F90)
#else()
file(GLOB sources_Fortran ${src_home}/*.F90)
#endif()
#list(REMOVE_ITEM sources_Fortran ${src_home}/fesom_partition_init.F90)
#file(GLOB sources_C ${src_home}/*.c)

# generate a custom file from fesom_version_info.F90 which includes the current git SHA
#set(FESOM_ORIGINAL_VERSION_FILE ${src_home}/fesom_version_info.F90)
#set(FESOM_GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/fesom_version_info-generated.F90)
#list(REMOVE_ITEM sources_Fortran ${FESOM_ORIGINAL_VERSION_FILE}) # we want to compile the generated file instead
#list(APPEND sources_Fortran ${FESOM_GENERATED_VERSION_FILE})
#add_custom_command(OUTPUT 5303B6F4_E4F4_45B2_A6E5_8E2B9FB5CDC4 ${FESOM_GENERATED_VERSION_FILE} # the first arg to OUTPUT is a name for a file we never create to make sure this command will run on every re-build (let our file be the second arg, as the first file is inadvertently removed by make)
#                   COMMAND ${CMAKE_COMMAND} -DFESOM_ORIGINAL_VERSION_FILE=${FESOM_ORIGINAL_VERSION_FILE} -DFESOM_GENERATED_VERSION_FILE=${FESOM_GENERATED_VERSION_FILE} -P GitRepositoryInfo.cmake
#                   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
#                   COMMENT "determining ${PROJECT_NAME} git SHA ...")

#if(${FESOM_STANDALONE})
#   list(REMOVE_ITEM sources_Fortran ${src_home}/cpl_driver.F90)
#endif()
#list(REMOVE_ITEM sources_Fortran ${src_home}/fvom_init.F90)
#list(REMOVE_ITEM sources_C       ${src_home}/fort_part.c)

# depends on the metis library
#add_subdirectory(../lib/metis-5.1.0 ${PROJECT_BINARY_DIR}/metis)
#include_directories(../lib/metis-5.1.0/include)

#add_subdirectory(async_threads_cpp)

#include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindNETCDF.cmake)

#add_library(${PROJECT_NAME}_C ${sources_C})
#target_compile_definitions(${PROJECT_NAME}_C PRIVATE USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)

# Determine which precision versions to build
if(BUILD_BOTH_PRECISIONS)
   set(BUILD_SINGLE_PRECISION ON)
   set(BUILD_DOUBLE_PRECISION ON)
endif()

# Function to create executable with specific precision
function(create_precision_executable PRECISION_NAME PRECISION_FLAG)
   set(EXEC_NAME ${PROJECT_NAME}_${PRECISION_NAME})
   add_executable(${EXEC_NAME} ${sources_Fortran})
   
   # Add precision-specific preprocessor definition
   if(PRECISION_FLAG)
      target_compile_definitions(${EXEC_NAME} PRIVATE USE_SINGLE_PRECISION)
   endif()
   
   if(${VERBOSE})
      target_compile_definitions(${EXEC_NAME} PRIVATE VERBOSE)
   endif()
   
   # Store executable name for later configuration
   set(${PRECISION_NAME}_EXECUTABLE ${EXEC_NAME} PARENT_SCOPE)
endfunction()

# Create executables based on options
if(BUILD_DOUBLE_PRECISION)
   create_precision_executable("dp" OFF)
   message(STATUS "Building double precision version: ${PROJECT_NAME}_dp")
endif()

if(BUILD_SINGLE_PRECISION)
   create_precision_executable("sp" ON)
   message(STATUS "Building single precision version: ${PROJECT_NAME}_sp")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# Function to apply compiler flags to an executable
function(apply_compiler_flags TARGET_NAME PRECISION_NAME)
   # CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
   if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
      # For Intel, -r8 forces double precision, -r4 for single precision
      if(PRECISION_NAME STREQUAL "sp")
         target_compile_options(${TARGET_NAME} PRIVATE 
                                $<$<CONFIG:DEBUG>:-r4 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit>
                                $<$<CONFIG:RELEASE>:-r4 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -ip -g -traceback>
                               )
      else()
         target_compile_options(${TARGET_NAME} PRIVATE 
                                $<$<CONFIG:DEBUG>:-r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit>
                                $<$<CONFIG:RELEASE>:-r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -ip -g -traceback>
                               )
      endif()
      if(${FESOM_PLATFORM_STRATEGY} STREQUAL levante.dkrz.de)
         target_compile_options(${TARGET_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
      else()
         target_compile_options(${TARGET_NAME} PRIVATE -xHost)
      endif()
   elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
      # For GNU, -fdefault-real-8 forces double precision
      # Single precision relies on WP=4 in source code (no compiler flag needed)
      if(PRECISION_NAME STREQUAL "sp")
         target_compile_options(${TARGET_NAME} PRIVATE -O3 -finit-local-zero -finline-functions -march=native -fimplicit-none -ffree-line-length-none)
      else()
         target_compile_options(${TARGET_NAME} PRIVATE -O3 -finit-local-zero -finline-functions -march=native -fimplicit-none -fdefault-real-8 -ffree-line-length-none)
      endif()
      if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
         target_compile_options(${TARGET_NAME} PRIVATE -fallow-argument-mismatch)
      endif()
   elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray)
      # For Cray, -s real64 forces double precision, -s real32 for single
      if(PRECISION_NAME STREQUAL "sp")
         target_compile_options(${TARGET_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real32 -hnoomp)
      else()
         target_compile_options(${TARGET_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -hnoomp)
      endif()
   elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL NVHPC)
      # For PGI/NVHPC, -Mr8 forces double precision, -Mr4 for single
      if(${ENABLE_OPENACC})
         if(PRECISION_NAME STREQUAL "sp")
            target_compile_options(${TARGET_NAME} PRIVATE
                                   $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr4 -pgf90libs -Mnofma -Minfo=all -acc=verystrict -gpu=math_uniform,cuda11.7,cc80>
                                   $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr4 -pgf90libs -Minfo=all -acc=verystrict -gpu=cc80>
                                  )
         else()
            target_compile_options(${TARGET_NAME} PRIVATE
                                   $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr8 -pgf90libs -Mnofma -Minfo=all -acc=verystrict -gpu=math_uniform,cuda11.7,cc80>
                                   $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr8 -pgf90libs -Minfo=all -acc=verystrict -gpu=cc80>
                                  )
         endif()
         set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-acc=verystrict -Mnofma -gpu=math_uniform,cuda11.7,cc80")
         set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-acc=verystrict -gpu=cc80")
         if(${DISABLE_OPENACC_ATOMICS})
            target_compile_definitions(${TARGET_NAME} PRIVATE DISABLE_OPENACC_ATOMICS)
         endif()
      else()
         if(PRECISION_NAME STREQUAL "sp")
            target_compile_options(${TARGET_NAME} PRIVATE 
                                   $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr4 -pgf90libs -Mnofma>
                                   $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr4 -pgf90libs>
                                  )
         else()
            target_compile_options(${TARGET_NAME} PRIVATE 
                                   $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr8 -pgf90libs -Mnofma>
                                   $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr8 -pgf90libs>
                                  )
         endif()
      endif()
      if(${ENABLE_NVTX})
         target_link_libraries(${TARGET_NAME} nvhpcwrapnvtx)
      endif()
   endif()
endfunction()

# Apply compiler flags to each built executable
if(BUILD_DOUBLE_PRECISION)
   apply_compiler_flags(${PROJECT_NAME}_dp "dp")
endif()

if(BUILD_SINGLE_PRECISION)
   apply_compiler_flags(${PROJECT_NAME}_sp "sp")
endif()


#elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
#       target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -hnoomp)
#endif()
#target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
#target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} async_threads_cpp)
#set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

# Install executables
set(FESOM_INSTALL_PATH "${CMAKE_CURRENT_LIST_DIR}/../bin" CACHE PATH "directory where the FESOM binaries should be installed")

if(BUILD_DOUBLE_PRECISION)
   install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_dp DESTINATION ${FESOM_INSTALL_PATH} RENAME fesom_tracer_mesh_init_dp)
endif()

if(BUILD_SINGLE_PRECISION)
   install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_sp DESTINATION ${FESOM_INSTALL_PATH} RENAME fesom_tracer_mesh_init_sp)
endif()
