set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

if(APPLE)
  set(BOOST_LIBS "-lpthread -lboost_thread-mt -lboost_program_options -lboost_system -lboost_filesystem")
else()
  set(BOOST_LIBS "-lpthread -lboost_thread -lboost_program_options -lboost_system -lboost_filesystem")
endif()

link_directories(parse2wig+ "/usr/local/lib")
link_directories(drompa+ "/usr/local/lib")

#find_package(Boost 1.53.0 REQUIRED thread program_options system filesystem)
#include_directories(${Boost_INCLUDE_DIRS})

add_executable(parse2wig+ pw_main.cpp)

target_link_libraries(parse2wig+
                      pw_func
                      common
                      ssp_func
                      ssp_common
                      ${PROJECT_SOURCE_DIR}/submodules/SSP/src/htslib-1.10.2/libhts.a
                      ${BOOST_LIBS}
#                      ${Boost_LIBRARIES}
                      -lcurl -llzma -lbz2 -lz
                      -lgsl -lgslcblas
                      ${GTKMM_LIBRARIES}
)

add_executable(drompa+ dd_main.cpp)

target_link_libraries(drompa+
                      dd_func
                      common
                      ssp_func
                      ssp_common
                      ${PROJECT_SOURCE_DIR}/submodules/SSP/src/htslib-1.10.2/libhts.a
                      ${BOOST_LIBS}
#                      ${Boost_LIBRARIES}
                      -lcurl -llzma -lbz2 -lz
                      -lgsl -lgslcblas
                      ${GTKMM_LIBRARIES}
)
