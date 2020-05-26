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

add_subdirectory(parse2wig)
add_subdirectory(drompa)
