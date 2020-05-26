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
