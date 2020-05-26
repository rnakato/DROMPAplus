add_library(dd_func
  STATIC
  dd_init.cpp dd_draw_dataframe.cpp dd_classfunc_draw.cpp dd_command.cpp dd_readfile.cpp dd_draw.cpp dd_chiadrop.cpp dd_drawgenes.cpp dd_sample_definition.cpp dd_profile.cpp
  )

find_package(PkgConfig REQUIRED)

pkg_check_modules(GTKMM gtkmm-3.0)
if(NOT GTKMM_FOUND)
    pkg_check_modules(GTKMM gtkmm-2.4 REQUIRED)
endif()

include_directories(${GTKMM_INCLUDE_DIRS})

target_include_directories(dd_func
	 PUBLIC ${PROJECT_SOURCE_DIR}/src
	 PUBLIC ${PROJECT_SOURCE_DIR}/src/drompa
	 PUBLIC ${PROJECT_SOURCE_DIR}/src/common
)
