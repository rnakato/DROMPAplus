add_library(pw_func
  STATIC
pw_makefile.cpp GenomeCoverage.cpp GCnormalization.cpp ReadMpbldata.cpp pw_strShiftProfile.cpp
  )

target_include_directories(pw_func
	 PUBLIC ${PROJECT_SOURCE_DIR}/src
	 PUBLIC ${PROJECT_SOURCE_DIR}/src/parse2wig
	 PUBLIC ${PROJECT_SOURCE_DIR}/src/common
)
