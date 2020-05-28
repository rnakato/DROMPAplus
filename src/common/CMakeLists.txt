add_library(common
  STATIC
  util.cpp WigStats.cpp significancetest.cpp statistics.cpp
  )

target_include_directories(common
	 PUBLIC ${PROJECT_SOURCE_DIR}/src
	 PUBLIC ${PROJECT_SOURCE_DIR}/src/common
)
