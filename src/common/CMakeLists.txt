add_library(common
  STATIC
  ReadAnnotation.cpp color.cpp util.cpp WigStats.cpp significancetest.cpp statistics.cpp
  )

target_include_directories(common
	 PUBLIC ${PROJECT_SOURCE_DIR}/src
	 PUBLIC ${PROJECT_SOURCE_DIR}/src/common
)
