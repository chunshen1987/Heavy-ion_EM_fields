add_executable (EM_fields.e
  main.cpp
  EM_fields.cpp
  ParameterReader.cpp
  )
target_link_libraries (EM_fields.e ${LIBS})

install(TARGETS EM_fields.e DESTINATION ${CMAKE_HOME_DIRECTORY})
