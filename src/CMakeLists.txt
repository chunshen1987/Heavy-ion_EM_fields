add_executable (EM_fields.e
  main.cpp
  EM_fields.cpp
  gauss_quadrature.cpp
  )
target_link_libraries (EM_fields.e ${LIBS})

install(TARGETS EM_fields.e DESTINATION ${CMAKE_HOME_DIRECTORY})
