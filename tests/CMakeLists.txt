add_test(
  NAME 1D_checkDirichlet
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/1D_checkDirichlet"
  )
add_test(
  NAME 2D_rectangleTransportAdvection
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/2D_rectangleTransportAdvection"
  )
