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
add_test(
  NAME boundaryDetection
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/boundaryDetection"
  )
add_test(
  NAME 2D_checkConvection
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/2D_checkConvection"
  )
add_test(
  NAME 1D_checkNeumann
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/1D_checkNeumann"
  )
add_test(
  NAME 2D_checkNeumannWactivation
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/2D_checkNeumannWactivation"
  )
add_test(
  NAME 1D_checkAdvection
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/1D_checkAdvection"
  )
