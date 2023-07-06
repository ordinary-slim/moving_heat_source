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
  NAME 2D_triaBounDetection
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/2D_triaBounDetection"
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
add_test(
  NAME 2D_nd_coupling
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/2D_nd_coupling"
  )
add_test(
  NAME 3D_cubeLaplace
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/3D_cubeLaplace"
  )
add_test(
  NAME 3D_printing_5on5
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/3D_printing_5on5"
  )
