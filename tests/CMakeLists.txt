add_test(
  NAME 1D_checkDirichlet
  COMMAND ${PYTHON_EXECUTABLE} -m pytest test.py
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/1D_checkDirichlet"
  )
