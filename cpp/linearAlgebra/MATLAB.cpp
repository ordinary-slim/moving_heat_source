#ifdef hasMATLAB
#include "MATLAB.h"
#include "Solver.h"
#include "LinearSystem.h"

namespace mdat = matlab::data;

void MatLabSolver::solve(LinearSystem *ls) {
  /*
   * Copy lhs and rhs to MatLab DSs
   * Solve
   * Copy solution back to Eigen DS
   */
  // Create buffers for the data
  std::cout << "Using MATLAB linear system solver" << std::endl;

  unsigned long nnz = ls->lhs.nonZeros();
  mdat::ArrayFactory factory;
  mdat::buffer_ptr_t<double> data_sparse_p = factory.createBuffer<double>(nnz);
  mdat::buffer_ptr_t<size_t> rows_p = factory.createBuffer<size_t>(nnz);
  mdat::buffer_ptr_t<size_t> cols_p = factory.createBuffer<size_t>(nnz);

  // Write data into the buffers
  double* dataSparsePtr = data_sparse_p.get();
  size_t* rowsPtr = rows_p.get();
  size_t* colsPtr = cols_p.get();
  int idx = 0;
  for (int k = 0; k < ls->lhs.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(ls->lhs, k); it; ++it) {
      rowsPtr[idx] = it.row();
      colsPtr[idx] = it.col();
      dataSparsePtr[idx] = it.value();
      ++idx;
    }
  }
  // Use the buffers to create the sparse array
  const mdat::SparseArray<double> lhsMatlab =
      factory.createSparseArray<double>({static_cast<unsigned long>(ls->lhs.rows()),
                                         static_cast<unsigned long>(ls->lhs.cols())},
                                         nnz, 
                                         std::move(data_sparse_p),
                                         std::move(rows_p),
                                         std::move(cols_p));
  const mdat::TypedArray<double> rhsMatlab =
      factory.createArray({static_cast<unsigned long>(ls->rhs.size())},
                           ls->rhs.begin(),
                           ls->rhs.end()
          );
  mdat::TypedArray<double> solMatlab = matlab->ptr->feval("mldivide",
                                                          {lhsMatlab, rhsMatlab});
  for (int i = 0; i < ls->sol.size(); ++i) {
    ls->sol[i] = solMatlab[i];
  }
}
void MatLabSolver::setSession(MatLabSession* ms) {
  matlab = ms;
}
#endif
