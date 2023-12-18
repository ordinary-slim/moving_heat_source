#ifndef MYSOLVER
#define MYSOLVER

class LinearSystem;

class Solver {
  public:
    virtual void solve(LinearSystem* ls) {
    }
};

class EigenBiCGSTAB_IncompleteLUT : public Solver {
  void solve(LinearSystem *ls);
};

class EigenBiCGSTAB_Jacobi : public Solver {
  void solve(LinearSystem *ls);
};

class EigenCG : public Solver {
  void solve(LinearSystem *ls);
};
#endif
