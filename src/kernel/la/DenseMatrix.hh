#ifndef __DENSE_MATRIX_HH
#define __DENSE_MATRIX_HH

#include <kw_constants.h>

class DirectSolver;
class SparseMatrix;

class DenseMatrix{
public:

  DenseMatrix(int m, int n);
  DenseMatrix(SparseMatrix *A);
  ~DenseMatrix();

  void Resize(int m, int n);
  
  void     Set        (int i, int j, real value);
  int      Size       (int dim);
  real     Get        (int i, int j);

  void Display();
  void DisplayAll();
  void DisplayRow(int i);
  void Write(const char *filename);
  
  friend class DirectSolver;
  
private:
  
  int m,n;
  real **values;
  int *permutation;
  
};

#endif
