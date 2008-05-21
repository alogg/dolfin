#include "assembly_tester.h"

using namespace dolfin;

void dolfin::assembly_tester(Mesh& mesh, Form& a, std::string name, int n_reasm)
{
  std::string* backends = new std::string[4];
  backends[0] = "uBLAS";
  backends[1] = "PETSc";
  backends[2] = "Epetra";
  backends[3] = "Assembly";
  for (int j = 0; j < 4; ++j){
    ufc::finite_element* finite_element = 
      a.form().create_finite_element(0);
    //printf("Element: %s\n",finite_element->signature() );
    delete finite_element;
    
    dolfin_set("linear algebra backend", backends[j].c_str());
    Matrix A;
    Assembler ass(mesh);
    std::string backend = dolfin_get("linear algebra backend");
    printf("Test: %s with backend: %s\n",name.c_str(), backend.c_str());
    
    tic();
    ass.assemble(A,a);
    real t = toc();
    
    Array<uint> cols;
    Array<real> vals; 
    uint nnz = 0;
    for (uint i = 0; i < A.size(0); ++i){
      A.getrow(i, cols, vals);
      nnz += cols.size();       
    }
    
    //  printf("Size: %d x %d\n",A.size(0), A.size(1));
    printf("Nonzeros: %d\n",nnz);
    
    printf("\tFull assembly: %f seconds\n",t);
    printf("\tRate: %g (non-zero matrix elements per second)\n\n",nnz/t);
    
    /*for(int i = 0; i<n_reasm; i++)
      {
      tic();
      ass.assemble(A,a,false);
      t = toc();
      printf("\tReassembly: %f seconds\t rate: %g\n",t, nnz/t);
      }
    */
    //  printf("\n");
  }
}
