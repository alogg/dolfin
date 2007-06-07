#include <dolfin.h>
#include <iostream>

#include <parmetis.h>
extern "C"
{
  #include <metis.h>
  #include <scotch.h>
}
using namespace dolfin;

void metisPartitioning(Graph& graph, int num_part, idxtype* vtx_part)
{
  int options[5];
  int num_vertices = graph.numVertices();
  int wgtflag = 0;
  int pnumflag = 0;
  int edgecut = 0;
  
  // This results in segfault
  //idxtype* adjncy = reinterpret_cast<idxtype*>(graph.connectivity());
  //idxtype* xadj = reinterpret_cast<idxtype*>(graph.offsets() + 1);

  // Copying values instead
  idxtype* adjncy = new idxtype[graph.numArches()];
  idxtype* xadj = new idxtype[graph.numVertices() + 1];
  for(unsigned int i=0; i<graph.numArches(); ++i)
    adjncy[i] = graph.connectivity()[i];
  for(unsigned int i=0; i<=graph.numVertices(); ++i)
    xadj[i] = graph.offsets()[i];

  options[0] = 0;
  METIS_PartGraphKway(&num_vertices, xadj, adjncy, NULL, NULL, &wgtflag, &pnumflag, &num_part, options, &edgecut, vtx_part);
}

void scotchPartitioning(Graph& graph, int num_part, int* vtx_part)
{
  SCOTCH_Graph grafdat;
  SCOTCH_Strat strat;

  if (SCOTCH_graphInit (&grafdat) != 0) {
  }
  if (SCOTCH_graphBuild (&grafdat, 0, static_cast<int>(graph.numVertices()), reinterpret_cast<int*>(graph.offsets()), NULL, NULL, NULL, static_cast<int>(graph.numArches()), reinterpret_cast<int*>(graph.connectivity()), NULL) != 0) {
  }

  SCOTCH_stratInit(&strat);

  /*
  if (SCOTCH_stratGraphBipart (&strat, "b{job=t,map=t,poli=S,strat=m{type=h,vert=80,low=h{pass=10}f{bal=0.0005,move=80},asc=b{bnd=d{dif=1,rem=1,pass=40}f{bal=0.005,move=80},org=f{bal=0.005,move=80}}}|m{type=h,vert=80,low=h{pass=10}f{bal=0.0005,move=80},asc=b{bnd=d{dif=1,rem=1,pass=40}f{bal=0.005,move=80},org=f{bal=0.005,move=80}}}}")) {
  }
  */

  if (SCOTCH_graphPart (&grafdat, num_part, &strat, vtx_part) != 0) {
  }

  SCOTCH_stratExit (&strat);
  SCOTCH_graphExit (&grafdat);
}
