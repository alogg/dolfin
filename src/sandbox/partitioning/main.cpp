#include <dolfin.h>
//#include <dolfin/GraphPartition.h>
#include <iostream>

/*
#include <parmetis.h>
#include <scotch.h>
extern "C"
{
  #include <metis.h>
}
*/
using namespace dolfin;

/*
void testMetisGraph(Graph& graph, int num_partitions)
{
  int options[5];
  int num_vertices = graph.numVertices();
  int wgtflag = 0;
  int pnumflag = 0;
  int edgecut = 0;
  idxtype* parts = new int[num_vertices];
  idxtype* adjncy = (idxtype*) graph.connectivity();
  idxtype* xadj = (idxtype*) graph.offsets();
  options[0] = 0;
  std::cout << "METIS_PartGraphKway" << std::endl;
  METIS_PartGraphRecursive(&num_vertices, xadj, adjncy, NULL, NULL, &wgtflag, &pnumflag, &num_partitions, options, &edgecut, parts);
}
//-----------------------------------------------------------------------------
void testMetisMesh(Mesh& mesh, int num_partitions)
{
  int num_cells     = mesh.numCells() ;
  int num_vertices  = mesh.numVertices();
  
  int index_base = 0;  // zero-based indexing
  int edges_cut  = 0;

  int cell_type = 0;
  idxtype* cell_partition   = new int[num_cells];
  idxtype* vertex_partition = new int[num_vertices];
  idxtype* mesh_data = 0;

  // Set cell type and allocate memory for METIS mesh structure
  if(mesh.type().cellType() == CellType::triangle)
  {
    cell_type = 1;
    mesh_data = new int[3*num_cells];
  }
  else if(mesh.type().cellType() == CellType::tetrahedron) 
  {
    cell_type = 2;
    mesh_data = new int[4*num_cells];
  }
  else
    dolfin_error("Do not know how to partition mesh of this type");
  
  if(num_partitions > 1)
  {
    // Create mesh structure for METIS
    dolfin::uint i = 0;
    for (CellIterator cell(mesh); !cell.end(); ++cell)
      for (VertexIterator vertex(cell); !vertex.end(); ++vertex)
        mesh_data[i++] = vertex->index();

      // Use METIS to partition mesh
    METIS_PartMeshNodal(&num_cells, &num_vertices, mesh_data, &cell_type, &index_base, 
                        &num_partitions, &edges_cut, cell_partition, vertex_partition);
  
  }
  // Clean up
  delete [] cell_partition;
  delete [] vertex_partition;
  delete [] mesh_data;
}
//-----------------------------------------------------------------------------

void testScotch(Graph& graph, int num_partitions)
{
  //SCOTCH_Graph grafdat;
  //FILE* fileptr;

  //if ((fileptr = fopen ("/home/magnus/doc/uio/master/dolfin_tests/fem_graphs/unitcube.grf", "r")) == NULL) {
  //}
  //if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) {
 // }
}
*/

void testDolfinGraph(Graph& graph, int num_part)
{
  std::cout << "Testing graph partitioning" << std::endl;
  dolfin::uint* parts = new dolfin::uint[graph.numVertices()];
  GraphPartition::partition(graph, num_part, parts);

  GraphPartition::check(graph, num_part, parts);
  GraphPartition::eval(graph, num_part, parts);
  GraphPartition::disp(graph, num_part, parts);
}

void testDolfinMesh(Mesh& mesh, int num_part)
{
  std::cout << "Testing mesh partitioning" << std::endl;
  MeshFunction<dolfin::uint> partitions;
  partitions.init(mesh, mesh.topology().dim());
  mesh.partition(num_part, partitions);

  //Graph graph(mesh, Graph::dual);
  //GraphPartition::check(graph, num_part, partitions.values());
  //GraphPartition::eval(graph, num_part, partitions.values());
  //GraphPartition::disp(graph, num_part, partitions.values());
  //Array<Mesh> meshes 
  //MeshPartitioning::partition(mesh, meshes, num_partitions);

  File file("mesh_partition.xml");
  file << partitions;
}

void createSimpleMesh(Mesh& mesh)
{
  /* 2---3
	* |\1 |\
	* | \ | \
	* |0 \|2 \
	* 0---1---4
	*
	*/
  MeshEditor editor;
  editor.open(mesh, "triangle", 2, 2);
  editor.initVertices(5);
  editor.addVertex(0, 0.4, 0.4);
  editor.addVertex(1, 0.7, 0.4);
  editor.addVertex(2, 0.4, 0.7);
  editor.addVertex(3, 0.7, 0.7);
  editor.addVertex(4, 1.0, 0.4);
  editor.initCells(3);
  editor.addCell(0, 0, 1, 2);
  editor.addCell(1, 2, 1, 3);
  editor.addCell(2, 3, 1, 4);
  editor.close();
}

int main(int argc, char* argv[])
{
  //testMetis(graph, 16);
  //testMeshPartition(mesh, 16);

  UnitSquare mesh(4, 4);

  Graph graph(mesh, "nodal");
  //testDolfinGraph(graph, 3);
  testDolfinMesh(mesh, 3);
}
