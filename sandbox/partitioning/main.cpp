#include "partition.h"
#include <iostream>
#include <fstream>

using namespace dolfin;

void partition(Mesh& mesh, int num_part, int max_cells, char* filename,
               void (partFunc(Graph&, int, int*)), char* plotname)
{
  MeshFunction<dolfin::uint> parts;
  int num_cells = mesh.numCells();
  std::ofstream resultfile(filename);
  resultfile << plotname << " with " << num_part << " partitions" << std::endl;
  resultfile << "Cells " << "Seconds " << "Edgecut " << "Balance" << std::endl;
  while(num_cells < max_cells)
  {
    Graph graph(mesh);
    parts.init(mesh, mesh.topology().dim());

    real time = 0;
    int runs = 0;

    std::cout << "\nPartitioning mesh with " << num_cells << " cells" << std::endl;
    do
    {
      tic();
      partFunc(graph, num_part, reinterpret_cast<int*>(parts.values()));
      time += toc();
      runs++;
    }
    while(time < 1);
    real runtime = time/runs;
    int edgecut = GraphPartition::edgecut(graph, num_part, parts.values());

    resultfile << num_cells << " " << runtime << " " << edgecut << std::endl;
    std::cout << "Time to partition: " << runtime << std::endl;
    std::cout << "Edge cut: " << edgecut << std::endl;

    std::cout << "Refining mesh... " << std::endl;
    mesh.refine();
    num_cells = mesh.numCells();
  }
  resultfile.close();

}
int main(int argc, char* argv[])
{
  if(argc != 6)
  {
    std::cerr << "Usage: " << argv[0] << 
      " <metis1|metis2|scotch> <meshfile> <num_part> <max_cells> <resultfile>\n";
    exit(1);
  }

  char* partitioner = argv[1];
  char* infilename = argv[2];
  int num_part = atoi(argv[3]);
  int max_cells = atoi(argv[4]);
  char* outfilename = argv[5];

  Mesh mesh(infilename);
  if(strcmp(partitioner, "metis1") == 0)
  {
    std::cout << "Using Metis Kway to partition mesh into " << num_part << " partitions " << std::endl;
    partition(mesh, num_part, max_cells, outfilename, 
              metisKwayPartitioning, "Metis Kway");
    std::cout << "\n\n";
  }
  else if(strcmp(partitioner, "metis2") == 0)
  {
    std::cout << "Using Metis recursive to partition mesh into " << num_part << " partitions " << std::endl;
    partition(mesh, num_part, max_cells, outfilename, 
              metisRecursivePartitioning, "Metis recursive");
    std::cout << "\n\n";
  }
  else if(strcmp(partitioner, "scotch") == 0)
  {
    std::cout << "Using Scotch to partition mesh into " << num_part << " partitions " << std::endl;
    partition(mesh, num_part, max_cells, outfilename,
                 scotchPartitioning, "Scotch");
    std::cout << "\n\n";
  }
  else
  {
    std::cerr << "Usage: " << argv[0] << 
      " <metis1|metis2|scotch> <meshfile> <num_part> <max_cells> <resultfile>\n";
    exit(1);
  }
}
