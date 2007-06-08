#include "partition.h"
#include <iostream>
#include <fstream>

using namespace dolfin;

void metis(Mesh& mesh, int num_part, int max_cells, char* filename)
{
  MeshFunction<dolfin::uint> parts;
  int num_cells = mesh.numCells();
  std::ofstream resultfile(filename);
  resultfile << "Metis with " << num_part << " partitions" << std::endl;
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
      metisPartitioning(graph, num_part, reinterpret_cast<int*>(parts.values()));
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
void scotch(Mesh& mesh, int num_part, int max_cells, char* filename)
{
  MeshFunction<dolfin::uint> parts;
  int num_cells = mesh.numCells();
  std::ofstream resultfile(filename);
  resultfile << "Scotch with " << num_part << " partitions" << std::endl;
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
      scotchPartitioning(graph, num_part, reinterpret_cast<int*>(parts.values()));
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
      " <metis|scotch> <meshfile> <num_part> <max_cells> <resultfile>\n";
    exit(1);
  }

  char* partitioner = argv[1];
  char* infilename = argv[2];
  int num_part = atoi(argv[3]);
  int max_cells = atoi(argv[4]);
  char* outfilename = argv[5];

  Mesh mesh(infilename);
  if(strcmp(partitioner, "metis") == 0)
  {
    std::cout << "Using Metis to partition mesh into " << num_part << " partitions " << std::endl;
    metis(mesh, num_part, max_cells, outfilename);
    std::cout << "\n\n";
  }
  else if(strcmp(partitioner, "scotch") == 0)
  {
    std::cout << "Using Scotch to partition mesh into " << num_part << " partitions " << std::endl;
    scotch(mesh, num_part, max_cells, outfilename);
    std::cout << "\n\n";
  }
  else
  {
    std::cerr << "Usage: " << argv[0] << 
      " <metis|scotch> <meshfile> <num_part> <max_cells> <resultfile>\n";
    exit(1);
  }
}
