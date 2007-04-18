#include <dolfin.h>
#include <iostream>

using namespace dolfin;

void testGraphEditor()
{
  /// Testing GraphEditor
  std::cout << "Testing GraphEditor" << std::endl;
  GraphEditor ge;
  Graph graph;
  ge.open(graph, Graph::undirected);
  ge.initVertices(4);
  ge.addVertex(0, 2);
  ge.addVertex(1, 3);
  ge.addVertex(2, 2);
  ge.addVertex(3, 3);

  ge.initEdges(5);
  ge.addEdge(0, 1);
  ge.addEdge(1, 2);
  ge.addEdge(2, 3);
  ge.addEdge(1, 3);
  ge.addEdge(0, 3);

  ge.close();

  graph.disp();
}
void testMeshToGraph()
{
  std::cout << "Testing Mesh to graph convertion" << std::endl;
  UnitSquare mesh(2, 2);

  Graph graph(mesh, "nodal");

}
void testGraphPartition()
{
  std::cout << "Testing graph partitioning" << std::endl;
  // Test graph
  /*
  
  Graph graph;
  GraphEditor ge;
  ge.open(graph, Graph::undirected);
  ge.initVertices(6);
  ge.addVertex(0, 3);
  ge.addVertex(1, 2);
  ge.addVertex(2, 3);
  ge.addVertex(3, 3);
  ge.addVertex(4, 2);
  ge.addVertex(5, 3);

  ge.initEdges(8);
  ge.addEdge(0, 1);
  ge.addEdge(0, 2);
  ge.addEdge(0, 5);
  ge.addEdge(1, 2);
  ge.addEdge(2, 3);
  ge.addEdge(3, 4);
  ge.addEdge(3, 5);
  ge.addEdge(4, 5);

  
  // -----
 
  */
  UndirectedClique graph(20);
  dolfin::uint* parts = new dolfin::uint[graph.numVertices()];
  dolfin::uint num_part = 6;
  GraphPartition::partition(graph, num_part, parts);

  GraphPartition::check(graph, num_part, parts);
  GraphPartition::eval(graph, num_part, parts);
  GraphPartition::disp(graph, num_part, parts);
}
void testMeshNodalGraphPartition(Mesh& mesh, dolfin::uint num_part)
{
  std::cout << "Testing Mesh partitioning" << std::endl;

  Graph graph(mesh, "nodal");
  dolfin::uint* parts = new dolfin::uint[graph.numVertices()];
  GraphPartition::partition(graph, num_part, parts);

  GraphPartition::check(graph, num_part, parts);
  GraphPartition::eval(graph, num_part, parts);
  GraphPartition::disp(graph, num_part, parts);
}
void testInputOutput()
{
  std::cout << "Testing InputOutput" << std::endl;
  Graph graph;
  GraphEditor editor;
  editor.open(graph, Graph::undirected);
  editor.initVertices(3);
  editor.addVertex(0, 2);
  editor.addVertex(1, 2);
  editor.addVertex(2, 2);
  editor.initEdges(3);
  editor.addEdge(0, 1);
  editor.addEdge(0, 2);
  editor.addEdge(1, 2);
  editor.close();
  graph.disp();

  File file("graph_test.xml");
  file << graph;

  Graph graph2;
  file >> graph2;
  graph2.disp();
}

void testCloseError()
{
  std::cout << "Testing editor closing" << std::endl;
  Graph graph;
  GraphEditor editor;
  editor.open(graph, Graph::undirected);
  editor.initVertices(3);
  editor.addVertex(0, 2);
  editor.addVertex(1, 2);
  editor.addVertex(2, 2);
  editor.initEdges(3);
  editor.addEdge(0, 1);
  editor.addEdge(0, 2);
  editor.close();
}

void testOutOfOrderError()
{
  std::cout << "Testing adding vertices out of order" << std::endl;
  Graph graph;
  GraphEditor editor;
  editor.open(graph, Graph::undirected);
  editor.initVertices(3);
  editor.addVertex(0, 2);
  editor.addVertex(2, 2);
  editor.addVertex(1, 2);
  editor.initEdges(3);
  editor.addEdge(0, 1);
  editor.addEdge(0, 2);
  editor.addEdge(1, 2);
  editor.close();
}

void testTooManyVerticesError()
{
  std::cout << "Testing adding too many vertices" << std::endl;
  Graph graph;
  GraphEditor editor;
  editor.open(graph, Graph::undirected);
  editor.initVertices(3);
  editor.addVertex(0, 2);
  editor.addVertex(1, 2);
  editor.addVertex(2, 2);
  editor.addVertex(3, 2);
  editor.initEdges(3);
  editor.addEdge(0, 1);
  editor.addEdge(0, 2);
  editor.addEdge(1, 2);
  editor.close();
}

void testInitEdgesError1()
{
  std::cout << "Testing inititializing too few edges" << std::endl;
  Graph graph;
  GraphEditor editor;
  editor.open(graph, Graph::undirected);
  editor.initVertices(4);
  editor.addVertex(0, 2);
  editor.addVertex(1, 3);
  editor.addVertex(2, 2);
  editor.addVertex(3, 2);
  editor.initEdges(3);
}

void testInitEdgesError2()
{
  std::cout << "Testing inititializing too many edges" << std::endl;
  Graph graph;
  GraphEditor editor;
  editor.open(graph, Graph::undirected);
  editor.initVertices(4);
  editor.addVertex(0, 2);
  editor.addVertex(1, 3);
  editor.addVertex(2, 2);
  editor.addVertex(3, 2);
  editor.initEdges(5);
}

int main(int argc, char* argv[])
{
  //testMeshToGraph();
  //testGraphPartition();

  UnitSquare mesh(10, 10);
  testMeshNodalGraphPartition(mesh, 9);
  //testGraphEditor();
  //testInputOutput(); 
  //testCloseError();
  //testTooManyVerticesError();
  //testInitEdgesError1();
  //testInitEdgesError2();
}
