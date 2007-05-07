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
  editor.close();
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
  editor.close();
}

void testIteration()
{
  /* 2---3
   * |\1 |\
   * | \ | \
   * |0 \|2 \
   * 0---1---4
   *
   *
   *
   */
  
  Mesh mesh;
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

  mesh.init(2, 2);
  mesh.disp();

  for (CellIterator c0(mesh); !c0.end(); ++c0)
  {
    cout << *c0 << endl;
    for (CellIterator c1(*c0); !c1.end(); ++c1)
      cout << "  " << *c1 << endl;
  }
}

int main(int argc, char* argv[])
{
  //testMeshToGraph();
  //testGraphEditor();
  //testInputOutput(); 
  //testCloseError();
  //testTooManyVerticesError();
  //testInitEdgesError1();
  //testInitEdgesError2();

  testIteration();
}
