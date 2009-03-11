#include <dolfin.h>
#include <vector>
#include <iostream>

template<typename X> void test(X& x, const char* filename)
{
  dolfin::File outfile(filename, true); outfile << x;

  x.clear();
  dolfin::File infile(filename, true);
  infile >> x;

  std::cout << "File content:" << std::endl;
  std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
}

template<typename M> void test2(M& map, const char* filename)
{
  dolfin::File outfile(filename, true);
  outfile << map;

  map.clear();
  dolfin::File infile(filename, true);
  infile >> map;

  typedef typename M::iterator It;

  std::cout << "File content:" << std::endl;
  for (It iter = map.begin(); iter != map.end(); ++iter)
    std::cout << (*iter).first << ": " << (*iter).second << std::endl;
}

template<typename M> void test3(M& map, const char* filename)
{
  dolfin::File outfile(filename, true);
  outfile << map;

  map.clear();
  dolfin::File infile(filename, true);
  infile >> map;

  typedef typename M::iterator It;

  std::cout << "File content:" << std::endl;
  for (It iter = map.begin(); iter != map.end(); ++iter)
  {
    std::cout << "Key: " << (*iter).first << ", values: ";
    const dolfin::uint size = (*iter).second.size();
      for (dolfin::uint i = 0; i < size; ++i)
        std::cout << (*iter).second[i] << " ";
    std::cout << std::endl;
  }


}


int main()
{

  {
    std::vector<int> x;
    x.push_back(-10);
    x.push_back(-5);
    x.push_back(5);
    test(x, "int_array.xml");
  }

  {
    std::vector<dolfin::uint> x;
    x.push_back(1);
    x.push_back(2);
    x.push_back(3);
    test(x, "uint_array.xml");
  }

  {
    std::vector<double> x;
    x.push_back(1.1);
    x.push_back(2.2);
    x.push_back(3.3);
    test(x, "double_array.xml");
  }

  {
    std::map<dolfin::uint, int> map;
    map[0] = -10;
    map[1] = -5;
    map[2] = 5;
    test2(map, "int_map.xml");
  }

  {
    std::map<dolfin::uint, dolfin::uint> map;
    map[0] = 1;
    map[1] = 3;
    map[2] = 7;
    test2(map, "uint_map.xml");
  }

  {
    std::map<dolfin::uint, double> map;
    map[0] = 1.1;
    map[1] = 3.3;
    map[2] = 7.77;
    test2(map, "double_map.xml");
  }

  {
    std::map<dolfin::uint, std::vector<int> > array_map;
    array_map[0].push_back(-10);
    array_map[0].push_back(-5);
    array_map[0].push_back(5);
    array_map[1].push_back(-11);
    array_map[1].push_back(-6);
    array_map[1].push_back(6);
    test3(array_map, "int_array_map.xml");
  }

  {
    std::map<dolfin::uint, std::vector<dolfin::uint> > array_map;
    array_map[0].push_back(1);
    array_map[0].push_back(2);
    array_map[0].push_back(3);
    array_map[1].push_back(11);
    array_map[1].push_back(6);
    array_map[1].push_back(16);
    test3(array_map, "uint_array_map.xml");
  }

  {
    std::map<dolfin::uint, std::vector<double> > array_map;
    array_map[0].push_back(1.1);
    array_map[0].push_back(2.2);
    array_map[0].push_back(3.3);
    array_map[1].push_back(4.3);
    array_map[1].push_back(5.5);
    array_map[1].push_back(10.01);
    test3(array_map, "double_array_map.xml");
  }

  {
    dolfin::UnitSquare mesh(5,5);
    dolfin::File f("mesh.xml", true);
    f << mesh;
  }

  {
    dolfin::Mesh mesh;
    dolfin::File f("mesh.xml", true);
    f >> mesh;
    dolfin::File f2("mesh_copy.xml", true);
    f2 << mesh;

    dolfin::BoundaryMesh bmesh(mesh);
    dolfin::uint size = bmesh.numVertices();
    bmesh.data().create_array("spam", size);
    std::vector<dolfin::uint>* arr = bmesh.data().array("spam");
    for (dolfin::uint i = 0; i < size; ++i)
      (*arr)[i] = i+2;
    dolfin::File fb("bmesh.xml", true);
    fb << bmesh;

    dolfin::Mesh bmesh_copy;
    dolfin::File fb2("bmesh.xml", true);
    fb2 >> bmesh_copy;

    // Write file to standard error :)
    dolfin::File fb3(std::cerr);
    fb3 << bmesh_copy;
  }

  return 0;
}
