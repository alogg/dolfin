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
    std::vector<uint> x;
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
    std::map<dolfin::uint, uint> map;
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


  /*
  std::map<dolfin::uint, std::vector<int> > iam;
  iam[0].push_back(-10);
  iam[0].push_back(-5);
  iam[0].push_back(5);


  std::cout << iam[0][0] << std::endl;
  std::cout << iam[0][1] << std::endl;
  std::cout << iam[1][0] << std::endl;
  std::cout << iam[1][1] << std::endl;
  std::cout << iam[1][2] << std::endl;
  */

  return 0;
}
