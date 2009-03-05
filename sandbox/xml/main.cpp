#include <dolfin.h>
#include <vector>
#include <iostream>

template<typename X> void test(X& x, const char* filename)
{
  {
    dolfin::File _f(filename, true);
    _f << x;
  }
 
  x.clear();

  {
    dolfin::File _f(filename, true);
    _f >> x;
  }

  std::cout << "File content:" << std::endl;
  std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl;

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



  std::map<dolfin::uint, std::vector<int> > iam;
  dolfin::File f4("int_array_map.xml", true);
  f4 >> iam;


  std::cout << iam[0][0] << std::endl;
  std::cout << iam[0][1] << std::endl;
  std::cout << iam[1][0] << std::endl;
  std::cout << iam[1][1] << std::endl;
  std::cout << iam[1][2] << std::endl;

  return 0;
}
