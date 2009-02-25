// Place for random tests

#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

class Image {};

class ImageMesh : public Mesh
{
public:
  
  ImageMesh(const Image& image) : Mesh()
  {
    // Compute nx, ny etc
    unsigned int nx = 10;
    unsigned int ny = 10;

    // Create unit square mesh and copy it
    UnitSquare mesh(nx, ny);
    static_cast<Mesh&>(*this) = mesh;
  }
};

int main()
{
  Image image;
  ImageMesh mesh(image);

  mesh.disp();
}
