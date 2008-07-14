// Copyright (c) 2008 Nuno Lopes
// Licensed under the GNU LGPL Version 2.1
//
// First added:  2008-06-20


#include <dolfin.h>

using namespace dolfin;

int main()
{

  // rectangle (-1,2)x(-0.5,0.5) with nx=20 and ny=10
  Rectangle mesh1(-1.,3.,-1.,1.,20,10);
  File fmesh1("rectangle.pvd");
  fmesh1<<mesh1;

  //Rectagular prysm or block (2,10)x(1,3)x(0,2) with nx=70,
  // ny=30 and nx=15
  Box mesh2(2.0,10.0,1.0,3.0,0.0,2.0,70,30,15);
  File fmesh2("box.pvd");
  fmesh2<<mesh2;
  
  //UnitCircle
  UnitCircle mesh3(40);
  for (dolfin::uint i=0;i<5;i++)
    mesh3.smooth();
  File fmesh3("circle.pvd");
  fmesh3<<mesh3;

 
  //UnitSphere
  UnitSphere mesh4(20);
  File fmesh4("sphere.pvd");
  for (dolfin::uint i=0;i<5;i++)
    mesh4.smooth();
  fmesh4<<mesh4;

  
 
  return (EXIT_SUCCESS);
}
