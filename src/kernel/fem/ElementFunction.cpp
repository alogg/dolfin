#include <iostream>

#include <dolfin/ShapeFunction.h>
#include <dolfin/Product.h>
#include <dolfin/ElementFunction.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction()
{
  n = 1;
  
  a = new real[n];
  v = new Product[n]();

  a[0] = 0.0;
}

//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(real a)
{
  n = 1;

  this->a = new real[n];
  v = new Product[n];

  this->a[0] = a;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const ShapeFunction &v)
{
  n = 1;

  a = new real[n];
  this->v = new Product[n](v);

  a[0] = 1.0;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const Product &v)
{
  n = 1;

  a = new real[n];
  this->v = new Product[n](v);

  a[0] = 1.0;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(real a, const ShapeFunction &v)
{
  n = 1;

  this->a = new real[n];
  this->v = new Product[n](v);

  this->a[0] = a;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const ElementFunction &v)
{
  n = v.n;

  a = new real[n];
  this->v = new Product[n];

  for (int i = 0; i < n; i++) {
	 a[i] = v.a[i];
	 this->v[i] = v.v[i];
  }
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(real a, const Product &v)
{
  n = 1;

  this->a = new real[n];
  this->v = new Product[n](v);

  this->a[0] = a;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(real a, const ElementFunction &v)
{
  n = v.n;
  
  this->a = new real[n];
  this->v = new Product[n];
  
  for (int i = 0; i < n; i++) {
	 this->a[i] = a * v.a[i];
	 this->v[i] = v.v[i];
  }
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction
(real a0, const ShapeFunction &v0, real a1, const ShapeFunction &v1)
{
  n = 2;

  a = new real[n];
  a[0] = a0;
  a[1] = a1;

  v = new Product[n];
  v[0] = v0;
  v[1] = v1;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction
(real a0, const Product &v0, real a1, const Product &v1)
{
  n = 2;

  a = new real[n];
  a[0] = a0;
  a[1] = a1;

  v = new Product[n];
  v[0] = v0;
  v[1] = v1;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction
(real a0, const ElementFunction &v0, real a1, const ElementFunction &v1)
{
  n = v0.n + v1.n;

  a = new real[n];
  v = new Product[n];

  for (int i = 0; i < v0.n; i++) {
	 a[i] = v0.a[i];
	 v[i] = v0.v[i];
  }

  for (int i = 0; i < v1.n; i++) {
	 a[v0.n + i] = v1.a[i];
	 v[v0.n + i] = v1.v[i];
  }
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction
(real a0, const ShapeFunction &v0, real a1, const Product &v1)
{
  n = 2;

  a = new real[n];
  a[0] = a0;
  a[1] = a1;

  v = new Product[n];
  v[0] = v0;
  v[1] = v1;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction
(real a0, const ShapeFunction &v0, real a1, const ElementFunction &v1)
{
  n = 1 + v1.n;

  a = new real[n];
  v = new Product[n];
  
  a[0] = a0;
  v[0] = v0;

  for (int i = 0; i < v1.n; i++) {
	 a[1 + i] = a1 * v1.a[i];
	 v[1 + i] = v1.v[i];
  }
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction
(real a0, const Product &v0, real a1, const ElementFunction &v1)
{
  n = 1 + v1.n;

  a = new real[n];
  v = new Product[n];

  a[0] = a0;
  v[0] = v0;

  for (int i = 0; i < v1.n; i++) {
	 a[1 + i] = a1 * v1.a[i];
	 v[1 + i] = v1.v[i];
  }
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const ShapeFunction &v0,
																const ShapeFunction &v1)
{
  n = 1;

  a = new real[n];
  v = new Product[n](v0, v1);

  a[0] = 1.0;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const Product &v0,
																const Product &v1)
{
  n = 1;

  a = new real[n];
  v = new Product[n](v0, v1);

  a[0] = 1.0;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const ElementFunction &v0,
																const ElementFunction &v1)
{
  n = v0.n + v1.n;

  a = new real[n];
  v = new Product[n];

  for (int i = 0; i < v0.n; i++)
	 for (int j = 0; j < v1.n; j++) {
		a[i*v0.n + j] = v0.a[i] * v1.a[j];
		v[i*v0.n + j].set(v0.v[i],v1.v[j]);
	 }	 
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const ShapeFunction &v0,
																const Product &v1)
{
  n = 1;

  a = new real[n];
  v = new Product[n](v0, v1);

  a[0] = 1.0;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const ShapeFunction   &v0,
																const ElementFunction &v1)
{
  n = v1.n;

  a = new real[n];
  v = new Product[n];

  for (int i = 0; i < n; i++) {
	 a[i] = v1.a[i];
	 v[i].set(v0, v1.v[i]); 
  }
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::ElementFunction(const Product         &v0,
																const ElementFunction &v1)
{
  n = v1.n;

  a = new real[n];
  v = new Product[n];

  for (int i = 0; i < n; i++) {
	 a[i] = v1.a[i];
	 v[i].set(v0, v1.v[i]);
  }
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction::~ElementFunction()
{
  if ( n > 0 ) {
	 delete [] a;
	 delete [] v;
  }
}
//-----------------------------------------------------------------------------
real FunctionSpace::ElementFunction::operator()
(real x, real y, real z, real t)
{
  real value = 0.0;
  
  for (int i = 0; i < n; i++)
	 value += a[i] * v[i](x,y,z,t);

  return value;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction& FunctionSpace::ElementFunction::operator=
(real a)
{
  if ( n > 0 ) {
	 delete [] this->a;
	 delete [] v;
  }

  n = 1;

  this->a = new real[n];
  v = new Product[n]();

  this->a[0] = a;

  return *this;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction& FunctionSpace::ElementFunction::operator=
(const ShapeFunction &v)
{
  if ( n > 0 ) {
	 delete [] a;
	 delete [] this->v;
  }

  n = 1;

  a = new real[n];
  this->v = new Product[n](v);

  a[0] = 1.0;

  return *this;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction& FunctionSpace::ElementFunction::operator=
(const Product &v)
{
  if ( n > 0 ) {
	 delete [] a;
	 delete [] this->v;
  }

  n = 1;

  a = new real[n];
  this->v = new Product[n](v);

  a[0] = 1.0;

  return *this;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction& FunctionSpace::ElementFunction::operator=
(const ElementFunction &v)
{
  if ( n > 0 ) {
	 delete [] a;
	 delete [] this->v;
  }
  
  n = v.n;

  a = new real[n];
  this->v = new Product[n];

  for (int i = 0; i < n; i++) {
	 a[i] = v.a[i];
	 this->v[i] = v.v[i];
  }

  return *this;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction FunctionSpace::ElementFunction::operator+
(const ShapeFunction &v) const
{
  ElementFunction w(1.0, v, 1.0, *this);
  return w;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction FunctionSpace::ElementFunction::operator+
(const Product &v) const
{
  ElementFunction w(1.0, v, 1.0, *this);
  return w;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction FunctionSpace::ElementFunction::operator+
(const ElementFunction &v) const
{
  ElementFunction w(1.0, v, 1.0, *this);
  return w;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction FunctionSpace::ElementFunction::operator-
(const ShapeFunction &v) const
{
  ElementFunction w(-1.0, v, 1.0, *this);
  return w;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction FunctionSpace::ElementFunction::operator-
(const Product &v) const
{
  ElementFunction w(-1.0, v, 1.0, *this);
  return w;
}
//-----------------------------------------------------------------------------
FunctionSpace::ElementFunction FunctionSpace::ElementFunction::operator-
(const ElementFunction &v) const
{
  ElementFunction w(1.0, *this, -1.0, v);
  return w;
}
//-----------------------------------------------------------------------------
