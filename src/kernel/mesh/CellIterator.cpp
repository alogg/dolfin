#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/GenericCell.h>
#include <dolfin/CellIterator.h>
#include <dolfin/MeshData.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
// CellIterator
//-----------------------------------------------------------------------------
CellIterator::CellIterator(const Mesh& mesh)
{
  c = new MeshCellIterator(mesh);
}
//-----------------------------------------------------------------------------
CellIterator::CellIterator(const Mesh* mesh)
{
  c = new MeshCellIterator(*mesh);
}
//-----------------------------------------------------------------------------
CellIterator::CellIterator(const Cell& cell)
{
  c = new CellCellIterator(cell);
}
//-----------------------------------------------------------------------------
CellIterator::CellIterator(const CellIterator& cellIterator)
{
  c = new CellCellIterator(*cellIterator);
}
//-----------------------------------------------------------------------------
CellIterator::CellIterator(const Node& node)
{
  c = new NodeCellIterator(node);
}
//-----------------------------------------------------------------------------
CellIterator::CellIterator(const NodeIterator& nodeIterator)
{
  c = new NodeCellIterator(*nodeIterator);
}
//-----------------------------------------------------------------------------
CellIterator::~CellIterator()
{
  delete c;
}
//-----------------------------------------------------------------------------
CellIterator::operator CellPointer() const
{
  return c->pointer();
}
//-----------------------------------------------------------------------------
CellIterator::operator Cell&() const
{
  return *(*c);
}
//-----------------------------------------------------------------------------
CellIterator& CellIterator::operator++()
{
  ++(*c);

  return *this;
}
//-----------------------------------------------------------------------------
bool CellIterator::end()
{
  return c->end();
}
//-----------------------------------------------------------------------------
bool CellIterator::last()
{
  return c->last();
}
//-----------------------------------------------------------------------------
int CellIterator::index()
{
  return c->index();
}
//-----------------------------------------------------------------------------
Cell& CellIterator::operator*() const
{
  return *(*c);
}
//-----------------------------------------------------------------------------
Cell* CellIterator::operator->() const
{
  return c->pointer();
}
//-----------------------------------------------------------------------------
bool CellIterator::operator==(const CellIterator& c) const
{
  return this->c->pointer() == c.c->pointer();
}
//-----------------------------------------------------------------------------
bool CellIterator::operator!=(const CellIterator& c) const
{
  return this->c->pointer() != c.c->pointer();
}
//-----------------------------------------------------------------------------
bool CellIterator::operator==(const Cell& c) const
{
  return this->c->pointer() == &c;
}
//-----------------------------------------------------------------------------
bool CellIterator::operator!=(const Cell& c) const
{
  return this->c->pointer() != &c;
}
//-----------------------------------------------------------------------------
// CellIterator::MeshCellIterator
//-----------------------------------------------------------------------------
CellIterator::MeshCellIterator::MeshCellIterator(const Mesh &mesh)
{
  cell_iterator = mesh.gd->cells.begin();
  at_end = mesh.gd->cells.end();
}
//-----------------------------------------------------------------------------
void CellIterator::MeshCellIterator::operator++()
{
  ++cell_iterator;
}
//-----------------------------------------------------------------------------
bool CellIterator::MeshCellIterator::end()
{
  return cell_iterator == at_end;
}
//-----------------------------------------------------------------------------
bool CellIterator::MeshCellIterator::last()
{
  return cell_iterator.last();
}
//-----------------------------------------------------------------------------
int CellIterator::MeshCellIterator::index()
{
  return cell_iterator.index();
}
//-----------------------------------------------------------------------------
Cell& CellIterator::MeshCellIterator::operator*() const
{
  return *cell_iterator;
}
//-----------------------------------------------------------------------------
Cell* CellIterator::MeshCellIterator::operator->() const
{
  return cell_iterator.pointer();
}
//-----------------------------------------------------------------------------
Cell* CellIterator::MeshCellIterator::pointer() const
{
  return cell_iterator.pointer();
}
//-----------------------------------------------------------------------------
// CellIterator::NodeCellIterator
//-----------------------------------------------------------------------------
CellIterator::NodeCellIterator::NodeCellIterator(const Node& node)
{
  cell_iterator = node.nc.begin();
}
//-----------------------------------------------------------------------------
void CellIterator::NodeCellIterator::operator++()
{
  ++cell_iterator;
}
//-----------------------------------------------------------------------------
bool CellIterator::NodeCellIterator::end()
{
  return cell_iterator.end();
}
//-----------------------------------------------------------------------------
bool CellIterator::NodeCellIterator::last()
{
  return cell_iterator.last();
}
//-----------------------------------------------------------------------------
int CellIterator::NodeCellIterator::index()
{
  return cell_iterator.index();
}
//-----------------------------------------------------------------------------
Cell& CellIterator::NodeCellIterator::operator*() const
{
  return **cell_iterator;
}
//-----------------------------------------------------------------------------
Cell* CellIterator::NodeCellIterator::operator->() const
{
  return *cell_iterator;
}
//-----------------------------------------------------------------------------
Cell* CellIterator::NodeCellIterator::pointer() const
{
  return *cell_iterator;
}
//-----------------------------------------------------------------------------
// CellIterator::CellCellIterator
//-----------------------------------------------------------------------------
CellIterator::CellCellIterator::CellCellIterator(const Cell& cell)
{
  cell_iterator = cell.c->cc.begin();
}
//-----------------------------------------------------------------------------
void CellIterator::CellCellIterator::operator++()
{
  ++cell_iterator;
}
//-----------------------------------------------------------------------------
bool CellIterator::CellCellIterator::end()
{
  return cell_iterator.end();
}
//-----------------------------------------------------------------------------
bool CellIterator::CellCellIterator::last()
{
  return cell_iterator.last();
}
//-----------------------------------------------------------------------------
int CellIterator::CellCellIterator::index()
{
  return cell_iterator.index();
}
//-----------------------------------------------------------------------------
Cell& CellIterator::CellCellIterator::operator*() const
{
  return **cell_iterator;
}
//-----------------------------------------------------------------------------
Cell* CellIterator::CellCellIterator::operator->() const
{
  return *cell_iterator;
}
//-----------------------------------------------------------------------------
Cell* CellIterator::CellCellIterator::pointer() const
{
  return *cell_iterator;
}
//-----------------------------------------------------------------------------
