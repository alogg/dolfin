// Copyright (C) 2002 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#ifndef __TETRAHEDRON_HH
#define __TETRAHEDRON_HH

#include "Cell.hh"

class Tetrahedron : public Cell{
public:

  void Set(int n1, int n2, int n3, int n4, int material);

  int  GetSize();
  int  GetNode(int node);

  real ComputeVolume       (Grid *grid);
  real ComputeCircumRadius (Grid *grid);
  real ComputeCircumRadius (Grid *grid, real volume);
  
  /// Give access to the special functions below
  friend class Grid;
  friend class Node;
  
private:

  void CountCell            (Node *node_list);
  void AddCell              (Node *node_list, int *current, int thiscell);
  void AddNodes             (int exclude_node, int *new_nodes, int *pos);
  void ComputeCellNeighbors (Node *node_list, int thiscell);
  
  int nodes[4];
};

#endif
