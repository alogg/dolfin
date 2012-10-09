"Unit tests for BoundaryMesh library"

# Copyright (C) 2012 Garth N. Wells
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2011-10-09
# Last changed:

import unittest
import numpy
from dolfin import *

class BoundaryMeshConstruction(unittest.TestCase):

    def xtest_1D_mesh(self):
        mesh = UnitInterval(10)

    def test_2D_mesh(self):
        mesh = UnitSquare(4, 4)

        # Create global boundary mesh
        bmesh1 = BoundaryMesh()
        bmesh1.init_exterior_boundary(mesh)
        self.assertEqual(MPI.sum(bmesh1.num_cells()), 16)
        self.assertEqual(bmesh1.topology().dim(), 1)

if __name__ == "__main__":
    unittest.main()
