"Unit tests for the MeshData class"

# Copyright (C) 2011 Anders Logg
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
# First added:  2011-08-22
# Last changed: 2011-08-22

import unittest
from dolfin import *

class MeshData(unittest.TestCase):

    def test_meshfunction(self):
        "Test input/output"

        mesh = UnitCubeMesh(3, 3, 3)

        f = mesh.data().create_mesh_function("foo")
        f.init(0)
        g = mesh.data().mesh_function("foo")

        self.assertEqual(g.size(), mesh.num_vertices())

if __name__ == "__main__":
    unittest.main()
