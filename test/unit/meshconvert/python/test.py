""" Tests for the meshconvert module.
"""
from unittest import TestCase as _TestCase
import unittest
import os
import tempfile

from dolfin_utils import meshconvert
from dolfin_utils.meshconvert import DataHandler

class TestCase(_TestCase):
    def _get_tempfname(self, suffix=None):
        fd, fname = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        return fname

class _TestHandler(DataHandler):
    def __init__(self, cell_type, mesh_dim, test):
        DataHandler.__init__(self)
        self.vertices, self.cells, self.functions = [], [], {}
        self.vertices_ended = self.cells_ended = self.closed = False
        self.__type, self.__dim = cell_type, mesh_dim
        self.__test = test
        test.assertEqual(self._state, self.State_Invalid)
        self.test_warnings = []

    def set_mesh_type(self, *args):
        DataHandler.set_mesh_type(self, *args)
        test = self.__test
        test.assertEqual(self._state, self.State_Init)
        test.assertEqual(self._cell_type, self.__type)
        test.assertEqual(self._dim, self.__dim)

    def start_vertices(self, num_vertices):
        DataHandler.start_vertices(self, num_vertices)
        self.vertices = []
        for i in range(num_vertices):
            self.vertices.append(None)

    def add_vertex(self, vertex, coords):
        DataHandler.add_vertex(self, vertex, coords)
        self.vertices[vertex] = tuple(coords)

    def end_vertices(self):
        DataHandler.end_vertices(self)
        self.vertices_ended = True

    def start_cells(self, num_cells):
        DataHandler.start_cells(self, num_cells)
        for i in range(num_cells):
            self.cells.append(None)

    def add_cell(self, cell, nodes):
        DataHandler.add_cell(self, cell, nodes)
        self.cells[cell] = tuple(nodes)

    def end_cells(self):
        DataHandler.end_cells(self)
        self.cells_ended = True

    def start_meshfunction(self, name, dim, sz):
        DataHandler.start_meshfunction(self, name, dim, sz)
        entries = []
        for i in range(sz):
            entries.append(None)
        self.__curfunc = self.functions[name] = [dim, sz, entries, False]

    def add_entity_meshfunction(self, index, value):
        DataHandler.add_entity_meshfunction(self, index, value)
        self.__curfunc[2][index] = value

    def end_meshfunction(self):
        DataHandler.end_meshfunction(self)
        self.__curfunc[-1] = True

    def warn(self, msg):
        self.test_warnings.append(msg)

    def close(self):
        DataHandler.close(self)
        self.__test.assertEqual(self._state, self.State_Invalid)
        self.closed = True

class _ConverterTest(TestCase):
    """ Base converter test class.
    """

class AbaqusTest(_ConverterTest):
    """ Test AbaqusConverter.
    """
    def test_success(self):
        """ Test success case.
        """
        handler = self.__convert("abaqus.inp")
        # Verify vertices
        self.assertEqual(handler.vertices, [
            (0, 0, 0),
            (1, 0., 0),
            (10, 10, 10),
            (10, 10, 11),
            (0, 1, 0),
            (1, 1, 0),
            (0., 0., 0),
            (1., 0., 0),
            (0., 1., 0),
            (1., 1., 0),
            (10., 10., 10),
            (10., 10., 11),
            ])
        self.assert_(handler.vertices_ended)

        # Verify cells
        self.assertEqual(handler.cells, [
            (0, 1, 4, 5),
            (0, 1, 2, 3),
            (4, 5, 2, 3),
            ])
        self.assert_(handler.cells_ended)

        # Verify materials
        self.assertEqual(handler.functions.keys(), ["material"])
        dim, sz, entries, ended = handler.functions["material"]
        self.assertEqual(dim, 3)
        self.assertEqual(sz, 2)
        # Cell 0 should have material 0, cell 1 material 1
        self.assertEqual(entries, [0, 1])
        self.assert_(ended)

        self.assert_(handler.closed)

    def test_error(self):
        """ Test various cases of erroneus input.
        """
        def convert(fname, text, error=False):
            f = file(fname, "w")
            f.write(text)
            f.close()
            if not error:
                handler = self.__convert(fname)
                self.assert_(handler.test_warnings)
            else:
                self.assertRaises(meshconvert.ParseError, self.__convert, fname)
            os.remove(fname)

        fname = self._get_tempfname(suffix=".inp")

        # Too few coordinates
        convert(fname, """*NODE
1, 0, 0""")
        # Non-numeric index
        convert(fname, """*NODE
a, 0, 0, 0""")
        # Non-numeric coordinate
        convert(fname, """*NODE
1, 0, 0, a""")
        # Unsupported element type, also the body should be ignored
        convert(fname, """*ELEMENT, TYPE=sometype
0
""")
        # Bad parameter syntax
        convert(fname, "*ELEMENT, TYPE=sometype, BAD")
        # Missing type specifier
        convert(fname, "*ELEMENT", error=True)
        # Non-existent node
        convert(fname, """*NODE
1, 0, 0, 0
2, 0, 0, 0
3, 0, 0, 0
*ELEMENT, TYPE=C3D4
1, 1, 2, 3, 4
""", error=True)
        # Too few nodes
        convert(fname, """*NODE
1, 0, 0, 0
2, 0, 0, 0
3, 0, 0, 0
*ELEMENT, TYPE=C3D4
1, 1, 2, 3
""", error=True)
        # Non-existent element set
        convert(fname, """*MATERIAL, NAME=MAT
*SOLID SECTION, ELSET=NONE, MATERIAL=MAT""", error=True)

    def __convert(self, fname):
        handler = _TestHandler(DataHandler.CellType_Tetrahedron, 3, self)
        if not os.path.isabs(fname):
            fname = os.path.join("data", fname)
        meshconvert.convert(fname, handler)
        return handler

class GmshTest(_ConverterTest):
    """ Test Gmsh convertor.
    """
    def test_success(self):
        """ Test success case.
        """
        handler = self.__convert("gmsh.msh")
        # Verify vertices
        self.assertEqual(handler.vertices, [
            (0,    0,    0),
            (1,    0,    1),
            (-0,    0.8,  0.6),
            (0.3,  0.8, -0.1),
            (0.6,  0.3, -0.4),
            (0.5,  0,    0.5),
            (0.5,  0.4,  0.8),
            (0.76, 0.26, 0.63),
            (0.53, 0.53, 0.26),
            (0.8,  0.15, 0.3)
            ])
        self.assert_(handler.vertices_ended)

        # Verify cells
        self.assertEqual(handler.cells, [
            (9, 5, 1, 7),
            (4, 8, 0, 9),
            (8, 5, 0, 9),
            (8, 5, 9, 7),
            (8, 3, 4, 0),
            (8, 2, 3, 0),
            (7, 5, 6, 8),
            (5, 2, 6, 8),
            (1, 7, 5, 6),
            (5, 2, 8, 0)
            ])
        self.assert_(handler.cells_ended)

        # Verify physical regions
        self.assertEqual(handler.functions.keys(), ["physical_region"])
        dim, sz, entries, ended = handler.functions["physical_region"]
        self.assertEqual(dim, 3)
        self.assertEqual(sz, 10)        # There are 10 cells
        # Cells 0 thru 4 should be in region 1000, cells 5 thru 9 in
        # region 2000
        self.assertEqual(entries, [1000]*5 + [2000]*5)
        self.assert_(ended)
        self.assert_(handler.closed)
    
    def test_2D_facet_markings_1 (self):
        """
        Test to see if the 2D facet markings behave as expected.
        1 edge marked
        """
        marked_facets = [7]
        self._facet_marker_driver(2, 1, marked_facets, 8)
    
    def test_2D_facet_markings_2 (self):
        """
        Test to see if the 2D facet markings behave as expected.
        2 edges marked
        """
        marked_facets = [2,5]
        self._facet_marker_driver(2, 2, marked_facets, 8)
        
    def test_2D_facet_markings_3 (self):
        """
        Test to see if the 2D facet markings behave as expected.
        3 edges marked
        """
        marked_facets = [5,6,7]
        self._facet_marker_driver(2, 3, marked_facets, 8)
        
    def test_2D_facet_markings_4 (self):
        """
        Test to see if the 2D facet markings behave as expected.
        4 edges marked
        """
        marked_facets = [2,5,6,7]
        self._facet_marker_driver(2, 4, marked_facets, 8)
    
    def test_3D_facet_markings_1 (self):
        """
        Test the marking of 3D facets
        Unit cube, 1 Face marked
        """
#         [0, 0, 0, 999, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 999, 0, 0, 0, 0, 0, 0, 0, 0, 999, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 999, 
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        marked_facets = [3, 15, 24, 39,]
        self._facet_marker_driver(3, 1, marked_facets, 60)
    
    def _facet_marker_driver (self, dim, id, marked_facets, size ):
        if dim == 2:
            cell_type = DataHandler.CellType_Triangle
        elif dim == 3:
            cell_type = DataHandler.CellType_Tetrahedron
            
        handler = self.__convert("gmsh_test_facet_regions_%dD_%d.msh" % (dim, id), cell_type, dim)
        
        free_facets = range(size)

        for i in marked_facets:
            free_facets.remove(i)
        
        function_dim, sz, entries, ended = handler.functions["facet_region"]
        
        # the dimension of the meshfunction should be 1
        self.assertEqual(function_dim, dim-1)
        # There should be 8 edges in the mesh function
        self.assertEqual(len(entries), size)
        self.assertEqual(sz, size)
        # marked
        self.assert_( all ( entries[i] == 999 for i in marked_facets ) )
        # all other edges should be zero
        self.assert_( all ( entries[i] == 0 for i in free_facets ) )
        
        self.assert_(ended)
        self.assert_(handler.closed)
        
    def __convert(self, fname, cell_type=DataHandler.CellType_Tetrahedron, dim=3):
        handler = _TestHandler(cell_type, dim, self)
        if not os.path.isabs(fname):
            fname = os.path.join("data", fname)
        meshconvert.convert(fname, handler)
        return handler

if __name__ == "__main__":
    unittest.main()
