%feature("docstring")  dolfin::Point "

    A Point represents a point in R^3 with coordinates x, y, z, or,
    alternatively, a vector in R^3, supporting standard operations like
    the norm, distances, scalar and vector products etc.

    C++ includes: Point.h 
    
";

%feature("docstring")  dolfin::Point::Point "

        __init__(self, double x = 0.0, double y = 0.0, double z = 0.0) -> Point
        __init__(self, double x = 0.0, double y = 0.0) -> Point
        __init__(self, double x = 0.0) -> Point
        __init__(self) -> Point
        __init__(self, uint dim, double x) -> Point
        __init__(self, Point p) -> Point

        Copy constructor. 
        
";

%feature("docstring")  dolfin::Point::distance "

        distance(self, Point p) -> double

        Compute distance to given point. 
        
";

%feature("docstring")  dolfin::Point::dot "

        dot(self, Point p) -> double

        Compute dot product with given vector. 
        
";

%feature("docstring")  dolfin::Point::cross "

        cross(self, Point p) -> Point

        Compute cross product with given vector. 
        
";

%feature("docstring")  dolfin::Point::coordinates "

        coordinates(self) -> double

        Return coordinate array. 
        
";

%feature("docstring")  dolfin::Point::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::Point::y "

        y(self) -> double

        Return y-coordinate. 
        
";

%feature("docstring")  dolfin::Point::x "

        x(self) -> double

        Return x-coordinate. 
        
";

%feature("docstring")  dolfin::Point::z "

        z(self) -> double

        Return z-coordinate. 
        
";

%feature("docstring")  dolfin::Point::norm "

        norm(self) -> double

        Compute norm of point representing a vector from the origin. 
        
";

%feature("docstring")  dolfin::Mesh "

    A Mesh consists of a set of connected and numbered mesh entities.

    Both the representation and the interface are dimension-independent,
    but a concrete interface is also provided for standard named mesh
    entities:

    .. tabularcolumns:: |c|c|c|

    +--------+-----------+-------------+
    | Entity | Dimension | Codimension |
    +========+===========+=============+
    | Vertex |  0        |             |
    +--------+-----------+-------------+
    | Edge   |  1        |             | 
    +--------+-----------+-------------+
    | Face   |  2        |             | 
    +--------+-----------+-------------+
    | Facet  |           |      1      |
    +--------+-----------+-------------+
    | Cell   |           |        0    |
    +--------+-----------+-------------+

    When working with mesh iterators, all entities and connectivity
    are precomputed automatically the first time an iterator is
    created over any given topological dimension or connectivity.

    Note that for efficiency, only entities of dimension zero
    (vertices) and entities of the maximal dimension (cells) exist
    when creating a Mesh. Other entities must be explicitly created
    by calling init().

    For example, all edges in a mesh may be created by a call to mesh.init(1).
    Similarly, connectivities such as all edges connected to a given vertex
    must also be explicitly created (in this case by a call to
    mesh.init(0, 1)).
    
";

%feature("docstring")  dolfin::Mesh::ordered "


        *Returns*
            bool
                Return True iff topology is ordered according to the UFC
                numbering.
        
";

%feature("docstring")  dolfin::Mesh::smooth "

        Smooth internal vertices of mesh by local averaging.

        *Arguments*
            num_iterations
                An integer, number of iterations to perform smoothing, default
                value is 1.
        
";

%feature("docstring")  dolfin::Mesh::move "

        Move coordinates of Mesh.

        **Overloaded versions**

        * move\ **(boundary, method=hermite)**

          Move coordinates of mesh according to new boundary coordinates.

          *Arguments*
              boundary
                  A :py:class:`BoundaryMesh` instance.

              method
                  An :py:class:`ALEType` (enum).
                  Method which defines how the coordinates should be moved,
                  default is *hermite*.

        * move\ **(mesh, method=hermite)**

          Move coordinates of mesh according to adjacent mesh with common
          global vertices.

          *Arguments*
              mesh
                  A :py:class:`Mesh` instance.

              method
                  An :py:class:`ALEType` (enum).
                  Method which defines how the coordinates should be moved,
                  default is *hermite*.

        * move\ **(function)**

          Move coordinates of mesh according to displacement function. 

          *Arguments*
              function
                  A :py:class:`Function` instance.
        
";

%feature("docstring")  dolfin::Mesh::num_edges "

        Get number of edges in mesh.

        *Returns*
            integer
                Number of edges.


        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_edges()
            0
            >>> mesh.init(1)
            16
            >>> mesh.num_edges()
            16
        
";

%feature("docstring")  dolfin::Mesh::closest_point "

        Computes the point inside the mesh which is closest to the point query.

        *Arguments*
            point
                A :py:class:`Point` instance.

        *Returns*
            :py:class:`Point`
                The point inside the mesh which is closest to the point.
        
";

%feature("docstring")  dolfin::Mesh::all_intersected_entities "

        This function computes the cell ids of all cells of the current mesh
        which intersects with a given mesh entity. The result is stored in
        the last argument to the function which might be a vector or a set
        depending on which version is called.

        **Overloaded versions**

        * all_intersected_entities\ **(point, ids_result)**

          Compute all ids of all cells which are intersected by the given
          point.

          *Arguments*
              point
                  A :py:class:`Point` instance.

              ids_result
                  A set of integers.
                  The cell ids which are intersected are stored in a set for
                  efficiency reasons, to avoid to sort out duplicates later on.

        * all_intersected_entities\ **(points, ids_result)**

          Compute all ids of all cells which are intersected by any point in
          points.

          *Arguments*
              points
                  A list of :py:class:`Point` instances.

              ids_result
                  A set of integers.
                  The cell ids which are intersected are stored in a set for
                  efficiency reasons, to avoid to sort out duplicates later on.

        * all_intersected_entities\ **(entity, ids_result)**

          Compute all ids of all cells which are intersected by the given
          entity.

          *Arguments*
              entity
                  A :py:class:`MeshEntity` instance.

              ids_result
                  A list of integers.
                  The ids of the intersected cells are saved in a list. This is
                  more efficent than using a set and allows a map between the
                  (external) cell and the intersected cell of the mesh.

        * all_intersected_entities\ **(entities, ids_result)**

          Compute all id of all cells which are intersected by any entity in the
          list entities.

          *Arguments*
              entities
                  A list of :py:class:`MeshEntity` instances.

              ids_result
                  A set of integers.
                  The cell ids which are intersected are stored in a set for
                  efficiency reasons, to avoid to sort out duplicates later on.

        * all_intersected_entities\ **(another_mesh, ids_result)**

          Compute all ids of all cells which are intersected by another_mesh.

          *Arguments*
              another_mesh
                  A :py:class:`Mesh` instance.

              ids_result
                  A set of integers.
                  The cell ids which are intersected are stored in a set for
                  efficiency reasons, to avoid to sort out duplicates later on.
        
";

%feature("docstring")  dolfin::Mesh::smooth_boundary "

        Smooth boundary vertices of mesh by local averaging.

        *Arguments*
            num_iterations
                An integer, number of iterations to perform smoothing, default
                value is 1.

            harmonic_smoothing
                A bool, flag to turn on harmonics smoothing, default value is
                True.
        
";

%feature("docstring")  dolfin::Mesh::num_vertices "

        Get number of vertices in mesh.

        *Returns*
            integer
                Number of vertices.

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_vertices()
            9
        
";

%feature("docstring")  dolfin::Mesh::size "

        Get number of entities of given topological dimension.

        *Arguments*
            d
                An integer, topological dimension.

        *Returns*
            integer
                Number of entities of topological dimension d.

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.init(0,1)
            >>> mesh.num_entities(0)
            9
            >>> mesh.num_entities(1)
            16
            >>> mesh.num_entities(2)
            8
        
";

%feature("docstring")  dolfin::Mesh::coordinates "

        *Returns*
            numpy.ndarray
                Coordinates of all vertices.

        *Example*
            >>> mesh = dolfin.UnitSquare(1,1)
            >>> mesh.coordinates()
            array([[ 0.,  0.],
                   [ 1.,  0.],
                   [ 0.,  1.],
                   [ 1.,  1.]])
        
";

%feature("docstring")  dolfin::Mesh::init "

        Initialise mesh entities and connectivity.

        **Overloaded versions**

        * init\ **()**

          Compute all entities and connectivity.

        * init\ **(d)**

          Compute entities of given topological dimension.

          *Arguments*
              d
                  An integer, topological dimension.

          *Returns*
              integer
                  Number of created entities.

        * init\ **(d0, d1)**

          Compute connectivity between given pair of dimensions.

          *Arguments*
              d0
                  An integer, topological dimension.

              d1
                  An integer, topological dimension.
        
";

%feature("docstring")  dolfin::Mesh::Mesh "

        **Overloaded versions**

        * Mesh\ **()**

          Create empty mesh.

        * Mesh\ **(mesh)**

          Copy constructor.

          *Arguments*
              mesh
                  A :py:class:`Mesh` instance.

        * Mesh\ **(filename)**

          Create mesh from data file.

          *Arguments*
              filename
                  A string, name of file to load.
        
";

%feature("docstring")  dolfin::Mesh::any_intersected_entity "

        Computes only the first id  of the entity, which contains the point.

        *Arguments*
            point
                A :py:class:`Point` instance.

        *Returns*
            integer
                The first id of the cell, which contains the point, returns -1
                if no cell is intersected.
        
";

%feature("docstring")  dolfin::Mesh::type "

        *Returns*
            :py:class:`CellType`
                The cell type object associated with the mesh.
        
";

%feature("docstring")  dolfin::Mesh::num_faces "

        Get number of faces in mesh.

        *Returns*
            integer
                Number of faces.

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_faces()
            8
        
";

%feature("docstring")  dolfin::Mesh::intersection_operator "

        *Returns*
            :py:class:`IntersectionOperator`
                The intersection operator object associated with the mesh.
        
";

%feature("docstring")  dolfin::Mesh::num_entities "

        Get number of entities of given topological dimension.

        *Arguments*
            d
                An integer, topological dimension.

        *Returns*
            integer
                Number of entities of topological dimension d.

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.init(0,1)
            >>> mesh.num_entities(0)
            9
            >>> mesh.num_entities(1)
            16
            >>> mesh.num_entities(2)
            8
        
";

%feature("docstring")  dolfin::Mesh::closest_point_and_cell "

        Computes the point inside the mesh and the corresponding cell index
        which are closest to the point query.

        *Arguments*
            point
                A :py:class:`Point` instance.

        *Returns*
            std::pair<(dolfin::Point,dolfin::uint)>
                The point inside the mesh and the corresponding cell index
                which is closest to the point query.

        .. warning::

            Incomplete documentation: Don't know what the return value translates into.
        
";

%feature("docstring")  dolfin::Mesh::num_cells "

        Get number of cells in mesh.

        *Returns*
            integer
                Number of cells.

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_cells()
            8
        
";

%feature("docstring")  dolfin::Mesh::data "

        *Returns*
            :py:class:`MeshData`
                The mesh data object associated with the mesh.
        
";

%feature("docstring")  dolfin::Mesh::topology "

        *Returns*
            :py:class:`MeshTopology`
                The topology object associated with the mesh.
        
";

%feature("docstring")  dolfin::Mesh::closest_cell "

        Computes the index of the cell in the mesh which is closest to the
        point query.

        *Arguments*
            point
                A :py:class:`Point` instance.

        *Returns*
            integer
                The index of the cell in the mesh which is closest to point.

        *Example*
            >>> mesh = dolfin.UnitSquare(1,1)
            >>> point = dolfin.Point(0.0, 2.0)
            >>> mesh.closest_cell(point)
            1
        
";

%feature("docstring")  dolfin::Mesh::hmax "

        Compute maximum cell diameter.

        *Returns*
            float
                The maximum cell diameter, the diameter is computed as two
                times the circumradius (http://mathworld.wolfram.com).

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.hmax()
            0.70710678118654757
        
";

%feature("docstring")  dolfin::Mesh::geometry "

        *Returns*
            :py:class:`MeshGeometry`
                The geometry object associated with the mesh.
        
";

%feature("docstring")  dolfin::Mesh::num_facets "

        Get number of facets in mesh.

        *Returns*
            integer
                Number of facets.

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_facets()
            0
            >>> mesh.init(0,1)
            >>> mesh.num_facets()
            16
        
";

%feature("docstring")  dolfin::Mesh::clear "

        Clear all mesh data. 
        
";

%feature("docstring")  dolfin::Mesh::snap_boundary "

        Snap boundary vertices of mesh to match given sub domain.

        *Arguments*
            sub_domain
                A :py:class:`SubDomain` instance.

            harmonic_smoothing
                A bool, flag to turn on harmonics smoothing, default value is
                True.
        
";

%feature("docstring")  dolfin::Mesh::order "

        Order all mesh entities (not needed if 'mesh order entities' is
        set).

        .. seealso::

            UFC documentation (put link here!)
        
";

%feature("docstring")  dolfin::Mesh::hmin "

        Compute minimum cell diameter.

        *Returns*
            float
                The minimum cell diameter, the diameter is computed as two
                times the circumradius (http://mathworld.wolfram.com).

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.hmin()
            0.70710678118654757
        
";

%feature("docstring")  dolfin::Mesh::str "

        Informal string representation.

        *Arguments*
            verbose
                A bool, flag to turn on additional output.

        *Returns*
            string
                An informal representation of the mesh.

        *Example*
            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.str(False)
            '<Mesh of topological dimension 2 (triangles) with 9 vertices and 8 cells, ordered>'
        
";

%feature("docstring")  dolfin::Mesh::cells "

        *Returns*
            numpy.ndarray
                Connectivity for all cells.

        *Example*
            >>> mesh = dolfin.UnitSquare(1,1)
            >>> mesh.coordinates()
            array([[0, 1, 3],
                   [0, 2, 3]])
        
";

%feature("docstring")  dolfin::Variable "

    Common base class for DOLFIN variables.

    C++ includes: Variable.h 
    
";

%feature("docstring")  dolfin::Variable::rename "

        rename(self, string name, string label)

        Rename variable. 
        
";

%feature("docstring")  dolfin::Variable::disp "

        .. warning::

            Deprecated, to be removed.
        
";

%feature("docstring")  dolfin::Variable::name "

        name(self) -> string

        Return name. 
        
";

%feature("docstring")  dolfin::Variable::label "

        label(self) -> string

        Return label (description). 
        
";

%feature("docstring")  dolfin::Variable::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::Variable::Variable "

        __init__(self) -> Variable
        __init__(self, string name, string label) -> Variable

        Create variable with given name and label. 
        
";

