%feature("docstring")  dolfin::TableEntry "

    This class represents an entry in a Table.

    C++ includes: Table.h 
    
";

%feature("docstring")  dolfin::TableEntry::TableEntry "

        __init__(self, string row, string col, Table table) -> TableEntry

        Create table entry. 
        
";

%feature("docstring")  dolfin::KrylovSolver "

    This class defines an interface for a Krylov solver. The approproiate
    solver is chosen on the basis of the matrix/vector type.

    C++ includes: KrylovSolver.h 
    
";

%feature("docstring")  dolfin::KrylovSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::KrylovSolver::solve "

        solve(self, GenericVector x, GenericVector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint

        Solve linear system Ax = b. 
        
";

%feature("docstring")  dolfin::KrylovSolver::KrylovSolver "

        __init__(self, string solver_type = "default", string pc_type = "default") -> KrylovSolver
        __init__(self, string solver_type = "default") -> KrylovSolver
        __init__(self) -> KrylovSolver

        Create Krylov solver. 
        
";

%feature("docstring")  dolfin::MeshConnectivity "

    Mesh connectivity stores a sparse data structure of connections
    (incidence relations) between mesh entities for a fixed pair of
    topological dimensions.

    The connectivity can be specified either by first giving the number of
    entities and the number of connections for each entity, which may
    either be equal for all entities or different, or by giving the entire
    (sparse) connectivity pattern.

    C++ includes: MeshConnectivity.h 
    
";

%feature("docstring")  dolfin::MeshConnectivity::set "

        set(self, uint entity, uint connection, uint pos)
        set(self, uint entity, std::vector<(dolfin::uint)> connections)
        set(self, uint entity, uint connections)
        set(self, std::vector<(std::vector<(dolfin::uint)>)> connectivity)

        Set all connections for all entities. 
        
";

%feature("docstring")  dolfin::MeshConnectivity::clear "

        clear(self)

        Clear all data. 
        
";

%feature("docstring")  dolfin::MeshConnectivity::init "

        init(self, uint num_entities, uint num_connections)
        init(self, std::vector<(dolfin::uint)> num_connections)

        Initialize number of entities and number of connections
        (individually). 
        
";

%feature("docstring")  dolfin::MeshConnectivity::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::MeshConnectivity::MeshConnectivity "

        __init__(self, uint d0, uint d1) -> MeshConnectivity
        __init__(self, MeshConnectivity connectivity) -> MeshConnectivity

        Copy constructor. 
        
";

%feature("docstring")  dolfin::MeshConnectivity::size "

        size(self) -> uint
        size(self, uint entity) -> uint

        Return number of connections for given entity. 
        
";

%feature("docstring")  dolfin::MeshGeometry "
Proxy of C++ dolfin::MeshGeometry class
";

%feature("docstring")  dolfin::MeshGeometry::set "

        set(self, uint n, uint i, double x)

        Set value of coordinate n in direction i. 
        
";

%feature("docstring")  dolfin::MeshGeometry::point "

        point(self, uint n) -> Point

        Return coordinate n as a 3D point value. 
        
";

%feature("docstring")  dolfin::MeshGeometry::init_higher_order_cells "

        init_higher_order_cells(self, uint num_cells, uint num_dof)

        Initialize higher order cell data list to given number of cells and
        dofs. 
        
";

%feature("docstring")  dolfin::MeshGeometry::num_higher_order_vertices_per_cell "

        num_higher_order_vertices_per_cell(self) -> uint

        Return number of vertices used (per cell) to represent the higher
        order geometry. 
        
";

%feature("docstring")  dolfin::MeshGeometry::affine_cell_bool "

        affine_cell_bool(self) -> bool

        Return pointer to boolean affine indicator array. 
        
";

%feature("docstring")  dolfin::MeshGeometry::MeshGeometry "

        __init__(self) -> MeshGeometry
        __init__(self, MeshGeometry geometry) -> MeshGeometry

        Copy constructor. 
        
";

%feature("docstring")  dolfin::MeshGeometry::size "

        size(self) -> uint

        Return number of coordinates. 
        
";

%feature("docstring")  dolfin::MeshGeometry::dim "

        dim(self) -> uint

        Return Euclidean dimension of coordinate system. 
        
";

%feature("docstring")  dolfin::MeshGeometry::higher_order_cell "

        higher_order_cell(self, uint c) -> uint
        higher_order_cell(self, uint c) -> uint

        Return array of higher order vertex indices for a specific higher
        order cell. 
        
";

%feature("docstring")  dolfin::MeshGeometry::set_higher_order_cell_data "

        set_higher_order_cell_data(self, uint N, std::vector<(dolfin::uint)> vector_cell_data)

        Set higher order cell data for cell # N in direction i. 
        
";

%feature("docstring")  dolfin::MeshGeometry::higher_order_x "

        higher_order_x(self, uint n) -> double
        higher_order_x(self, uint n) -> double
        higher_order_x(self) -> double
        higher_order_x(self) -> double

        Return array of values for all higher order coordinates. 
        
";

%feature("docstring")  dolfin::MeshGeometry::init_higher_order_vertices "

        init_higher_order_vertices(self, uint dim, uint size_higher_order)

        Initialize higher order coordinate list to given dimension and size.

        
";

%feature("docstring")  dolfin::MeshGeometry::clear "

        clear(self)

        Clear all data. 
        
";

%feature("docstring")  dolfin::MeshGeometry::init_affine_indicator "

        init_affine_indicator(self, uint num_cells)

        Initialize the affine indicator array. 
        
";

%feature("docstring")  dolfin::MeshGeometry::set_affine_indicator "

        set_affine_indicator(self, uint i, bool value)

        set affine indicator at index i 
        
";

%feature("docstring")  dolfin::MeshGeometry::init "

        init(self, uint dim, uint size)

        Initialize coordinate list to given dimension and size. 
        
";

%feature("docstring")  dolfin::MeshGeometry::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::MeshGeometry::set_higher_order_coordinates "

        set_higher_order_coordinates(self, uint N, uint i, double x)

        Set value of higher order coordinate N in direction i. 
        
";

%feature("docstring")  dolfin::MeshGeometry::higher_order_cells "

        higher_order_cells(self) -> uint
        higher_order_cells(self) -> uint

        Return array of values for all higher order cell data. 
        
";

%feature("docstring")  dolfin::MeshGeometry::x "

        x(self, uint n, uint i) -> double
        x(self, uint n) -> double
        x(self, uint n) -> double
        x(self) -> double
        x(self) -> double

        Return array of values for all coordinates. 
        
";

%feature("docstring")  dolfin::BoolParameter "

    Parameter with value type bool.

    C++ includes: Parameter.h 
    
";

%feature("docstring")  dolfin::BoolParameter::BoolParameter "

        __init__(self, string key, bool value) -> BoolParameter

        Create bool-valued parameter. 
        
";

%feature("docstring")  dolfin::BoolParameter::_assign_bool "
_assign_bool(self, bool value) -> BoolParameter
";

%feature("docstring")  dolfin::PETScLUSolver "
Proxy of C++ dolfin::PETScLUSolver class
";

%feature("docstring")  dolfin::PETScLUSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::PETScLUSolver::solve "

        solve(self, GenericVector x, GenericVector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint
        solve(self, PETScMatrix A, PETScVector x, PETScVector b) -> uint

        Solve linear system Ax = b. 
        
";

%feature("docstring")  dolfin::PETScLUSolver::ksp "
ksp(self) -> boost::shared_ptr<(KSP)>
";

%feature("docstring")  dolfin::PETScLUSolver::set_operator "

        set_operator(self, GenericMatrix A)
        set_operator(self, PETScMatrix A)

        Set operator (matrix). 
        
";

%feature("docstring")  dolfin::PETScLUSolver::PETScLUSolver "

        __init__(self, string lu_package = "default") -> PETScLUSolver
        __init__(self) -> PETScLUSolver
        __init__(self, GenericMatrix A, string lu_package = "default") -> PETScLUSolver
        __init__(self, GenericMatrix A) -> PETScLUSolver
        __init__(self, boost::shared_ptr<(q(const).dolfin::PETScMatrix)> A, 
            string lu_package = "default") -> PETScLUSolver
        __init__(self, boost::shared_ptr<(q(const).dolfin::PETScMatrix)> A) -> PETScLUSolver
        
";

%feature("docstring")  dolfin::PrimitiveIntersector "

    This class implements an intersection detection, detecting whether two
    given (arbitrary) meshentities intersect.

    C++ includes: PrimitiveIntersector.h 
    
";

%feature("docstring")  dolfin::PrimitiveIntersector::do_intersect "

        do_intersect(MeshEntity entity_1, MeshEntity entity_2) -> bool
        do_intersect(MeshEntity entity_1, Point point) -> bool
        
";

%feature("docstring")  dolfin::PrimitiveIntersector::do_intersect_exact "

        do_intersect_exact(MeshEntity entity_1, MeshEntity entity_2) -> bool
        do_intersect_exact(MeshEntity entity_1, Point point) -> bool
        
";

%feature("docstring")  dolfin::PrimitiveIntersector::PrimitiveIntersector "
__init__(self) -> PrimitiveIntersector
";

%feature("docstring")  dolfin::MeshFunctionBool "

    A MeshFunction is a function that can be evaluated at a set of mesh
    entities. A MeshFunction is discrete and is only defined at the set of
    mesh entities of a fixed topological dimension. A MeshFunction may for
    example be used to store a global numbering scheme for the entities of
    a (parallel) mesh, marking sub domains or boolean markers for mesh
    refinement.

    C++ includes: MeshFunction.h 
    
";

%feature("docstring")  dolfin::MeshFunctionBool::dim "

        dim(self) -> uint

        Return topological dimension. 
        
";

%feature("docstring")  dolfin::MeshFunctionBool::set_all "

        set_all(self, bool value)

        Set all values to given value. 
        
";

%feature("docstring")  dolfin::MeshFunctionBool::init "

        init(self, uint dim)
        init(self, uint dim, uint size)
        init(self, Mesh mesh, uint dim)
        init(self, Mesh mesh, uint dim, uint size)

        Initialize mesh function for given topological dimension of given
        size. 
        
";

%feature("docstring")  dolfin::MeshFunctionBool::mesh "

        mesh(self) -> Mesh

        Return mesh associated with mesh function. 
        
";

%feature("docstring")  dolfin::MeshFunctionBool::values "

        values(self) -> PyObject

        Return array of values. 
        
";

%feature("docstring")  dolfin::MeshFunctionBool::MeshFunctionBool "

        __init__(self) -> MeshFunctionBool
        __init__(self, Mesh mesh) -> MeshFunctionBool
        __init__(self, Mesh mesh, uint dim) -> MeshFunctionBool
        __init__(self, Mesh mesh, uint dim, bool value) -> MeshFunctionBool
        __init__(self, Mesh mesh, string filename) -> MeshFunctionBool
        __init__(self, MeshFunctionBool f) -> MeshFunctionBool

        Copy constructor. 
        
";

%feature("docstring")  dolfin::MeshFunctionBool::size "

        size(self) -> uint

        Return size (number of entities). 
        
";

%feature("docstring")  dolfin::Progress "

    This class provides a simple way to create and update progress bars
    during a computation. A progress bar may be used either in an
    iteration with a known number of steps:

    Progress p("Iterating...", n); for (int i = 0; i < n; i++) { ...
    p++; }

    or in an iteration with an unknown number of steps:

    Progress p("Iterating..."); while (t < T) { ... p = t / T; }

    C++ includes: Progress.h 
    
";

%feature("docstring")  dolfin::Progress::Progress "

        __init__(self, string title, unsigned int n) -> Progress
        __init__(self, string title) -> Progress

        Create progress bar with an unknown number of steps. 
        
";

%feature("docstring")  dolfin::Progress::update "
Missing docstring
";

%feature("docstring")  dolfin::uBLASILUPreconditioner "

    This class implements an incomplete LU factorization (ILU)
    preconditioner for the uBLAS Krylov solver.

    C++ includes: uBLASILUPreconditioner.h 
    
";

%feature("docstring")  dolfin::uBLASILUPreconditioner::uBLASILUPreconditioner "

        __init__(self, Parameters krylov_parameters) -> uBLASILUPreconditioner

        Constructor. 
        
";

%feature("docstring")  dolfin::IntParameter "

    Parameter with value type int.

    C++ includes: Parameter.h 
    
";

%feature("docstring")  dolfin::IntParameter::_assign "
_assign(self, int value) -> IntParameter
";

%feature("docstring")  dolfin::IntParameter::IntParameter "

        __init__(self, string key, int value) -> IntParameter

        Create int-valued parameter. 
        
";

%feature("docstring")  dolfin::UnitSphere "

    Triangular mesh of the 3D unit sphere.

    Given the number of cells (nx, ny, nz) in each direction, the total
    number of tetrahedra will be 6*nx*ny*nz and the total number of
    vertices will be (nx + 1)*(ny + 1)*(nz + 1).

    C++ includes: UnitSphere.h 
    
";

%feature("docstring")  dolfin::UnitSphere::UnitSphere "
__init__(self, uint nx) -> UnitSphere
";

%feature("docstring")  dolfin::ALE "

    This class provides functionality useful for implementation of ALE
    (Arbitrary Lagrangian-Eulerian) methods, in particular moving the
    boundary vertices of a mesh and then interpolating the new coordinates
    for the interior vertices accordingly.

    C++ includes: ALE.h 
    
";

%feature("docstring")  dolfin::ALE::move "

        move(Mesh mesh, BoundaryMesh new_boundary, ALEType method = lagrange)
        move(Mesh mesh, BoundaryMesh new_boundary)
        move(Mesh mesh0, Mesh mesh1, ALEType method = lagrange)
        move(Mesh mesh0, Mesh mesh1)
        move(Mesh mesh, Function displacement)
        
";

%feature("docstring")  dolfin::ALE::ALE "
__init__(self) -> ALE
";

%feature("docstring")  dolfin::Vertex "

    A Vertex is a MeshEntity of topological dimension 0.

    C++ includes: Vertex.h 
    
";

%feature("docstring")  dolfin::Vertex::x "

        x(self, uint i) -> double
        x(self) -> double

        Return array of vertex coordinates (const version). 
        
";

%feature("docstring")  dolfin::Vertex::Vertex "

        __init__(self, Mesh mesh, uint index) -> Vertex
        __init__(self, MeshEntity entity) -> Vertex

        Create vertex from mesh entity. 
        
";

%feature("docstring")  dolfin::Vertex::point "

        point(self) -> Point

        Return vertex coordinates as a 3D point value. 
        
";

%feature("docstring")  dolfin::entities "

    MeshEntityIterator provides a common iterator for mesh entities over
    meshes, boundaries and incidence relations. The basic use is
    illustrated below.

    The following example shows how to iterate over all mesh entities of a
    mesh of topological dimension dim:

    for ( MeshEntityIterator e(mesh, dim); !e. end(); ++e) { e->foo(); }

    The following example shows how to iterate over mesh entities of
    topological dimension dim connected (incident) to some mesh entity f:

    for ( MeshEntityIterator e(f, dim); !e. end(); ++e) { e->foo(); }

    In addition to the general iterator, a set of specific named iterators
    are provided for entities of type Vertex, Edge, Face, Facet and Cell.
    These iterators are defined along with their respective classes.

    C++ includes: MeshEntityIterator.h 
    
";

%feature("docstring")  dolfin::entities::end "

        end(self) -> bool

        Check if iterator has reached the end. 
        
";

%feature("docstring")  dolfin::entities::_dereference "
_dereference(self) -> MeshEntity
";

%feature("docstring")  dolfin::entities::_decrease "
_decrease(self) -> entities
";

%feature("docstring")  dolfin::entities::pos "

        pos(self) -> uint

        Return current position. 
        
";

%feature("docstring")  dolfin::entities::next "
Missing docstring
";

%feature("docstring")  dolfin::entities::end_iterator "

        end_iterator(self) -> entities

        Provide a safeguard iterator pointing beyond the end of an iteration
        process, either iterating over the mesh /or incident entities. Added
        to be bit more like STL iteratoren, since many algorithms rely on a
        kind of beyond iterator. 
        
";

%feature("docstring")  dolfin::entities::_increment "
_increment(self) -> entities
";

%feature("docstring")  dolfin::entities::entities "

        MeshEntityIterator() -> entities
        MeshEntityIterator(Mesh mesh, uint dim) -> entities
        MeshEntityIterator(MeshEntity entity, uint dim) -> entities
        __init__(self, entities it) -> entities

        Copy Constructor. 
        
";

%feature("docstring")  dolfin::uBLASDenseFactory "
Proxy of C++ dolfin::uBLASFactory<(dolfin::ublas_dense_matrix)> class
";

%feature("docstring")  dolfin::uBLASDenseFactory::create_matrix "

        create_matrix(self) -> uBLASDenseMatrix

        Create empty matrix. 
        
";

%feature("docstring")  dolfin::uBLASDenseFactory::create_vector "

        create_vector(self) -> uBLASVector

        Create empty vector. 
        
";

%feature("docstring")  dolfin::uBLASDenseFactory::instance "
instance() -> uBLASDenseFactory
";

%feature("docstring")  dolfin::uBLASDenseFactory::create_lu_solver "

        create_lu_solver(self) -> UmfpackLUSolver

        Create LU solver. 
        
";

%feature("docstring")  dolfin::uBLASDenseFactory::create_local_vector "

        create_local_vector(self) -> uBLASVector

        Create empty vector (local). 
        
";

%feature("docstring")  dolfin::uBLASDenseFactory::create_pattern "

        create_pattern(self) -> SparsityPattern

        Create empty sparsity pattern. 
        
";

%feature("docstring")  dolfin::uBLASDenseFactory::uBLASDenseFactory "
No constructor defined
";

%feature("docstring")  dolfin::MeshFunctionInt "

    A MeshFunction is a function that can be evaluated at a set of mesh
    entities. A MeshFunction is discrete and is only defined at the set of
    mesh entities of a fixed topological dimension. A MeshFunction may for
    example be used to store a global numbering scheme for the entities of
    a (parallel) mesh, marking sub domains or boolean markers for mesh
    refinement.

    C++ includes: MeshFunction.h 
    
";

%feature("docstring")  dolfin::MeshFunctionInt::dim "

        dim(self) -> uint

        Return topological dimension. 
        
";

%feature("docstring")  dolfin::MeshFunctionInt::set_all "

        set_all(self, int value)

        Set all values to given value. 
        
";

%feature("docstring")  dolfin::MeshFunctionInt::init "

        init(self, uint dim)
        init(self, uint dim, uint size)
        init(self, Mesh mesh, uint dim)
        init(self, Mesh mesh, uint dim, uint size)

        Initialize mesh function for given topological dimension of given
        size. 
        
";

%feature("docstring")  dolfin::MeshFunctionInt::mesh "

        mesh(self) -> Mesh

        Return mesh associated with mesh function. 
        
";

%feature("docstring")  dolfin::MeshFunctionInt::values "

        values(self) -> PyObject

        Return array of values. 
        
";

%feature("docstring")  dolfin::MeshFunctionInt::MeshFunctionInt "

        __init__(self) -> MeshFunctionInt
        __init__(self, Mesh mesh) -> MeshFunctionInt
        __init__(self, Mesh mesh, uint dim) -> MeshFunctionInt
        __init__(self, Mesh mesh, uint dim, int value) -> MeshFunctionInt
        __init__(self, Mesh mesh, string filename) -> MeshFunctionInt
        __init__(self, MeshFunctionInt f) -> MeshFunctionInt

        Copy constructor. 
        
";

%feature("docstring")  dolfin::MeshFunctionInt::size "

        size(self) -> uint

        Return size (number of entities). 
        
";

%feature("docstring")  dolfin::BoundaryCondition "

    Common base class for boundary conditions.

    C++ includes: BoundaryCondition.h 
    
";

%feature("docstring")  dolfin::BoundaryCondition::_function_space "

        _function_space(self) -> __dummy_19__

        Return shared pointer to function space. 
        
";

%feature("docstring")  dolfin::BoundaryCondition::apply "

        apply(self, GenericMatrix A)
        apply(self, GenericVector b)
        apply(self, GenericMatrix A, GenericVector b)
        apply(self, GenericVector b, GenericVector x)
        apply(self, GenericMatrix A, GenericVector b, GenericVector x)

        Apply boundary condition to a linear system for a nonlinear problem.

        
";

%feature("docstring")  dolfin::BoundaryCondition::function_space "
Return the FunctionSpace
";

%feature("docstring")  dolfin::BoundaryCondition::BoundaryCondition "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::PETScFactory "
Proxy of C++ dolfin::PETScFactory class
";

%feature("docstring")  dolfin::PETScFactory::create_matrix "
create_matrix(self) -> PETScMatrix
";

%feature("docstring")  dolfin::PETScFactory::create_vector "
create_vector(self) -> PETScVector
";

%feature("docstring")  dolfin::PETScFactory::create_krylov_solver "

        create_krylov_solver(self, string method, string pc) -> PETScKrylovSolver

        Create Krylov solver. 
        
";

%feature("docstring")  dolfin::PETScFactory::instance "
instance() -> PETScFactory
";

%feature("docstring")  dolfin::PETScFactory::create_lu_solver "

        create_lu_solver(self) -> PETScLUSolver

        Create LU solver. 
        
";

%feature("docstring")  dolfin::PETScFactory::create_local_vector "

        create_local_vector(self) -> PETScVector

        Create empty vector (local). 
        
";

%feature("docstring")  dolfin::PETScFactory::create_pattern "
create_pattern(self) -> SparsityPattern
";

%feature("docstring")  dolfin::PETScFactory::PETScFactory "
No constructor defined
";

%feature("docstring")  dolfin::PETScVector "
Proxy of C++ dolfin::PETScVector class
";

%feature("docstring")  dolfin::PETScVector::sum "

        sum(self) -> double
        sum(self, dolfin::Array<(dolfin::uint)> rows) -> double

        Return sum of selected rows in vector. Repeated entries only summed
        once. 
        
";

%feature("docstring")  dolfin::PETScVector::_assign "

        _assign(self, GenericVector x) -> GenericVector
        _assign(self, double a) -> PETScVector
        _assign(self, PETScVector x) -> PETScVector
        
";

%feature("docstring")  dolfin::PETScVector::get_local "

        get_local(self, double block, uint m)
        get_local(self, DoubleArray values)

        Get all values on local process. 
        
";

%feature("docstring")  dolfin::PETScVector::vec "
vec(self) -> boost::shared_ptr<(Vec)>
";

%feature("docstring")  dolfin::PETScVector::copy "
copy(self) -> PETScVector
";

%feature("docstring")  dolfin::PETScVector::PETScVector "

        __init__(self, string type = "global") -> PETScVector
        __init__(self) -> PETScVector
        __init__(self, uint N, string type = "global") -> PETScVector
        __init__(self, uint N) -> PETScVector
        __init__(self, PETScVector x) -> PETScVector
        __init__(self, boost::shared_ptr<(Vec)> x) -> PETScVector
        
";

%feature("docstring")  dolfin::PETScMatrixDeleter "
Proxy of C++ dolfin::PETScMatrixDeleter class
";

%feature("docstring")  dolfin::PETScMatrixDeleter::PETScMatrixDeleter "
__init__(self) -> PETScMatrixDeleter
";

%feature("docstring")  dolfin::MeshTopology "

    MeshTopology stores the topology of a mesh, consisting of mesh
    entities and connectivity (incidence relations for the mesh entities).
    Note that the mesh entities don't need to be stored, only the number
    of entities and the connectivity. Any numbering scheme for the mesh
    entities is stored separately in a MeshFunction over the entities.

    A mesh entity e may be identified globally as a pair e = (dim, i),
    where dim is the topological dimension and i is the index of the
    entity within that topological dimension.

    C++ includes: MeshTopology.h 
    
";

%feature("docstring")  dolfin::MeshTopology::dim "

        dim(self) -> uint

        Return topological dimension. 
        
";

%feature("docstring")  dolfin::MeshTopology::clear "

        clear(self)

        Clear all data. 
        
";

%feature("docstring")  dolfin::MeshTopology::init "

        init(self, uint dim)
        init(self, uint dim, uint size)

        Set number of entities (size) for given topological dimension. 
        
";

%feature("docstring")  dolfin::MeshTopology::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::MeshTopology::MeshTopology "

        __init__(self) -> MeshTopology
        __init__(self, MeshTopology topology) -> MeshTopology

        Copy constructor. 
        
";

%feature("docstring")  dolfin::MeshTopology::size "

        size(self, uint dim) -> uint

        Return number of entities for given dimension. 
        
";

%feature("docstring")  dolfin::dGqMethod "

    Contains all numeric constants, such as nodal points and nodal
    weights, needed for the dG(q) method. The order q must be at least 0.
    Note that q refers to the polynomial order and not the order of
    convergence for the method, which is 2q + 1.

    C++ includes: dGqMethod.h 
    
";

%feature("docstring")  dolfin::dGqMethod::ueval "

        ueval(self, real x0, real values, real tau) -> real
        ueval(self, real x0, real values, uint i) -> real

        Evaluate solution at given node (inline optimized). 
        
";

%feature("docstring")  dolfin::dGqMethod::dGqMethod "
__init__(self, unsigned int q) -> dGqMethod
";

%feature("docstring")  dolfin::PETScUserPreconditioner "
Proxy of C++ dolfin::PETScUserPreconditioner class
";

%feature("docstring")  dolfin::PETScUserPreconditioner::setup "
setup(KSP ksp, PETScUserPreconditioner pc)
";

%feature("docstring")  dolfin::PETScUserPreconditioner::solve "
solve(self, PETScVector x, PETScVector b)
";

%feature("docstring")  dolfin::PETScUserPreconditioner::PETScUserPreconditioner "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::GlobalParameters "

    This class defines the global DOLFIN parameter database.

    C++ includes: GlobalParameters.h 
    
";

%feature("docstring")  dolfin::GlobalParameters::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::GlobalParameters::GlobalParameters "

        __init__(self) -> GlobalParameters

        Constructor. 
        
";

%feature("docstring")  dolfin::MPICommunicator "
Proxy of C++ dolfin::MPICommunicator class
";

%feature("docstring")  dolfin::MPICommunicator::MPICommunicator "
__init__(self) -> MPICommunicator
";

%feature("docstring")  dolfin::RadauQuadrature "

    Radau (Gauss-Radau) quadrature on the interval [-1,1]. The n
    quadrature points are given by the zeros of

    ( Pn-1(x) + Pn(x) ) / (1+x)

    where Pn is the n:th Legendre polynomial.

    The quadrature points are computed using Newton's method, and the
    quadrature weights are computed by solving a linear system determined
    by the condition that Radau quadrature with n points should be exact
    for polynomials of degree 2n-2.

    C++ includes: RadauQuadrature.h 
    
";

%feature("docstring")  dolfin::RadauQuadrature::RadauQuadrature "

        __init__(self, unsigned int n) -> RadauQuadrature

        Create Radau quadrature with n points. 
        
";

%feature("docstring")  dolfin::SubVector "
Proxy of C++ dolfin::SubVector class
";

%feature("docstring")  dolfin::SubVector::SubVector "
__init__(self, uint n, BlockVector bv) -> SubVector
";

%feature("docstring")  dolfin::cells "

    A CellIterator is a MeshEntityIterator of topological codimension 0.

    C++ includes: Cell.h 
    
";

%feature("docstring")  dolfin::cells::_dereference "
_dereference(self) -> Cell
";

%feature("docstring")  dolfin::cells::cells "

        CellIterator() -> cells
        CellIterator(Mesh mesh) -> cells
        __init__(self, MeshEntity entity) -> cells
        
";

%feature("docstring")  dolfin::MeshFunctionDouble "

    A MeshFunction is a function that can be evaluated at a set of mesh
    entities. A MeshFunction is discrete and is only defined at the set of
    mesh entities of a fixed topological dimension. A MeshFunction may for
    example be used to store a global numbering scheme for the entities of
    a (parallel) mesh, marking sub domains or boolean markers for mesh
    refinement.

    C++ includes: MeshFunction.h 
    
";

%feature("docstring")  dolfin::MeshFunctionDouble::dim "

        dim(self) -> uint

        Return topological dimension. 
        
";

%feature("docstring")  dolfin::MeshFunctionDouble::set_all "

        set_all(self, double value)

        Set all values to given value. 
        
";

%feature("docstring")  dolfin::MeshFunctionDouble::init "

        init(self, uint dim)
        init(self, uint dim, uint size)
        init(self, Mesh mesh, uint dim)
        init(self, Mesh mesh, uint dim, uint size)

        Initialize mesh function for given topological dimension of given
        size. 
        
";

%feature("docstring")  dolfin::MeshFunctionDouble::mesh "

        mesh(self) -> Mesh

        Return mesh associated with mesh function. 
        
";

%feature("docstring")  dolfin::MeshFunctionDouble::values "

        values(self) -> PyObject

        Return array of values. 
        
";

%feature("docstring")  dolfin::MeshFunctionDouble::MeshFunctionDouble "

        __init__(self) -> MeshFunctionDouble
        __init__(self, Mesh mesh) -> MeshFunctionDouble
        __init__(self, Mesh mesh, uint dim) -> MeshFunctionDouble
        __init__(self, Mesh mesh, uint dim, double value) -> MeshFunctionDouble
        __init__(self, Mesh mesh, string filename) -> MeshFunctionDouble
        __init__(self, MeshFunctionDouble f) -> MeshFunctionDouble

        Copy constructor. 
        
";

%feature("docstring")  dolfin::MeshFunctionDouble::size "

        size(self) -> uint

        Return size (number of entities). 
        
";

%feature("docstring")  dolfin::Scalar "

    This class represents a real-valued scalar quantity and implements the
    GenericTensor interface for scalars.

    C++ includes: Scalar.h 
    
";

%feature("docstring")  dolfin::Scalar::getval "

        getval(self) -> double

        Get value. 
        
";

%feature("docstring")  dolfin::Scalar::copy "

        copy(self) -> Scalar

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::Scalar::Scalar "

        __init__(self) -> Scalar

        Create zero scalar. 
        
";

%feature("docstring")  dolfin::Scalar::assign "
assign(self, double value) -> Scalar
";

%feature("docstring")  dolfin::SparsityPattern "

    This class implements the GenericSparsityPattern interface. It is used
    by most linear algebra backends, except for Epetra which uses a
    special/native implementation.

    C++ includes: SparsityPattern.h 
    
";

%feature("docstring")  dolfin::SparsityPattern::off_diagonal_pattern "

        off_diagonal_pattern(self) -> std::vector<(dolfin::Set<(dolfin::uint)>)>

        Return underlying sparsity pattern (off-diagional). 
        
";

%feature("docstring")  dolfin::SparsityPattern::SparsityPattern "

        __init__(self, Type type) -> SparsityPattern

        Create empty sparsity pattern. 
        
";

%feature("docstring")  dolfin::SparsityPattern::str "

        str(self) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::SparsityPattern::diagonal_pattern "

        diagonal_pattern(self) -> std::vector<(dolfin::Set<(dolfin::uint)>)>

        Return underlying sparsity pattern (diagonal). 
        
";

%feature("docstring")  dolfin::MeshCoordinates "

    This Function represents the mesh coordinates on a given mesh.

    C++ includes: SpecialFunctions.h 
    
";

%feature("docstring")  dolfin::MeshCoordinates::MeshCoordinates "

        __init__(self, Mesh mesh) -> MeshCoordinates

        Constructor. 
        
";

%feature("docstring")  dolfin::LinearAlgebraFactory "
Proxy of C++ dolfin::LinearAlgebraFactory class
";

%feature("docstring")  dolfin::LinearAlgebraFactory::create_matrix "

        create_matrix(self) -> GenericMatrix

        Create empty matrix. 
        
";

%feature("docstring")  dolfin::LinearAlgebraFactory::create_vector "

        create_vector(self) -> GenericVector

        Create empty vector (global). 
        
";

%feature("docstring")  dolfin::LinearAlgebraFactory::create_krylov_solver "

        create_krylov_solver(self, string method, string pc) -> GenericLinearSolver

        Create Krylov solver. 
        
";

%feature("docstring")  dolfin::LinearAlgebraFactory::create_pattern "

        create_pattern(self) -> GenericSparsityPattern

        Create empty sparsity pattern (returning zero if not used/needed). 
        
";

%feature("docstring")  dolfin::LinearAlgebraFactory::create_lu_solver "

        create_lu_solver(self) -> GenericLinearSolver

        Create LU solver. 
        
";

%feature("docstring")  dolfin::LinearAlgebraFactory::create_local_vector "

        create_local_vector(self) -> GenericVector

        Create empty vector (local). 
        
";

%feature("docstring")  dolfin::LinearAlgebraFactory::LinearAlgebraFactory "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::MeshFunction "
Missing docstring
";

%feature("docstring")  dolfin::Event "

    A event is a string message which is displayed only a limited number
    of times.

    Example of usage:

    Event event("System is stiff, damping is needed."); while () { ...
    if ( ... ) { event(); ... } }

    C++ includes: Event.h 
    
";

%feature("docstring")  dolfin::Event::count "

        count(self) -> unsigned int

        Display count. 
        
";

%feature("docstring")  dolfin::Event::maxcount "

        maxcount(self) -> unsigned int

        Maximum display count. 
        
";

%feature("docstring")  dolfin::Event::Event "

        __init__(self, string msg, unsigned int maxcount = 1) -> Event
        __init__(self, string msg) -> Event

        Constructor. 
        
";

%feature("docstring")  dolfin::MeshFunctionUInt "

    A MeshFunction is a function that can be evaluated at a set of mesh
    entities. A MeshFunction is discrete and is only defined at the set of
    mesh entities of a fixed topological dimension. A MeshFunction may for
    example be used to store a global numbering scheme for the entities of
    a (parallel) mesh, marking sub domains or boolean markers for mesh
    refinement.

    C++ includes: MeshFunction.h 
    
";

%feature("docstring")  dolfin::MeshFunctionUInt::dim "

        dim(self) -> uint

        Return topological dimension. 
        
";

%feature("docstring")  dolfin::MeshFunctionUInt::set_all "

        set_all(self, unsigned int value)

        Set all values to given value. 
        
";

%feature("docstring")  dolfin::MeshFunctionUInt::init "

        init(self, uint dim)
        init(self, uint dim, uint size)
        init(self, Mesh mesh, uint dim)
        init(self, Mesh mesh, uint dim, uint size)

        Initialize mesh function for given topological dimension of given
        size. 
        
";

%feature("docstring")  dolfin::MeshFunctionUInt::mesh "

        mesh(self) -> Mesh

        Return mesh associated with mesh function. 
        
";

%feature("docstring")  dolfin::MeshFunctionUInt::values "

        values(self) -> PyObject

        Return array of values. 
        
";

%feature("docstring")  dolfin::MeshFunctionUInt::MeshFunctionUInt "

        __init__(self) -> MeshFunctionUInt
        __init__(self, Mesh mesh) -> MeshFunctionUInt
        __init__(self, Mesh mesh, uint dim) -> MeshFunctionUInt
        __init__(self, Mesh mesh, uint dim, unsigned int value) -> MeshFunctionUInt
        __init__(self, Mesh mesh, string filename) -> MeshFunctionUInt
        __init__(self, MeshFunctionUInt f) -> MeshFunctionUInt

        Copy constructor. 
        
";

%feature("docstring")  dolfin::MeshFunctionUInt::size "

        size(self) -> uint

        Return size (number of entities). 
        
";

%feature("docstring")  dolfin::facets "

    A FacetIterator is a MeshEntityIterator of topological codimension 1.

    C++ includes: Facet.h 
    
";

%feature("docstring")  dolfin::facets::_dereference "
_dereference(self) -> Facet
";

%feature("docstring")  dolfin::facets::facets "

        FacetIterator(Mesh mesh) -> facets
        __init__(self, MeshEntity entity) -> facets
        
";

%feature("docstring")  dolfin::PETScBaseMatrix "
Proxy of C++ dolfin::PETScBaseMatrix class
";

%feature("docstring")  dolfin::PETScBaseMatrix::resize "
resize(self, uint m, uint n)
";

%feature("docstring")  dolfin::PETScBaseMatrix::mat "
mat(self) -> boost::shared_ptr<(Mat)>
";

%feature("docstring")  dolfin::PETScBaseMatrix::PETScBaseMatrix "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::PETScBaseMatrix::size "
size(self, uint dim) -> uint
";

%feature("docstring")  dolfin::DoubleArray "

    This class provides a simple wrapper for a pointer to an array. A
    purpose of this class is to enable the simple and safe exchange of
    data between C++ and Python.

    C++ includes: Array.h 
    
";

%feature("docstring")  dolfin::DoubleArray::min "

        min(self) -> double

        Return minimum value of array. 
        
";

%feature("docstring")  dolfin::DoubleArray::zero_eps "

        zero_eps(self, double eps = 3.0e-16)
        zero_eps(self)
        
";

%feature("docstring")  dolfin::DoubleArray::update "

        update(self, uint N, double _x)

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::DoubleArray::zero "

        zero(self)

        Zero array. 
        
";

%feature("docstring")  dolfin::DoubleArray::resize "

        resize(self, uint N)

        Resize array to size N. If size changes, contents will be destroyed.

        
";

%feature("docstring")  dolfin::DoubleArray::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). Note that the
        Array class is not a subclass of Variable (for efficiency) which means
        that one needs to call str() directly instead of using the info()
        function on Array objects. 
        
";

%feature("docstring")  dolfin::DoubleArray::max "

        max(self) -> double

        Return maximum value of array. 
        
";

%feature("docstring")  dolfin::DoubleArray::array "
array(self) -> PyObject
";

%feature("docstring")  dolfin::DoubleArray::data "

        data(self) -> boost::shared_array<(double)>
        data(self) -> boost::shared_array<(double)>

        Return pointer to data (non-const version). 
        
";

%feature("docstring")  dolfin::DoubleArray::DoubleArray "

        __init__(self) -> DoubleArray
        __init__(self, uint N) -> DoubleArray
        __init__(self, uint N) -> DoubleArray

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::DoubleArray::size "

        size(self) -> uint

        Return size of array. 
        
";

%feature("docstring")  dolfin::GenericVector "

    This class defines a common interface for vectors.

    C++ includes: GenericVector.h 
    
";

%feature("docstring")  dolfin::GenericVector::gather "

        gather(self, GenericVector x, dolfin::Array<(dolfin::uint)> indices)

        Gather entries into local vector x. 
        
";

%feature("docstring")  dolfin::GenericVector::_scale "
_scale(self, double a)
";

%feature("docstring")  dolfin::GenericVector::set_local "

        set_local(self, DoubleArray values)

        Set all values on local process. 
        
";

%feature("docstring")  dolfin::GenericVector::local_range "

        local_range(self) -> std::pair<(dolfin::uint,dolfin::uint)>

        Return local ownership range of a vector. 
        
";

%feature("docstring")  dolfin::GenericVector::array "
Return a numpy array representation of Vector
";

%feature("docstring")  dolfin::GenericVector::GenericVector "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::GenericVector::size "

        size(self, uint dim) -> uint
        size(self) -> uint

        Return global size of vector. 
        
";

%feature("docstring")  dolfin::GenericVector::_vec_mul "
_vec_mul(self, GenericVector other)
";

%feature("docstring")  dolfin::GenericVector::sum "

        sum(self) -> double
        sum(self, dolfin::Array<(dolfin::uint)> rows) -> double

        Return sum of selected rows in vector. Repeated entries only summed
        once. 
        
";

%feature("docstring")  dolfin::GenericVector::add "

        add(self, double block, uint m)

        Add block of values. 
        
";

%feature("docstring")  dolfin::GenericVector::local_size "

        local_size(self) -> uint

        Return local size of vector. 
        
";

%feature("docstring")  dolfin::GenericVector::axpy "

        axpy(self, double a, GenericVector x)

        Add multiple of given vector (AXPY operation). 
        
";

%feature("docstring")  dolfin::GenericVector::inner "

        inner(self, GenericVector x) -> double

        Return inner product with given vector. 
        
";

%feature("docstring")  dolfin::GenericVector::max "

        max(self) -> double

        Return maximum value of vector. 
        
";

%feature("docstring")  dolfin::GenericVector::_assign "

        _assign(self, GenericVector x) -> GenericVector
        _assign(self, double a) -> GenericVector
        
";

%feature("docstring")  dolfin::GenericVector::copy "

        copy(self) -> GenericVector

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::GenericVector::resize "

        resize(self, uint rank, uint dims)
        resize(self, uint N)

        Resize vector to size N. 
        
";

%feature("docstring")  dolfin::GenericVector::min "

        min(self) -> double

        Return minimum value of vector. 
        
";

%feature("docstring")  dolfin::GenericVector::get_local "

        get_local(self, double block, uint m)
        get_local(self, DoubleArray values)

        Get all values on local process. 
        
";

%feature("docstring")  dolfin::GenericVector::norm "

        norm(self, string norm_type) -> double

        Return norm of vector. 
        
";

%feature("docstring")  dolfin::GenericVector::add_local "

        add_local(self, DoubleArray values)

        Add values to each entry on local process. 
        
";

%feature("docstring")  dolfin::UnitInterval "

    Interval mesh of the 1D unit line (0,1). Given the number of cells
    (nx) in the axial direction, the total number of intervals will be nx
    and the total number of vertices will be (nx + 1).

    C++ includes: UnitInterval.h 
    
";

%feature("docstring")  dolfin::UnitInterval::UnitInterval "
__init__(self, uint nx) -> UnitInterval
";

%feature("docstring")  dolfin::Parameters "

    This class stores a set of parameters. Each parameter is identified by
    a unique string (the key) and a value of some given value type.
    Parameter sets can be nested at arbitrary depths.

    A parameter may be either int, double, string or boolean valued.

    Parameters may be added as follows:

    Parameters p("my_parameters"); p.add("relative_tolerance", 1e-15);
    p.add("absolute_tolerance", 1e-15); p.add("gmres_restart", 30);
    p.add("monitor_convergence", false);

    Parameters may be changed as follows:

    p("gmres_restart") = 50;

    Parameter values may be retrieved as follows:

    int gmres_restart = p("gmres_restart");

    Parameter sets may be nested as follows:

    Parameters q("nested_parameters"); p.add(q);

    Nested parameters may then be accessed by

    p["nested_parameters"]("...")

    Parameters may be nested at arbitrary depths.

    Parameters may be parsed from the command-line as follows:

    p.parse(argc, argv);

    Note: spaces in parameter keys are not allowed (to simplify usage from
    command-line).

    C++ includes: Parameters.h 
    
";

%feature("docstring")  dolfin::Parameters::rename "

        rename(self, string key)

        Rename parameter set. 
        
";

%feature("docstring")  dolfin::Parameters::_get_parameter_set_keys "

        _get_parameter_set_keys(self)

        Return a vector of parameter set keys. 
        
";

%feature("docstring")  dolfin::Parameters::_add "

        _add(self, string key, int value)
        _add(self, string key, int value, int min_value, int max_value)
        _add(self, string key, double value)
        _add(self, string key, double value, double min_value, double max_value)
        _add(self, string key, real value)
        _add(self, string key, real value, real min_value, real max_value)
        _add(self, string key, string value)
        _add(self, string key, char value)
        _add(self, string key, string value, std::set<(std::string)> range)
        _add(self, string key, char value, std::set<(std::string)> range)
        _add(self, Parameters parameters)

        Add nested parameter set. 
        
";

%feature("docstring")  dolfin::Parameters::clear "

        clear(self)

        Clear parameter set. 
        
";

%feature("docstring")  dolfin::Parameters::has_key "

        has_key(self, string key) -> bool

        Check if parameter set has given key. 
        
";

%feature("docstring")  dolfin::Parameters::_get_parameter "

        _get_parameter(self, string key) -> ParameterValue
        _get_parameter(self, string key) -> ParameterValue
        
";

%feature("docstring")  dolfin::Parameters::Parameters "

        __init__(self, string key = "parameters") -> Parameters
        __init__(self) -> Parameters
        __init__(self, Parameters parameters) -> Parameters

        Copy constructor. 
        
";

%feature("docstring")  dolfin::Parameters::set_range "
Set the range for the given parameter
";

%feature("docstring")  dolfin::Parameters::_get_parameter_set "

        _get_parameter_set(self, string key) -> Parameters
        _get_parameter_set(self, string key) -> Parameters
        
";

%feature("docstring")  dolfin::Parameters::add "
Add a parameter to the parameter set
";

%feature("docstring")  dolfin::Parameters::itervalues "
Returns an iterator to the parameter values
";

%feature("docstring")  dolfin::Parameters::iterdata "
Returns an iterator of a tuple of a parameter key together with its value
";

%feature("docstring")  dolfin::Parameters::get "
Return all data available for a certain parameter

        The data is returned in a tuple:
        value, range, access_count, change_count = parameters.get('name')
        
";

%feature("docstring")  dolfin::Parameters::keys "
Returns a list of the parameter keys
";

%feature("docstring")  dolfin::Parameters::update "
A recursive update that handles parameter subsets correctly.
";

%feature("docstring")  dolfin::Parameters::iteritems "
Returns an iterator over the (key, value) items of the Parameters
";

%feature("docstring")  dolfin::Parameters::copy "
Return a copy of it self
";

%feature("docstring")  dolfin::Parameters::iterkeys "
Returns an iterator for the parameter keys
";

%feature("docstring")  dolfin::Parameters::parse "
Parse command line arguments
";

%feature("docstring")  dolfin::Parameters::option_string "
Return an option string representation of the Parameters
";

%feature("docstring")  dolfin::Parameters::name "

        name(self) -> string

        Return name for parameter set. 
        
";

%feature("docstring")  dolfin::Parameters::items "
Missing docstring
";

%feature("docstring")  dolfin::Parameters::_add_bool "

        _add_bool(self, string key, bool value)

        Add nested parameter set. 
        
";

%feature("docstring")  dolfin::Parameters::values "
Returns a list of the parameter values
";

%feature("docstring")  dolfin::Parameters::to_dict "
Convert the Parameters to a dict
";

%feature("docstring")  dolfin::Parameters::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::Parameters::_parse "
_parse(self, PyObject op)
";

%feature("docstring")  dolfin::Parameters::assign "
assign(self, Parameters parameters) -> Parameters
";

%feature("docstring")  dolfin::Parameters::_get_parameter_keys "

        _get_parameter_keys(self)

        Return a vector of parameter keys. 
        
";

%feature("docstring")  dolfin::BoundaryMesh "

    A BoundaryMesh is a mesh over the boundary of some given mesh.

    C++ includes: BoundaryMesh.h 
    
";

%feature("docstring")  dolfin::BoundaryMesh::BoundaryMesh "

        __init__(self) -> BoundaryMesh
        __init__(self, Mesh mesh) -> BoundaryMesh

        Create (interior) boundary mesh from given mesh. 
        
";

%feature("docstring")  dolfin::BoundaryMesh::init_interior_boundary "

        init_interior_boundary(self, Mesh mesh)

        Initialize interior boundary of given mesh. 
        
";

%feature("docstring")  dolfin::BoundaryMesh::init_exterior_boundary "

        init_exterior_boundary(self, Mesh mesh)

        Initialize exterior boundary of given mesh. 
        
";

%feature("docstring")  dolfin::MeshPartitioning "

    This class partitions and distributes a mesh based on partitioned
    local mesh data. Note that the local mesh data will also be
    repartitioned and redistributed during the computation of the mesh
    partitioning.

    After partitioning, each process has a local mesh and set of mesh data
    that couples the meshes together.

    The following mesh data is created:

    1. "global entity indices 0" ( MeshFunction<uint>)

    This maps each local vertex to its global index.

    2. "overlap" (std::map<uint, std::vector<uint> >)

    This maps each shared vertex to a list of the processes sharing the
    vertex.

    3. "global entity indices %d" ( MeshFunction<uint>)

    After partitioning, the function number_entities() may be called to
    create global indices for all entities of a given topological
    dimension. These are stored as mesh data ( MeshFunction<uint>) named

    "global entity indices 1" "global entity indices 2" etc

    4. "num global entities" (std::vector<uint>)

    The function number_entities also records the number of global
    entities for the dimension of the numbered entities in the array named
    "num global entities". This array has size D + 1, where D is the
    topological dimension of the mesh. This array is initially created by
    the mesh and then contains only the number entities of dimension 0
    (vertices) and dimension D (cells).

    C++ includes: MeshPartitioning.h 
    
";

%feature("docstring")  dolfin::MeshPartitioning::number_entities "
number_entities(Mesh mesh, uint d)
";

%feature("docstring")  dolfin::MeshPartitioning::partition "

        partition(Mesh mesh)
        partition(Mesh mesh, LocalMeshData data)
        
";

%feature("docstring")  dolfin::MeshPartitioning::MeshPartitioning "
__init__(self) -> MeshPartitioning
";

%feature("docstring")  dolfin::UnitSquare "

    Triangular mesh of the 2D unit square (0,1) x (0,1). Given the number
    of cells (nx, ny) in each direction, the total number of triangles
    will be 2*nx*ny and the total number of vertices will be (nx + 1)*(ny
    + 1). std::string diagonal ("left", "right", "right/left" or
    "crossed") indicates the direction of the diagonals.

    C++ includes: UnitSquare.h 
    
";

%feature("docstring")  dolfin::UnitSquare::UnitSquare "

        __init__(self, uint nx, uint ny, string diagonal = "right") -> UnitSquare
        __init__(self, uint nx, uint ny) -> UnitSquare
        
";

%feature("docstring")  dolfin::DefaultFactory "
Proxy of C++ dolfin::DefaultFactory class
";

%feature("docstring")  dolfin::DefaultFactory::DefaultFactory "

        __init__(self) -> DefaultFactory

        Constructor. 
        
";

%feature("docstring")  dolfin::Quadrature "
Proxy of C++ dolfin::Quadrature class
";

%feature("docstring")  dolfin::Quadrature::weight "

        weight(self, unsigned int i) -> real

        Return quadrature weight. 
        
";

%feature("docstring")  dolfin::Quadrature::point "

        point(self, unsigned int i) -> real

        Return quadrature point. 
        
";

%feature("docstring")  dolfin::Quadrature::measure "

        measure(self) -> real

        Return sum of weights (length, area, volume). 
        
";

%feature("docstring")  dolfin::Quadrature::Quadrature "

        __init__(self, unsigned int n) -> Quadrature

        Constructor. 
        
";

%feature("docstring")  dolfin::Quadrature::size "

        size(self) -> int

        Return number of quadrature points. 
        
";

%feature("docstring")  dolfin::NewtonSolver "

    This class defines a Newton solver for equations of the form F(u) = 0.

    C++ includes: NewtonSolver.h 
    
";

%feature("docstring")  dolfin::NewtonSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::NewtonSolver::iteration "

        iteration(self) -> uint

        Return Newton iteration number. 
        
";

%feature("docstring")  dolfin::NewtonSolver::residual "

        residual(self) -> double

        Return current residual. 
        
";

%feature("docstring")  dolfin::NewtonSolver::relative_residual "

        relative_residual(self) -> double

        Return current relative residual. 
        
";

%feature("docstring")  dolfin::NewtonSolver::solve "

        solve(self, NonlinearProblem nonlinear_function, GenericVector x) -> std::pair<(dolfin::uint,bool)>

        Solve abstract nonlinear problem F(x) = 0 for given vector F and
        Jacobian dF/dx 
        
";

%feature("docstring")  dolfin::NewtonSolver::linear_solver "

        linear_solver(self) -> GenericLinearSolver

        Return the linear solver. 
        
";

%feature("docstring")  dolfin::NewtonSolver::NewtonSolver "

        __init__(self, string solver_type = "lu", string pc_type = "default") -> NewtonSolver
        __init__(self, string solver_type = "lu") -> NewtonSolver
        __init__(self) -> NewtonSolver
        __init__(self, GenericLinearSolver solver, LinearAlgebraFactory factory) -> NewtonSolver

        Create nonlinear solver using provided linear solver and linear
        algebra backend determined by factory 
        
";

%feature("docstring")  dolfin::uBLASSparseFactory "
Proxy of C++ dolfin::uBLASFactory<(dolfin::ublas_sparse_matrix)> class
";

%feature("docstring")  dolfin::uBLASSparseFactory::create_matrix "

        create_matrix(self) -> uBLASSparseMatrix

        Create empty matrix. 
        
";

%feature("docstring")  dolfin::uBLASSparseFactory::create_vector "

        create_vector(self) -> uBLASVector

        Create empty vector. 
        
";

%feature("docstring")  dolfin::uBLASSparseFactory::instance "
instance() -> uBLASSparseFactory
";

%feature("docstring")  dolfin::uBLASSparseFactory::create_lu_solver "

        create_lu_solver(self) -> UmfpackLUSolver

        Create LU solver. 
        
";

%feature("docstring")  dolfin::uBLASSparseFactory::create_local_vector "

        create_local_vector(self) -> uBLASVector

        Create empty vector (local). 
        
";

%feature("docstring")  dolfin::uBLASSparseFactory::create_pattern "

        create_pattern(self) -> SparsityPattern

        Create empty sparsity pattern. 
        
";

%feature("docstring")  dolfin::uBLASSparseFactory::uBLASSparseFactory "
No constructor defined
";

%feature("docstring")  dolfin::Form "

    Base class for UFC code generated by FFC for DOLFIN with option -l.

    C++ includes: Form.h 
    
";

%feature("docstring")  dolfin::Form::coefficient "

        coefficient(self, uint i) -> GenericFunction
        coefficient(self, string name) -> GenericFunction

        Return coefficient with given name. 
        
";

%feature("docstring")  dolfin::Form::set_coefficient "

        set_coefficient(self, uint i, GenericFunction coefficient)
        set_coefficient(self, uint i, __dummy_23__ coefficient)
        set_coefficient(self, string name, GenericFunction coefficient)
        set_coefficient(self, string name, __dummy_23__ coefficient)

        Set coefficient with given name (shared pointer version). 
        
";

%feature("docstring")  dolfin::Form::set_coefficients "

        set_coefficients(self, std::map<(std::string,p.q(const).dolfin::GenericFunction)> coefficients)
        set_coefficients(self, std::map<(std::string,boost::shared_ptr<(q(const).dolfin::GenericFunction)>)> coefficients)

        Set all coefficients in given map, possibly a subset (shared pointer
        version). 
        
";

%feature("docstring")  dolfin::Form::rank "

        rank(self) -> uint

        Return rank of form (bilinear form = 2, linear form = 1, functional =
        0, etc). 
        
";

%feature("docstring")  dolfin::Form::coefficient_number "

        coefficient_number(self, string name) -> uint

        Return the number of the coefficient with this name. 
        
";

%feature("docstring")  dolfin::Form::mesh "

        mesh(self) -> Mesh

        Return mesh. 
        
";

%feature("docstring")  dolfin::Form::ufc_form "

        ufc_form(self) -> form

        Return UFC form. 
        
";

%feature("docstring")  dolfin::Form::coefficient_name "

        coefficient_name(self, uint i) -> string

        Return the name of the coefficient with this number. 
        
";

%feature("docstring")  dolfin::Form::check "

        check(self)

        Check function spaces and coefficients. 
        
";

%feature("docstring")  dolfin::Form::coefficients "

        coefficients(self) -> std::vector<(p.q(const).dolfin::GenericFunction)>

        Return all coefficients. 
        
";

%feature("docstring")  dolfin::Form::set_mesh "

        set_mesh(self, Mesh mesh)
        set_mesh(self, __dummy_37__ mesh)

        Set mesh, necessary for functionals when there are no function spaces.

        
";

%feature("docstring")  dolfin::Form::function_space "

        function_space(self, uint i) -> __dummy_19__

        Return function space for given argument. 
        
";

%feature("docstring")  dolfin::Form::Form "

        __init__(self, uint rank, uint num_coefficients) -> Form
        __init__(self, form ufc_form, std::vector<(p.q(const).dolfin::FunctionSpace)> function_spaces, 
            std::vector<(p.q(const).dolfin::GenericFunction)> coefficients) -> Form

        Create form (constructor used from Python interface). 
        
";

%feature("docstring")  dolfin::Form::function_spaces "

        function_spaces(self) -> std::vector<(boost::shared_ptr<(q(const).dolfin::FunctionSpace)>)>

        Return function spaces for arguments. 
        
";

%feature("docstring")  dolfin::Form::num_coefficients "

        num_coefficients(self) -> uint

        Return number of coefficients. 
        
";

%feature("docstring")  dolfin::Matrix "

    This class provides the default DOLFIN matrix class, based on the
    default DOLFIN linear algebra backend.

    C++ includes: Matrix.h 
    
";

%feature("docstring")  dolfin::Matrix::zero "

        zero(self)
        zero(self, uint m)

        Set given rows to zero. 
        
";

%feature("docstring")  dolfin::Matrix::copy "

        copy(self) -> Matrix

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::Matrix::assign "

        assign(self, GenericMatrix A) -> GenericMatrix
        assign(self, Matrix A) -> Matrix
        
";

%feature("docstring")  dolfin::Matrix::Matrix "

        __init__(self) -> Matrix
        __init__(self, uint M, uint N) -> Matrix
        __init__(self, Matrix A) -> Matrix

        Copy constructor. 
        
";

%feature("docstring")  dolfin::IntersectionOperator "
Proxy of C++ dolfin::IntersectionOperator class
";

%feature("docstring")  dolfin::IntersectionOperator::closest_cell "

        closest_cell(self, Point point) -> uint

        Computes the index of the cell inside the mesh which are closest to
        the point query. 
        
";

%feature("docstring")  dolfin::IntersectionOperator::clear "

        clear(self)

        Clears search structure. Should be used if the mesh has changed. 
        
";

%feature("docstring")  dolfin::IntersectionOperator::closest_point_and_cell "

        closest_point_and_cell(self, Point point) -> std::pair<(dolfin::Point,dolfin::uint)>

        Computes the point inside the mesh and the corresponding cell index
        which are closest to the point query. 
        
";

%feature("docstring")  dolfin::IntersectionOperator::mesh "
mesh(self) -> Mesh
";

%feature("docstring")  dolfin::IntersectionOperator::all_intersected_entities "

        all_intersected_entities(self, Point point)
        all_intersected_entities(self, std::vector<(dolfin::Point)> points)
        all_intersected_entities(self, MeshEntity entity, std::vector<(dolfin::uint)> ids_result)
        all_intersected_entities(self, std::vector<(dolfin::MeshEntity)> entities)
        all_intersected_entities(self, Mesh another_mesh)

        Compute all id of all cells which are intersects by the given mesh
        another_mesh;

        Parameters:
        -----------

        ids_result:  The ids of the intersected entities are saved in a set
        for efficienty reasons, to avoid to sort out duplicates later on. 
        
";

%feature("docstring")  dolfin::IntersectionOperator::closest_point "

        closest_point(self, Point point) -> Point

        Computes the point inside the mesh which are closest to the point
        query. 
        
";

%feature("docstring")  dolfin::IntersectionOperator::any_intersected_entity "

        any_intersected_entity(self, Point point) -> int

        Computes only the first id of the entity, which contains the point.
        Returns -1 if no cell is intersected. 
        
";

%feature("docstring")  dolfin::IntersectionOperator::reset_kernel "

        reset_kernel(self, string kernel_type = "SimpleCartesian")
        reset_kernel(self)

        Rebuilds the underlying search structure from scratch and uses the
        kernel kernel_type underlying CGAL Geometry kernel. 
        
";

%feature("docstring")  dolfin::IntersectionOperator::IntersectionOperator "

        __init__(self, Mesh _mesh, string kernel_type = "SimpleCartesian") -> IntersectionOperator
        __init__(self, Mesh _mesh) -> IntersectionOperator
        __init__(self, __dummy_37__ _mesh, string kernel_type = "SimpleCartesian") -> IntersectionOperator
        __init__(self, __dummy_37__ _mesh) -> IntersectionOperator
        
";

%feature("docstring")  dolfin::CholmodCholeskySolver "

    This class implements the direct solution (Cholesky factorization) of
    linear systems of the form Ax = b. Sparse matrices are solved using
    CHOLMODhttp://www.cise.ufl.edu/research/sparse/cholmod/ if installed.

    C++ includes: CholmodCholeskySolver.h 
    
";

%feature("docstring")  dolfin::CholmodCholeskySolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::CholmodCholeskySolver::factorize "

        factorize(self, GenericMatrix A) -> uint

        Cholesky-factor sparse matrix A if CHOLMOD is installed. 
        
";

%feature("docstring")  dolfin::CholmodCholeskySolver::factorized_solve "

        factorized_solve(self, GenericVector x, GenericVector b) -> uint

        Solve factorized system (CHOLMOD). 
        
";

%feature("docstring")  dolfin::CholmodCholeskySolver::CholmodCholeskySolver "

        __init__(self) -> CholmodCholeskySolver
        __init__(self, GenericMatrix A) -> CholmodCholeskySolver
        __init__(self, boost::shared_ptr<(q(const).dolfin::GenericMatrix)> A) -> CholmodCholeskySolver

        Constructor. 
        
";

%feature("docstring")  dolfin::PETScObject "

    This class calls SubSystemsManger to initialise PETSc.

    All PETSc objects must be derived from this class.

    C++ includes: PETScObject.h 
    
";

%feature("docstring")  dolfin::PETScObject::PETScObject "
__init__(self) -> PETScObject
";

%feature("docstring")  dolfin::LUSolver "
Proxy of C++ dolfin::LUSolver class
";

%feature("docstring")  dolfin::LUSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::LUSolver::solve "

        solve(self, GenericVector x, GenericVector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint

        Solve linear system. 
        
";

%feature("docstring")  dolfin::LUSolver::LUSolver "

        __init__(self, string type = "lu") -> LUSolver
        __init__(self) -> LUSolver
        __init__(self, GenericMatrix A, string type = "lu") -> LUSolver
        __init__(self, GenericMatrix A) -> LUSolver

        Constructor. 
        
";

%feature("docstring")  dolfin::CellType "

    This class provides a common interface for different cell types. Each
    cell type implements mesh functionality that is specific to a certain
    type of cell.

    C++ includes: CellType.h 
    
";

%feature("docstring")  dolfin::CellType::diameter "

        diameter(self, MeshEntity entity) -> double

        Compute diameter of mesh entity. 
        
";

%feature("docstring")  dolfin::CellType::ordered "

        ordered(self, Cell cell, MeshFunctionUInt global_vertex_indices) -> bool

        Check if entities are ordered. 
        
";

%feature("docstring")  dolfin::CellType::orientation "

        orientation(self, Cell cell) -> uint

        Return orientation of the cell. 
        
";

%feature("docstring")  dolfin::CellType::normal "

        normal(self, Cell cell, uint facet, uint i) -> double
        normal(self, Cell cell, uint facet) -> Point

        Compute of given facet with respect to the cell. 
        
";

%feature("docstring")  dolfin::CellType::num_entities "

        num_entities(self, uint dim) -> uint

        Return number of entitites of given topological dimension. 
        
";

%feature("docstring")  dolfin::CellType::volume "

        volume(self, MeshEntity entity) -> double

        Compute (generalized) volume of mesh entity. 
        
";

%feature("docstring")  dolfin::CellType::string2type "
string2type(string type) -> Type
";

%feature("docstring")  dolfin::CellType::num_vertices "

        num_vertices(self, uint dim) -> uint

        Return number of vertices for entity of given topological dimension.

        
";

%feature("docstring")  dolfin::CellType::type2string "
type2string(Type type) -> string
";

%feature("docstring")  dolfin::CellType::dim "

        dim(self) -> uint

        Return topological dimension of cell. 
        
";

%feature("docstring")  dolfin::CellType::cell_type "

        cell_type(self) -> Type

        Return type of cell. 
        
";

%feature("docstring")  dolfin::CellType::facet_area "

        facet_area(self, Cell cell, uint facet) -> double

        Compute the area/length of given facet with respect to the cell. 
        
";

%feature("docstring")  dolfin::CellType::create "

        create(Type type) -> CellType
        create(string type) -> CellType
        
";

%feature("docstring")  dolfin::CellType::description "

        description(self, bool plural) -> string

        Return description of cell type. 
        
";

%feature("docstring")  dolfin::CellType::refine_cell "

        refine_cell(self, Cell cell, MeshEditor editor, uint current_cell)

        Refine cell uniformly. 
        
";

%feature("docstring")  dolfin::CellType::facet_type "

        facet_type(self) -> Type

        Return type of cell for facets. 
        
";

%feature("docstring")  dolfin::CellType::create_entities "

        create_entities(self, uint e, uint dim, uint v)

        Create entities e of given topological dimension from vertices v. 
        
";

%feature("docstring")  dolfin::CellType::CellType "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::CellType::order "

        order(self, Cell cell, MeshFunctionUInt global_vertex_indices)

        Order entities locally. 
        
";

%feature("docstring")  dolfin::EqualityBC "

    This class specifies the interface for setting equality boundary
    conditions for partial differential equations,

    u(x) = u(y), for all x and y on G,

    where G is subdomain of the mesh.

    The sub domain G may be specified in two different ways. Both of them
    produce a set of unknowns (dofs) with should be equal.

    The simplest approach is to specify a SubDomain object, using the
    inside() function to specify on which facets the boundary condition
    should be applied.

    Alternatively, the boundary may be specified by the boundary
    indicators included in the mesh.

    Current implementation assume that the problem is scalar, so in case
    of mixed systems (vector-valued and mixed elements) all compoments
    will be set equal.

    C++ includes: EqualityBC.h 
    
";

%feature("docstring")  dolfin::EqualityBC::init_from_mesh "
init_from_mesh(self, uint sub_domain)
";

%feature("docstring")  dolfin::EqualityBC::apply "

        apply(self, GenericMatrix A)
        apply(self, GenericVector b)
        apply(self, GenericMatrix A, GenericVector b)
        apply(self, GenericVector b, GenericVector x)
        apply(self, GenericMatrix A, GenericVector b, GenericVector x)

        Apply boundary condition to a linear system for a nonlinear problem.

        
";

%feature("docstring")  dolfin::EqualityBC::init_from_sub_domain "
init_from_sub_domain(self, SubDomain sub_domain)
";

%feature("docstring")  dolfin::EqualityBC::EqualityBC "

        __init__(self, FunctionSpace V, SubDomain sub_domain) -> EqualityBC
        __init__(self, __dummy_19__ V, SubDomain sub_domain) -> EqualityBC
        __init__(self, __dummy_19__ V, uint sub_domain) -> EqualityBC
        
";

%feature("docstring")  dolfin::BlockVector "
Proxy of C++ dolfin::BlockVector class
";

%feature("docstring")  dolfin::BlockVector::set "

        set(self, uint i, Vector v)

        Set function. 
        
";

%feature("docstring")  dolfin::BlockVector::get "

        get(self, uint i) -> Vector
        get(self, uint arg0) -> Vector
        
";

%feature("docstring")  dolfin::BlockVector::max "

        max(self) -> double

        Return maximum value of vector. 
        
";

%feature("docstring")  dolfin::BlockVector::copy "

        copy(self) -> BlockVector

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::BlockVector::axpy "

        axpy(self, double a, BlockVector x)

        Add multiple of given vector (AXPY operation). 
        
";

%feature("docstring")  dolfin::BlockVector::size "

        size(self) -> uint

        Number of vectors. 
        
";

%feature("docstring")  dolfin::BlockVector::min "

        min(self) -> double

        Return minimum value of vector. 
        
";

%feature("docstring")  dolfin::BlockVector::BlockVector "

        __init__(self, uint n_ = 0, bool owner = False) -> BlockVector
        __init__(self, uint n_ = 0) -> BlockVector
        __init__(self) -> BlockVector

        Constructor. 
        
";

%feature("docstring")  dolfin::BlockVector::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::BlockVector::norm "

        norm(self, string norm_type) -> double

        Return norm of vector. 
        
";

%feature("docstring")  dolfin::BlockVector::inner "

        inner(self, BlockVector x) -> double

        Return inner product with given vector. 
        
";

%feature("docstring")  dolfin::LinearSolver "

    This class provides a general solver for linear systems Ax = b.

    C++ includes: LinearSolver.h 
    
";

%feature("docstring")  dolfin::LinearSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::LinearSolver::solve "

        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint
        solve(self, GenericVector x, GenericVector b) -> uint

        Solve linear system Ax = b. 
        
";

%feature("docstring")  dolfin::LinearSolver::LinearSolver "

        __init__(self, string solver_type = "lu", string pc_type = "ilu") -> LinearSolver
        __init__(self, string solver_type = "lu") -> LinearSolver
        __init__(self) -> LinearSolver

        Create linear solver. 
        
";

%feature("docstring")  dolfin::BasisFunction "

    This class represents a finite element basis function. It can be used
    for computation of basis function values and derivatives.

    Evaluation of basis functions is also possible through the use of the
    functions evaluate_basis and evaluate_basis_derivatives available in
    the FiniteElement class. The BasisFunction class relies on these
    functions for evaluation but also implements the ufc::function
    interface which allows evaluate_dof to be evaluated for a basis
    function (on a possibly different element).

    C++ includes: BasisFunction.h 
    
";

%feature("docstring")  dolfin::BasisFunction::BasisFunction "

        __init__(self, uint index, FiniteElement element, cell cell) -> BasisFunction

        Create basis function with given index on element on given cell. 
        
";

%feature("docstring")  dolfin::BasisFunction::eval_derivatives "

        eval_derivatives(self, double values, double x, uint n)

        Evaluate all order n derivatives at given point. 
        
";

%feature("docstring")  dolfin::BasisFunction::eval "

        eval(self, double values, double x)

        Evaluate basis function at given point. 
        
";

%feature("docstring")  dolfin::PointPrimitive "
Proxy of C++ dolfin::PointPrimitive class
";

%feature("docstring")  dolfin::PointPrimitive::PointPrimitive "
__init__(self) -> PointPrimitive
";

%feature("docstring")  dolfin::SubMatrix "
Proxy of C++ dolfin::SubMatrix class
";

%feature("docstring")  dolfin::SubMatrix::SubMatrix "
__init__(self, uint row, uint col, BlockMatrix bm) -> SubMatrix
";

%feature("docstring")  dolfin::Timer "

    A timer can be used for timing tasks. The basic usage is

    Timer timer("Assembling over cells");

    The timer is started at construction and timing ends when the timer is
    destroyed (goes out of scope). It is also possible to start and stop a
    timer explicitly by

    timer.start(); timer.stop();

    Timings are stored globally and a summary may be printed by calling

    summary();

    C++ includes: Timer.h 
    
";

%feature("docstring")  dolfin::Timer::start "

        start(self)

        Start timer. 
        
";

%feature("docstring")  dolfin::Timer::Timer "

        __init__(self, string task) -> Timer

        Create timer. 
        
";

%feature("docstring")  dolfin::Timer::stop "

        stop(self)

        Stop timer. 
        
";

%feature("docstring")  dolfin::Timer::value "

        value(self) -> double

        Return value of timer (or time at start if not stopped). 
        
";

%feature("docstring")  dolfin::UnitCircle "

    std:string transformation ("maxn", "sumn" or "rotsumn")

    Triangular mesh of the 2D unit circle. Given the number of cells (nx,
    ny) in each direction, the total number of triangles will be 2*nx*ny
    and the total number of vertices will be (nx + 1)*(ny + 1).
    std::string diagonal ("left", "right" or "crossed") indicates
    the direction of the diagonals.

    C++ includes: UnitCircle.h 
    
";

%feature("docstring")  dolfin::UnitCircle::UnitCircle "

        __init__(self, uint nx, string diagonal = "crossed", string transformation = "rotsumn") -> UnitCircle
        __init__(self, uint nx, string diagonal = "crossed") -> UnitCircle
        __init__(self, uint nx) -> UnitCircle
        
";

%feature("docstring")  dolfin::Facet "

    A Facet is a MeshEntity of topological codimension 1.

    C++ includes: Facet.h 
    
";

%feature("docstring")  dolfin::Facet::interior "

        interior(self) -> bool

        Determine whether or not facet is an interior facet. This is
        'relative' to the given partition of the mesh if the mesh is
        distributed 
        
";

%feature("docstring")  dolfin::Facet::Facet "

        __init__(self, Mesh mesh, uint index) -> Facet

        Constructor. 
        
";

%feature("docstring")  dolfin::Facet::adjacent_cells "

        adjacent_cells(self, MeshFunctionUInt facet_orientation = None) -> std::pair<(q(const).dolfin::Cell,q(const).dolfin::Cell)>
        adjacent_cells(self) -> std::pair<(q(const).dolfin::Cell,q(const).dolfin::Cell)>

        Return adjacent cells. An optional argument that lists for each facet
        the index of the first cell may be given to specify the ordering of
        the two cells. If not specified, the ordering will depend on the
        (arbitrary) ordering of the mesh connectivity. 
        
";

%feature("docstring")  dolfin::ODECollection "

    An ODECollection represents a collection of initial value problems of
    the form

    u'(t) = f(u(t), t) on (0, T],

    u(0) = u0,

    where u(t) is a vector of length N.

    Each ODE is governed by the same equation but a separate state is
    maintained for each ODE. Using ODECollection is recommended when
    solving a large number of ODEs and the overhead of instantiating a
    large number of ODE objects should be avoided.

    C++ includes: ODECollection.h 
    
";

%feature("docstring")  dolfin::ODECollection::set_state "

        set_state(self, uint system, real u)
        set_state(self, real u)

        Set states for all ODE systems. 
        
";

%feature("docstring")  dolfin::ODECollection::solve "

        solve(self, real t0, real t1)

        Solve ODE collection on [t0, t1]. 
        
";

%feature("docstring")  dolfin::ODECollection::get_state "

        get_state(self, uint system, real u)
        get_state(self, real u)

        Get states for all ODE systems. 
        
";

%feature("docstring")  dolfin::ODECollection::update "

        update(self, real u, real t, uint system)

        Optional user-defined update, called between solves. 
        
";

%feature("docstring")  dolfin::ODECollection::ODECollection "

        __init__(self, ODE ode, uint num_systems) -> ODECollection

        Create a collection of ODEs. 
        
";

%feature("docstring")  dolfin::DynamicMeshEditor "

    This class provides an interface for dynamic editing of meshes, that
    is, when the number of vertices and cells are not known a priori. If
    the number of vertices and cells are known a priori, it is more
    efficient to use the default editor MeshEditor.

    C++ includes: DynamicMeshEditor.h 
    
";

%feature("docstring")  dolfin::DynamicMeshEditor::add_cell "

        add_cell(self, uint c, std::vector<(dolfin::uint)> v)
        add_cell(self, uint c, uint v0, uint v1)
        add_cell(self, uint c, uint v0, uint v1, uint v2)
        add_cell(self, uint c, uint v0, uint v1, uint v2, uint v3)

        Add cell (tetrahedron) with given vertices. 
        
";

%feature("docstring")  dolfin::DynamicMeshEditor::close "

        close(self, bool order = False)
        close(self)

        Close mesh, finish editing, and order entities locally. 
        
";

%feature("docstring")  dolfin::DynamicMeshEditor::open "

        open(self, Mesh mesh, Type type, uint tdim, uint gdim)
        open(self, Mesh mesh, string type, uint tdim, uint gdim)

        Open mesh of given cell type, topological and geometrical dimension.

        
";

%feature("docstring")  dolfin::DynamicMeshEditor::DynamicMeshEditor "

        __init__(self) -> DynamicMeshEditor

        Constructor. 
        
";

%feature("docstring")  dolfin::DynamicMeshEditor::add_vertex "

        add_vertex(self, uint v, Point p)
        add_vertex(self, uint v, double x)
        add_vertex(self, uint v, double x, double y)
        add_vertex(self, uint v, double x, double y, double z)

        Add vertex v at given coordinate (x, y, z). 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern "

    Base class (interface) for generic tensor sparsity patterns.
    Currently, this interface is mostly limited to matrices.

    C++ includes: GenericSparsityPattern.h 
    
";

%feature("docstring")  dolfin::GenericSparsityPattern::num_nonzeros_diagonal "

        num_nonzeros_diagonal(self, uint num_nonzeros)

        Fill array with number of nonzeros for diagonal block in local_range
        for dimension 0. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::insert "

        insert(self, uint num_rows, uint rows)

        Insert non-zero entries. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::num_nonzeros "

        num_nonzeros(self) -> uint

        Return total number of nonzeros in local_range for dimension 0. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::init "

        init(self, uint rank, uint dims)

        Initialize sparsity pattern for a generic tensor. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::rank "

        rank(self) -> uint

        Return rank. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::local_range "

        local_range(self, uint dim) -> std::pair<(dolfin::uint,dolfin::uint)>

        Return local range for dimension dim. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::num_nonzeros_off_diagonal "

        num_nonzeros_off_diagonal(self, uint num_nonzeros)

        Fill array with number of nonzeros for off-diagonal block in
        local_range for dimension 0. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::apply "

        apply(self)

        Finalize sparsity pattern. 
        
";

%feature("docstring")  dolfin::GenericSparsityPattern::GenericSparsityPattern "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::GenericSparsityPattern::size "

        size(self, uint i) -> uint

        Return global size for dimension i. 
        
";

%feature("docstring")  dolfin::VariationalProblem "

    This class represents a (system of) partial differential equation(s)
    in variational form: Find u in V such that

    F_u(v) = 0 for all v in V'.

    The variational problem is defined in terms of a bilinear form a(v, u)
    and a linear for L(v).

    For a linear variational problem, F_u(v) = a(v, u) - L(v), the forms
    should correspond to the canonical formulation

    a(v, u) = L(v) for all v in V'.

    For a nonlinear variational problem, the forms should be given by

    a(v, u) = F_u'(v) u = F_u'(v, u), L(v) = F(v),

    that is, a(v, u) should be the Frechet derivative of F_u with respect
    to u, and L = F.

    Parameters:

    "linear solvers": "direct" or "iterative" (default: "direct")
    "symmetric": true or false (default: false)

    C++ includes: VariationalProblem.h 
    
";

%feature("docstring")  dolfin::VariationalProblem::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::VariationalProblem::newton_solver "

        newton_solver(self) -> NewtonSolver

        Return Newton solver (only useful when solving a nonlinear problem).

        
";

%feature("docstring")  dolfin::VariationalProblem::solve "

        solve(self, Function u)
        solve(self, Function u0, Function u1)
        solve(self, Function u0, Function u1, Function u2)

        Solve variational problem and extract sub functions. 
        
";

%feature("docstring")  dolfin::VariationalProblem::update "

        update(self, GenericVector x)

        Optional callback called before calls to F() and J(). 
        
";

%feature("docstring")  dolfin::VariationalProblem::VariationalProblem "

        __init__(self, Form a, Form L, bool nonlinear = False) -> VariationalProblem
        __init__(self, Form a, Form L) -> VariationalProblem
        __init__(self, Form a, Form L, BoundaryCondition bc, bool nonlinear = False) -> VariationalProblem
        __init__(self, Form a, Form L, BoundaryCondition bc) -> VariationalProblem
        __init__(self, Form a, Form L, std::vector<(p.q(const).dolfin::BoundaryCondition)> bcs, 
            bool nonlinear = False) -> VariationalProblem
        __init__(self, Form a, Form L, std::vector<(p.q(const).dolfin::BoundaryCondition)> bcs) -> VariationalProblem
        __init__(self, Form a, Form L, std::vector<(p.q(const).dolfin::BoundaryCondition)> bcs, 
            MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains, 
            bool nonlinear = False) -> VariationalProblem
        __init__(self, Form a, Form L, std::vector<(p.q(const).dolfin::BoundaryCondition)> bcs, 
            MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains) -> VariationalProblem

        Define variational problem with a list of Dirichlet boundary
        conditions and subdomains 
        
";

%feature("docstring")  dolfin::SystemAssembler "

    This class provides implements an assembler for systems of the form Ax
    = b. It differs from the default DOLFIN assembler in that it assembles
    both A and b and the same time (leading to better performance) and in
    that it applies boundary conditions at the time of assembly.

    C++ includes: SystemAssembler.h 
    
";

%feature("docstring")  dolfin::SystemAssembler::assemble "

        assemble(GenericMatrix A, GenericVector b, Form a, Form L, bool reset_sparsity = True, 
            bool add_values = False)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, bool reset_sparsity = True)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc, 
            bool reset_sparsity = True, bool add_values = True)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc, 
            bool reset_sparsity = True)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
            bool reset_sparsity = True, 
            bool add_values = False)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
            bool reset_sparsity = True)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
            MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains, 
            GenericVector x0, 
            bool reset_sparsity = True, bool add_values = False)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
            MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains, 
            GenericVector x0, 
            bool reset_sparsity = True)
        assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
            MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains, 
            GenericVector x0)
        
";

%feature("docstring")  dolfin::SystemAssembler::SystemAssembler "
__init__(self) -> SystemAssembler
";

%feature("docstring")  dolfin::NonlinearProblem "

    This is a base class for nonlinear problems which can return the
    nonlinear function F(u) and its Jacobian J = dF(u)/du.

    C++ includes: NonlinearProblem.h 
    
";

%feature("docstring")  dolfin::NonlinearProblem::form "

        form(self, GenericMatrix A, GenericVector b, GenericVector x)

        Function called by Newton solver before requesting F or J. This can be
        used to compute F and J together 
        
";

%feature("docstring")  dolfin::NonlinearProblem::F "

        F(self, GenericVector b, GenericVector x)

        Compute F at current point x. 
        
";

%feature("docstring")  dolfin::NonlinearProblem::J "

        J(self, GenericMatrix A, GenericVector x)

        Compute J = F' at current point x. 
        
";

%feature("docstring")  dolfin::NonlinearProblem::NonlinearProblem "

        __init__(self) -> NonlinearProblem

        Constructor. 
        
";

%feature("docstring")  dolfin::PETScPreconditioner "
Proxy of C++ dolfin::PETScPreconditioner class
";

%feature("docstring")  dolfin::PETScPreconditioner::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::PETScPreconditioner::set "
set(self, PETScKrylovSolver solver)
";

%feature("docstring")  dolfin::PETScPreconditioner::PETScPreconditioner "

        __init__(self, string type = "default") -> PETScPreconditioner
        __init__(self) -> PETScPreconditioner
        
";

%feature("docstring")  dolfin::MPI "

    This class provides utility functions for easy communcation with MPI.

    C++ includes: MPI.h 
    
";

%feature("docstring")  dolfin::MPI::is_receiver "
is_receiver() -> bool
";

%feature("docstring")  dolfin::MPI::barrier "
barrier()
";

%feature("docstring")  dolfin::MPI::distribute "

        distribute(std::vector<(dolfin::uint)> values, std::vector<(dolfin::uint)> partition)
        distribute(std::vector<(dolfin::uint)> partition)
        
";

%feature("docstring")  dolfin::MPI::is_broadcaster "
is_broadcaster() -> bool
";

%feature("docstring")  dolfin::MPI::local_range "

        local_range(uint N) -> std::pair<(dolfin::uint,dolfin::uint)>
        local_range(uint process, uint N) -> std::pair<(dolfin::uint,dolfin::uint)>
        
";

%feature("docstring")  dolfin::MPI::global_offset "
global_offset(uint range, bool exclusive) -> uint
";

%feature("docstring")  dolfin::MPI::num_processes "
num_processes() -> uint
";

%feature("docstring")  dolfin::MPI::MPI "
__init__(self) -> MPI
";

%feature("docstring")  dolfin::MPI::send_recv "

        send_recv(uint send_buffer, uint send_size, uint dest, uint recv_buffer, 
            uint recv_size, uint source) -> uint
        send_recv(double send_buffer, uint send_size, uint dest, double recv_buffer, 
            uint recv_size, uint source) -> uint
        
";

%feature("docstring")  dolfin::MPI::sum "

        sum(double value) -> double
        sum(uint value) -> uint
        
";

%feature("docstring")  dolfin::MPI::gather "

        gather(uint value) -> std::vector<(dolfin::uint)>
        gather(std::vector<(dolfin::uint)> values)
        gather()
        
";

%feature("docstring")  dolfin::MPI::global_maximum "
global_maximum(uint size) -> uint
";

%feature("docstring")  dolfin::MPI::index_owner "
index_owner(uint index, uint N) -> uint
";

%feature("docstring")  dolfin::MPI::scatter "

        scatter(std::vector<(dolfin::uint)> values, uint sending_process = 0)
        scatter(std::vector<(dolfin::uint)> values)
        scatter(std::vector<(std::vector<(dolfin::uint)>)> values, 
            uint sending_process = 0)
        scatter(std::vector<(std::vector<(dolfin::uint)>)> values)
        scatter(std::vector<(std::vector<(double)>)> values, uint sending_process = 0)
        scatter(std::vector<(std::vector<(double)>)> values)
        
";

%feature("docstring")  dolfin::MPI::process_number "
process_number() -> uint
";

%feature("docstring")  dolfin::Vector "

    This class provides the default DOLFIN vector class, based on the
    default DOLFIN linear algebra backend.

    C++ includes: Vector.h 
    
";

%feature("docstring")  dolfin::Vector::sum "

        sum(self) -> double
        sum(self, dolfin::Array<(dolfin::uint)> rows) -> double

        Return sum of selected rows in vector. Repeated entries only summed
        once. 
        
";

%feature("docstring")  dolfin::Vector::_assign "

        _assign(self, GenericVector x) -> GenericVector
        _assign(self, double a) -> Vector
        _assign(self, Vector x) -> Vector
        
";

%feature("docstring")  dolfin::Vector::get_local "

        get_local(self, double block, uint m)
        get_local(self, DoubleArray values)

        Get all values on local process. 
        
";

%feature("docstring")  dolfin::Vector::_data "
_data(self) -> PyObject
";

%feature("docstring")  dolfin::Vector::copy "

        copy(self) -> Vector

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::Vector::data "
Return an array to the underlaying data
";

%feature("docstring")  dolfin::Vector::Vector "

        __init__(self) -> Vector
        __init__(self, uint N) -> Vector
        __init__(self, Vector x) -> Vector
        __init__(self, GenericVector x) -> Vector

        Create a Vector from a GenericVetor. 
        
";

%feature("docstring")  dolfin::MeshEntity "

    A MeshEntity represents a mesh entity associated with a specific
    topological dimension of some mesh.

    C++ includes: MeshEntity.h 
    
";

%feature("docstring")  dolfin::MeshEntity::intersects "

        intersects(self, Point point) -> bool
        intersects(self, MeshEntity entity) -> bool

        Check if given entity intersects (using inexact but fast numerics). 
        
";

%feature("docstring")  dolfin::MeshEntity::dim "

        dim(self) -> uint

        Return topological dimension. 
        
";

%feature("docstring")  dolfin::MeshEntity::index "

        index(self) -> uint
        index(self, MeshEntity entity) -> uint

        Compute local index of given incident entity (error if not found). 
        
";

%feature("docstring")  dolfin::MeshEntity::num_entities "

        num_entities(self, uint dim) -> uint

        Return number of incident mesh entities of given topological
        dimension. 
        
";

%feature("docstring")  dolfin::MeshEntity::midpoint "

        midpoint(self) -> Point

        Compute midpoint of cell. 
        
";

%feature("docstring")  dolfin::MeshEntity::entities "
 Return number of incident mesh entities of given topological dimension
";

%feature("docstring")  dolfin::MeshEntity::incident "

        incident(self, MeshEntity entity) -> bool

        Check if given entity is indicent. 
        
";

%feature("docstring")  dolfin::MeshEntity::mesh "

        mesh(self) -> Mesh

        Return mesh associated with mesh entity. 
        
";

%feature("docstring")  dolfin::MeshEntity::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::MeshEntity::intersects_exactly "

        intersects_exactly(self, Point point) -> bool
        intersects_exactly(self, MeshEntity entity) -> bool

        Check if given entity intersects (using exact numerics). 
        
";

%feature("docstring")  dolfin::MeshEntity::MeshEntity "

        __init__(self) -> MeshEntity
        __init__(self, Mesh mesh, uint dim, uint index) -> MeshEntity

        Constructor. 
        
";

%feature("docstring")  dolfin::SLEPcEigenSolver "
Proxy of C++ dolfin::SLEPcEigenSolver class
";

%feature("docstring")  dolfin::SLEPcEigenSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::SLEPcEigenSolver::_get_eigenvalue "
_get_eigenvalue(self, int i) -> PyObject
";

%feature("docstring")  dolfin::SLEPcEigenSolver::get_iteration_number "
get_iteration_number(self) -> int
";

%feature("docstring")  dolfin::SLEPcEigenSolver::get_eigenpair "
Gets the i-th solution of the eigenproblem
";

%feature("docstring")  dolfin::SLEPcEigenSolver::_get_eigenpair "
_get_eigenpair(self, PETScVector r, PETScVector c, int i) -> PyObject
";

%feature("docstring")  dolfin::SLEPcEigenSolver::solve "

        solve(self, PETScMatrix A)
        solve(self, PETScMatrix A, uint n)
        solve(self, PETScMatrix A, PETScMatrix B)
        solve(self, PETScMatrix A, PETScMatrix B, uint n)
        
";

%feature("docstring")  dolfin::SLEPcEigenSolver::SLEPcEigenSolver "
__init__(self) -> SLEPcEigenSolver
";

%feature("docstring")  dolfin::SLEPcEigenSolver::get_number_converged "
get_number_converged(self) -> int
";

%feature("docstring")  dolfin::SLEPcEigenSolver::get_eigenvalue "
Gets the i-th eigenvalue of the eigenproblem
";

%feature("docstring")  dolfin::uBLASVector "

    This class provides a simple vector class based on uBLAS. It is a
    simple wrapper for a uBLAS vector implementing the GenericVector
    interface.

    The interface is intentionally simple. For advanced usage, access the
    underlying uBLAS vector and use the standard uBLAS interface which is
    documented athttp://www.boost.org/libs/numeric/ublas/doc/index.htm.

    C++ includes: uBLASVector.h 
    
";

%feature("docstring")  dolfin::uBLASVector::_assign "

        _assign(self, GenericVector x) -> GenericVector
        _assign(self, double a) -> uBLASVector
        _assign(self, uBLASVector x) -> uBLASVector
        
";

%feature("docstring")  dolfin::uBLASVector::vec "

        vec(self) -> ublas_vector
        vec(self) -> ublas_vector

        Return reference to uBLAS vector (non-const version). 
        
";

%feature("docstring")  dolfin::uBLASVector::_data "
_data(self) -> PyObject
";

%feature("docstring")  dolfin::uBLASVector::copy "

        copy(self) -> uBLASVector

        Create copy of tensor. 
        
";

%feature("docstring")  dolfin::uBLASVector::data "
Return an array to the underlaying data
";

%feature("docstring")  dolfin::uBLASVector::uBLASVector "

        __init__(self) -> uBLASVector
        __init__(self, uint N) -> uBLASVector
        __init__(self, uBLASVector x) -> uBLASVector
        __init__(self, boost::shared_ptr<(dolfin::ublas_vector)> x) -> uBLASVector

        Construct vector from a ublas_vector. 
        
";

%feature("docstring")  dolfin::FacetCell "

    This class represents a cell in a mesh incident to a facet on the
    boundary. It is useful in cases where one needs to iterate over a
    boundary mesh and access the corresponding cells in the original mesh.

    C++ includes: FacetCell.h 
    
";

%feature("docstring")  dolfin::FacetCell::FacetCell "

        __init__(self, Mesh mesh, Cell facet) -> FacetCell

        Create cell on mesh corresponding to given facet (cell) on boundary.

        
";

%feature("docstring")  dolfin::FacetCell::facet_index "

        facet_index(self) -> uint

        Return local index of facet with respect to the cell. 
        
";

%feature("docstring")  dolfin::Constant "

    This class represents a constant-valued expression.

    C++ includes: Constant.h 
    
";

%feature("docstring")  dolfin::Constant::assign "

        assign(self, Constant constant) -> Constant
        assign(self, double constant) -> Constant
        
";

%feature("docstring")  dolfin::Constant::Constant "

        __init__(self, double value) -> Constant
        __init__(self, double value0, double value1) -> Constant
        __init__(self, double value0, double value1, double value2) -> Constant
        __init__(self, std::vector<(double)> values) -> Constant
        __init__(self, std::vector<(dolfin::uint)> value_shape, std::vector<(double)> values) -> Constant
        __init__(self, Constant constant) -> Constant

        Copy constructor. 
        
";

%feature("docstring")  dolfin::Cell "

    A Cell is a MeshEntity of topological codimension 0.

    C++ includes: Cell.h 
    
";

%feature("docstring")  dolfin::Cell::diameter "

        diameter(self) -> double

        Compute diameter of cell. 
        
";

%feature("docstring")  dolfin::Cell::ordered "

        ordered(self, MeshFunctionUInt global_vertex_indices) -> bool

        Check if entities are ordered. 
        
";

%feature("docstring")  dolfin::Cell::orientation "

        orientation(self) -> double

        Compute orientation of cell (0 is right, 1 is left). 
        
";

%feature("docstring")  dolfin::Cell::normal "

        normal(self, uint facet, uint i) -> double
        normal(self, uint facet) -> Point

        Compute normal of given facet with respect to the cell. 
        
";

%feature("docstring")  dolfin::Cell::facet_area "

        facet_area(self, uint facet) -> double

        Compute the area/length of given facet with respect to the cell. 
        
";

%feature("docstring")  dolfin::Cell::volume "

        volume(self) -> double

        Compute (generalized) volume of cell. 
        
";

%feature("docstring")  dolfin::Cell::type "

        type(self) -> Type

        Return type of cell. 
        
";

%feature("docstring")  dolfin::Cell::order "

        order(self, MeshFunctionUInt global_vertex_indices)

        Order entities locally. 
        
";

%feature("docstring")  dolfin::Cell::Cell "

        __init__(self) -> Cell
        __init__(self, Mesh mesh, uint index) -> Cell

        Create cell on given mesh with given index. 
        
";

%feature("docstring")  dolfin::SingularSolver "

    This class provides a linear solver for singular linear systems Ax = b
    where A has a one-dimensional null-space (kernel). This may happen for
    example when solving Poisson's equation with pure Neumann boundary
    conditions.

    The solver attempts to create an extended non-singular system by
    adding the constraint [1, 1, 1, ...]^T x = 0.

    If an optional mass matrix M is supplied, the solver attempts to
    create an extended non-singular system by adding the constraint m^T x
    = 0 where m is the lumped mass matrix. This corresponds to setting the
    average (integral) of the finite element function with coefficients x
    to zero.

    The solver makes not attempt to check that the null-space is indeed
    one-dimensional. It is also assumed that the system Ax = b retains its
    sparsity pattern between calls to solve().

    C++ includes: SingularSolver.h 
    
";

%feature("docstring")  dolfin::SingularSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::SingularSolver::solve "

        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b, 
            GenericMatrix M) -> uint

        Solve linear system Ax = b using mass matrix M for setting constraint.

        
";

%feature("docstring")  dolfin::SingularSolver::SingularSolver "

        __init__(self, string solver_type = "lu", string pc_type = "ilu") -> SingularSolver
        __init__(self, string solver_type = "lu") -> SingularSolver
        __init__(self) -> SingularSolver

        Create linear solver. 
        
";

%feature("docstring")  dolfin::IntArray "

    This class provides a simple wrapper for a pointer to an array. A
    purpose of this class is to enable the simple and safe exchange of
    data between C++ and Python.

    C++ includes: Array.h 
    
";

%feature("docstring")  dolfin::IntArray::min "

        min(self) -> int

        Return minimum value of array. 
        
";

%feature("docstring")  dolfin::IntArray::zero_eps "

        zero_eps(self, double eps = 3.0e-16)
        zero_eps(self)
        
";

%feature("docstring")  dolfin::IntArray::update "

        update(self, uint N, int _x)

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::IntArray::zero "

        zero(self)

        Zero array. 
        
";

%feature("docstring")  dolfin::IntArray::resize "

        resize(self, uint N)

        Resize array to size N. If size changes, contents will be destroyed.

        
";

%feature("docstring")  dolfin::IntArray::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). Note that the
        Array class is not a subclass of Variable (for efficiency) which means
        that one needs to call str() directly instead of using the info()
        function on Array objects. 
        
";

%feature("docstring")  dolfin::IntArray::max "

        max(self) -> int

        Return maximum value of array. 
        
";

%feature("docstring")  dolfin::IntArray::array "
array(self) -> PyObject
";

%feature("docstring")  dolfin::IntArray::data "

        data(self) -> boost::shared_array<(int)>
        data(self) -> boost::shared_array<(int)>

        Return pointer to data (non-const version). 
        
";

%feature("docstring")  dolfin::IntArray::IntArray "

        __init__(self) -> IntArray
        __init__(self, uint N) -> IntArray
        __init__(self, uint N) -> IntArray

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::IntArray::size "

        size(self) -> uint

        Return size of array. 
        
";

%feature("docstring")  dolfin::DummyComplexODE "
Proxy of C++ dolfin::DummyComplexODE class
";

%feature("docstring")  dolfin::DummyComplexODE::DummyComplexODE "
__init__(self, uint n, real T) -> DummyComplexODE
";

%feature("docstring")  dolfin::MTL4Vector "
Proxy of C++ dolfin::MTL4Vector class
";

%feature("docstring")  dolfin::MTL4Vector::copy "
copy(self) -> MTL4Vector
";

%feature("docstring")  dolfin::MTL4Vector::_data "
_data(self) -> PyObject
";

%feature("docstring")  dolfin::MTL4Vector::_assign "

        _assign(self, double a) -> MTL4Vector
        _assign(self, GenericVector x) -> GenericVector
        _assign(self, MTL4Vector x) -> MTL4Vector
        
";

%feature("docstring")  dolfin::MTL4Vector::data "
Return an array to the underlaying data
";

%feature("docstring")  dolfin::MTL4Vector::MTL4Vector "

        __init__(self) -> MTL4Vector
        __init__(self, uint N) -> MTL4Vector
        __init__(self, MTL4Vector x) -> MTL4Vector
        
";

%feature("docstring")  dolfin::STLFactory "
Proxy of C++ dolfin::STLFactory class
";

%feature("docstring")  dolfin::STLFactory::instance "
instance() -> STLFactory
";

%feature("docstring")  dolfin::STLFactory::create_matrix "

        create_matrix(self) -> STLMatrix

        Create empty matrix. 
        
";

%feature("docstring")  dolfin::STLFactory::create_vector "

        create_vector(self) -> uBLASVector

        Create empty vector (global). 
        
";

%feature("docstring")  dolfin::STLFactory::create_local_vector "

        create_local_vector(self) -> uBLASVector

        Create empty vector (local). 
        
";

%feature("docstring")  dolfin::STLFactory::STLFactory "
No constructor defined
";

%feature("docstring")  dolfin::uBLASKrylovSolver "

    This class implements Krylov methods for linear systems of the form Ax
    = b using uBLAS data types.

    C++ includes: uBLASKrylovSolver.h 
    
";

%feature("docstring")  dolfin::uBLASKrylovSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::uBLASKrylovSolver::solve "

        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint
        solve(self, uBLASDenseMatrix A, uBLASVector x, uBLASVector b) -> uint
        solve(self, uBLASSparseMatrix A, uBLASVector x, uBLASVector b) -> uint
        solve(self, uBLASKrylovMatrix A, uBLASVector x, uBLASVector b) -> uint

        Solve linear system Ax = b and return number of iterations (virtual
        matrix). 
        
";

%feature("docstring")  dolfin::uBLASKrylovSolver::uBLASKrylovSolver "

        __init__(self, string solver_type = "default", string pc_type = "default") -> uBLASKrylovSolver
        __init__(self, string solver_type = "default") -> uBLASKrylovSolver
        __init__(self) -> uBLASKrylovSolver
        __init__(self, uBLASPreconditioner pc) -> uBLASKrylovSolver
        __init__(self, string solver_type, uBLASPreconditioner preconditioner) -> uBLASKrylovSolver

        Create Krylov solver for a particular method and uBLASPreconditioner.

        
";

%feature("docstring")  dolfin::GenericDofMap "

    This class provides a generic interface for dof maps.

    C++ includes: GenericDofMap.h 
    
";

%feature("docstring")  dolfin::GenericDofMap::needs_mesh_entities "

        needs_mesh_entities(self, unsigned int d) -> bool

        Return true iff mesh entities of topological dimension d are needed.

        
";

%feature("docstring")  dolfin::GenericDofMap::collapse "

        collapse(self, std::map<(dolfin::uint,dolfin::uint)> collapsed_map, 
            Mesh dolfin_mesh) -> GenericDofMap

        "Collapse" a sub dofmap 
        
";

%feature("docstring")  dolfin::GenericDofMap::extract_sub_dofmap "

        extract_sub_dofmap(self, std::vector<(dolfin::uint)> component, Mesh dolfin_mesh) -> GenericDofMap

        Extract sub dofmap component. 
        
";

%feature("docstring")  dolfin::GenericDofMap::tabulate_coordinates "

        tabulate_coordinates(self, double coordinates, cell ufc_cell)
        tabulate_coordinates(self, double coordinates, Cell cell)

        Tabulate the coordinates of all dofs on a cell (DOLFIN cell version).

        
";

%feature("docstring")  dolfin::GenericDofMap::num_facet_dofs "

        num_facet_dofs(self) -> unsigned int

        Return number of facet dofs. 
        
";

%feature("docstring")  dolfin::GenericDofMap::GenericDofMap "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::GenericDofMap::tabulate_facet_dofs "

        tabulate_facet_dofs(self, uint dofs, uint local_facet)

        Tabulate local-local facet dofs. 
        
";

%feature("docstring")  dolfin::GenericDofMap::global_dimension "

        global_dimension(self) -> unsigned int

        Return the dimension of the global finite element function space. 
        
";

%feature("docstring")  dolfin::GenericDofMap::geometric_dimension "
geometric_dimension(self) -> unsigned int
";

%feature("docstring")  dolfin::GenericDofMap::local_dimension "

        local_dimension(self, cell cell) -> unsigned int

        Return the dimension of the local finite element function space on a
        cell 
        
";

%feature("docstring")  dolfin::GenericDofMap::tabulate_dofs "

        tabulate_dofs(self, uint dofs, cell ufc_cell, uint cell_index)
        tabulate_dofs(self, uint dofs, Cell cell)

        Tabulate the local-to-global mapping of dofs on a cell (DOLFIN cell
        version) 
        
";

%feature("docstring")  dolfin::GenericDofMap::str "
Return a string representation of it self
";

%feature("docstring")  dolfin::GenericDofMap::signature "

        signature(self) -> string

        Return a string identifying the dof map. 
        
";

%feature("docstring")  dolfin::GenericDofMap::max_local_dimension "

        max_local_dimension(self) -> unsigned int

        Return the maximum dimension of the local finite element function
        space. 
        
";

%feature("docstring")  dolfin::GenericDofMap::dofs "

        dofs(self, Mesh mesh, bool sort = False) -> dolfin::Set<(dolfin::uint)>
        dofs(self, Mesh mesh) -> dolfin::Set<(dolfin::uint)>

        Return the set of dof indices. 
        
";

%feature("docstring")  dolfin::Method "

    Base class for cGqMethod and dGqMethod, which contain all numeric
    constants, such as nodal points and nodal weights, needed for the
    method.

    C++ includes: Method.h 
    
";

%feature("docstring")  dolfin::Method::qsize "

        qsize(self) -> unsigned int

        Return number of quadrature points (inline optimized). 
        
";

%feature("docstring")  dolfin::Method::degree "

        degree(self) -> unsigned int

        Return degree (inline optimized). 
        
";

%feature("docstring")  dolfin::Method::timestep "

        timestep(self, real r, real tol, real k0, real kmax) -> real

        Compute new time step based on the given residual. 
        
";

%feature("docstring")  dolfin::Method::update "

        update(self, real x0, real f, real k, real values)
        update(self, real x0, real f, real k, real values, real alpha)

        Update solution values using fixed-point iteration (damped version).

        
";

%feature("docstring")  dolfin::Method::get_trial "

        get_trial(self) -> Lagrange

        Get trial functions. 
        
";

%feature("docstring")  dolfin::Method::npoint "

        npoint(self, unsigned int i) -> real

        Return nodal point (inline optimized). 
        
";

%feature("docstring")  dolfin::Method::residual "

        residual(self, real x0, real values, real f, real k) -> real

        Compute residual at right end-point. 
        
";

%feature("docstring")  dolfin::Method::ueval "

        ueval(self, real x0, real values, real tau) -> real
        ueval(self, real x0, real values, uint i) -> real

        Evaluate solution at given node. 
        
";

%feature("docstring")  dolfin::Method::nsize "

        nsize(self) -> unsigned int

        Return number of nodal points (inline optimized). 
        
";

%feature("docstring")  dolfin::Method::eval "

        eval(self, unsigned int i, real tau) -> real

        Evaluation of trial space basis function i at given tau (inline
        optimized). 
        
";

%feature("docstring")  dolfin::Method::qweight "

        qweight(self, unsigned int i) -> real

        Return quadrature weight, including only quadrature (inline
        optimized). 
        
";

%feature("docstring")  dolfin::Method::Method "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::Method::get_quadrature_weights "

        get_quadrature_weights(self) -> real

        Get quadrature weights. 
        
";

%feature("docstring")  dolfin::Method::nweight "

        nweight(self, unsigned int i, unsigned int j) -> real

        Return nodal weight j for node i, including quadrature (inline
        optimized). 
        
";

%feature("docstring")  dolfin::Method::qpoint "

        qpoint(self, unsigned int i) -> real

        Return quadrature point (inline optimized). 
        
";

%feature("docstring")  dolfin::Method::error "

        error(self, real k, real r) -> real

        Compute error estimate (modulo stability factor). 
        
";

%feature("docstring")  dolfin::Method::derivative "

        derivative(self, unsigned int i) -> real

        Evaluation of derivative of basis function i at t = 1 (inline
        optimized). 
        
";

%feature("docstring")  dolfin::Method::order "

        order(self) -> unsigned int

        Return order (inline optimized). 
        
";

%feature("docstring")  dolfin::Method::type "

        type(self) -> Type

        Return type (inline optimized). 
        
";

%feature("docstring")  dolfin::Method::get_nodal_values "

        get_nodal_values(self, real x0, real x, real nodal_values)

        Get nodal values. 
        
";

%feature("docstring")  dolfin::Rectangle "

    Triangular mesh of the 2D rectangle (x0, y0) x (x1, y1). Given the
    number of cells (nx, ny) in each direction, the total number of
    triangles will be 2*nx*ny and the total number of vertices will be (nx
    + 1)*(ny + 1).

    std::string diagonal ("left", "right", "right/left" or
    "crossed") indicates the direction of the diagonals.

    C++ includes: Rectangle.h 
    
";

%feature("docstring")  dolfin::Rectangle::Rectangle "

        __init__(self, double x0, double y0, double x1, double y1, uint nx, 
            uint ny, string diagonal = "right") -> Rectangle
        __init__(self, double x0, double y0, double x1, double y1, uint nx, 
            uint ny) -> Rectangle
        
";

%feature("docstring")  dolfin::Lagrange "

    Lagrange polynomial (basis) with given degree q determined by n = q +
    1 nodal points.

    Example: q = 1 (n = 2)

    Lagrange p(1); p.set(0, 0.0); p.set(1, 1.0);

    This creates a Lagrange polynomial (actually two Lagrange
    polynomials):

    p(0,x) = 1 - x (one at x = 0, zero at x = 1) p(1,x) = x (zero at x =
    0, one at x = 1)

    C++ includes: Lagrange.h 
    
";

%feature("docstring")  dolfin::Lagrange::set "

        set(self, unsigned int i, real x)

        Specify point. 
        
";

%feature("docstring")  dolfin::Lagrange::dqdx "

        dqdx(self, unsigned int i) -> real

        Return derivative q (a constant) of polynomial. 
        
";

%feature("docstring")  dolfin::Lagrange::degree "

        degree(self) -> unsigned int

        Return degree. 
        
";

%feature("docstring")  dolfin::Lagrange::point "

        point(self, unsigned int i) -> real

        Return point. 
        
";

%feature("docstring")  dolfin::Lagrange::ddx "

        ddx(self, unsigned int i, real x) -> real

        Return derivate of polynomial i at given point x. 
        
";

%feature("docstring")  dolfin::Lagrange::eval "

        eval(self, unsigned int i, real x) -> real

        Return value of polynomial i at given point x. 
        
";

%feature("docstring")  dolfin::Lagrange::Lagrange "

        __init__(self, unsigned int q) -> Lagrange
        __init__(self, Lagrange p) -> Lagrange

        Copy constructor. 
        
";

%feature("docstring")  dolfin::Lagrange::size "

        size(self) -> unsigned int

        Return number of points. 
        
";

%feature("docstring")  dolfin::cGqMethod "

    Contains all numeric constants, such as nodal points and nodal
    weights, needed for the cG(q) method. The order q must be at least 1.
    Note that q refers to the polynomial order and not the order of
    convergence for the method, which is 2q.

    C++ includes: cGqMethod.h 
    
";

%feature("docstring")  dolfin::cGqMethod::ueval "

        ueval(self, real x0, real values, real tau) -> real
        ueval(self, real x0, real values, uint i) -> real

        Evaluate solution at given node (inline optimized). 
        
";

%feature("docstring")  dolfin::cGqMethod::cGqMethod "
__init__(self, unsigned int q) -> cGqMethod
";

%feature("docstring")  dolfin::LocalMeshData "

    This class stores mesh data on a local processor corresponding to a
    portion of a (larger) global mesh.

    Note that the data stored in this class does typically not correspond
    to a topologically connected mesh; it merely stores a list of vertex
    coordinates, a list of cell-vertex mappings and a list of global
    vertex numbers for the locally stored vertices.

    It is typically used for parsing meshes in parallel from mesh XML
    files. After local mesh data has been parsed on each processor, a
    subsequent repartitioning takes place: first a geometric partitioning
    of the vertices followed by a redistribution of vertex and cell data,
    and then a topological partitioning again followed by redistribution
    of vertex and cell data, at that point corresponding to topologically
    connected meshes instead of local mesh data.

    C++ includes: LocalMeshData.h 
    
";

%feature("docstring")  dolfin::LocalMeshData::LocalMeshData "

        __init__(self) -> LocalMeshData
        __init__(self, Mesh mesh) -> LocalMeshData

        Create local mesh data for given mesh. 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix "

    This class provides a simple matrix class based on uBLAS. It is a
    simple wrapper for a uBLAS matrix implementing the GenericMatrix
    interface.

    The interface is intentionally simple. For advanced usage, access the
    underlying uBLAS matrix and use the standard uBLAS interface which is
    documented athttp://www.boost.org/libs/numeric/ublas/doc/index.htm.

    Developer note: specialised member functions must be inlined to avoid
    link errors.

    C++ includes: uBLASMatrix.h 
    
";

%feature("docstring")  dolfin::uBLASSparseMatrix::mat "

        mat(self) -> dolfin::ublas::compressed_matrix<(double,dolfin::ublas::row_major)>
        mat(self) -> dolfin::ublas::compressed_matrix<(double,dolfin::ublas::row_major)>

        Return reference to uBLAS matrix (non-const version). 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::invert "

        invert(self)

        Compute inverse of matrix. 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::compress "

        compress(self)

        Compress matrix (eliminate all non-zeros from a sparse matrix). 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::lump "

        lump(self, uBLASVector m)

        Lump matrix into vector m. 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::solveInPlace "

        solveInPlace(self, uBLASVector x, uBLASVector b)

        Solve Ax = b in-place using uBLAS(A is destroyed). 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::zero "

        zero(self)
        zero(self, uint m)

        Set given rows to zero. 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::solve "

        solve(self, uBLASVector x, uBLASVector b)

        Solve Ax = b out-of-place using uBLAS (A is not destroyed). 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::copy "

        copy(self) -> uBLASSparseMatrix

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::assign "

        assign(self, GenericMatrix A) -> GenericMatrix
        assign(self, uBLASSparseMatrix A) -> uBLASSparseMatrix
        
";

%feature("docstring")  dolfin::uBLASSparseMatrix::uBLASSparseMatrix "

        __init__(self) -> uBLASSparseMatrix
        __init__(self, uint M, uint N) -> uBLASSparseMatrix
        __init__(self, uBLASSparseMatrix A) -> uBLASSparseMatrix

        Create matrix from given uBLAS matrix expression. 
        
";

%feature("docstring")  dolfin::PeriodicBC "

    This class specifies the interface for setting periodic boundary
    conditions for partial differential equations,

    u(x) = u(F^{-1}(x)) on G, u(x) = u(F(x)) on H,

    where F : H --> G is a map from a subdomain H to a subdomain G.

    A periodic boundary condition must be defined by the domain G and the
    map F pulling coordinates back from H to G. The domain and the map are
    both defined by a subclass of SubDomain which must overload both the
    inside() function, which specifies the points of G, and the map()
    function, which specifies the map from the points of H to the points
    of G.

    The implementation is based on matching degrees of freedom on G with
    degrees of freedom on H and only works when the mapping F is bijective
    between the sets of coordinates associated with the two domains. In
    other words, the nodes (degrees of freedom) must be aligned on G and
    H.

    The matching of degrees of freedom is done at the construction of the
    periodic boundary condition and is reused on subsequent applications
    to a linear system. The matching may be recomputed by calling the
    rebuild() function.

    C++ includes: PeriodicBC.h 
    
";

%feature("docstring")  dolfin::PeriodicBC::apply "

        apply(self, GenericMatrix A)
        apply(self, GenericVector b)
        apply(self, GenericMatrix A, GenericVector b)
        apply(self, GenericVector b, GenericVector x)
        apply(self, GenericMatrix A, GenericVector b, GenericVector x)

        Apply boundary condition to a linear system for a nonlinear problem.

        
";

%feature("docstring")  dolfin::PeriodicBC::rebuild "

        rebuild(self)

        Rebuild mapping between dofs. 
        
";

%feature("docstring")  dolfin::PeriodicBC::PeriodicBC "

        __init__(self, __dummy_19__ V, __dummy_59__ sub_domain) -> PeriodicBC

        Create periodic boundary condition for sub domain. 
        
";

%feature("docstring")  dolfin::GenericMatrix "

    This class defines a common interface for matrices.

    C++ includes: GenericMatrix.h 
    
";

%feature("docstring")  dolfin::GenericMatrix::set "

        set(self, double block, uint m, uint n)

        Set block of values. 
        
";

%feature("docstring")  dolfin::GenericMatrix::_scale "
_scale(self, double a)
";

%feature("docstring")  dolfin::GenericMatrix::transpmult "

        transpmult(self, GenericVector x, GenericVector y)

        Matrix-vector product, y = A^T x. 
        
";

%feature("docstring")  dolfin::GenericMatrix::mult "

        mult(self, GenericVector x, GenericVector y)

        Matrix-vector product, y = Ax. 
        
";

%feature("docstring")  dolfin::GenericMatrix::ident_zeros "

        ident_zeros(self)

        Insert one on the diagonal for all zero rows. 
        
";

%feature("docstring")  dolfin::GenericMatrix::zero "

        zero(self)
        zero(self, uint m)

        Set given rows to zero. 
        
";

%feature("docstring")  dolfin::GenericMatrix::array "
Return a numpy array representation of Matrix
";

%feature("docstring")  dolfin::GenericMatrix::setrow "

        setrow(self, uint row)

        Set values for given row on local process. 
        
";

%feature("docstring")  dolfin::GenericMatrix::GenericMatrix "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::GenericMatrix::add "

        add(self, double block, uint m, uint n)

        Add block of values. 
        
";

%feature("docstring")  dolfin::GenericMatrix::axpy "

        axpy(self, double a, GenericMatrix A, bool same_nonzero_pattern)

        Add multiple of given matrix (AXPY operation). 
        
";

%feature("docstring")  dolfin::GenericMatrix::norm "

        norm(self, string norm_type) -> double

        Return norm of matrix. 
        
";

%feature("docstring")  dolfin::GenericMatrix::get "

        get(self, double block, uint m, uint n)

        Get block of values. 
        
";

%feature("docstring")  dolfin::GenericMatrix::getrow "

        getrow(self, uint row)

        Get non-zero values of given row on local process. 
        
";

%feature("docstring")  dolfin::GenericMatrix::_data "
_data(self) -> PyObject
";

%feature("docstring")  dolfin::GenericMatrix::copy "

        copy(self) -> GenericMatrix

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::GenericMatrix::data "
 Return arrays to underlying compresssed row/column storage data 
";

%feature("docstring")  dolfin::GenericMatrix::resize "

        resize(self, uint rank, uint dims)
        resize(self, uint M, uint N)

        Resize matrix to M x N. 
        
";

%feature("docstring")  dolfin::GenericMatrix::ident "

        ident(self, uint m)

        Set given rows to identity matrix. 
        
";

%feature("docstring")  dolfin::GenericMatrix::assign "
assign(self, GenericMatrix x) -> GenericMatrix
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

%feature("docstring")  dolfin::FacetArea "

    This function represents the area/length of a cell facet on a given
    mesh.

    C++ includes: SpecialFunctions.h 
    
";

%feature("docstring")  dolfin::FacetArea::FacetArea "

        __init__(self, Mesh mesh) -> FacetArea

        Constructor. 
        
";

%feature("docstring")  dolfin::Box "

    Tetrahedral mesh of the 3D rectangular prism (x0, y0) x (x1, y1) x
    (x2, y2). Given the number of cells (nx, ny, nz) in each direction,
    the total number of tetrahedra will be 6*nx*ny*nz and the total number
    of vertices will be (nx + 1)*(ny + 1)*(nz + 1).

    C++ includes: Box.h 
    
";

%feature("docstring")  dolfin::Box::Box "

        __init__(self, double x0, double y0, double z0, double x1, double y1, 
            double z1, uint nx, uint ny, uint nz) -> Box
        
";

%feature("docstring")  dolfin::ODE "

    An ODE represents an initial value problem of the form

    u'(t) = f(u(t), t) on [0, T],

    u(0) = u0,

    where u(t) is a vector of length N.

    To define an ODE, a user must create a subclass of ODE and create the
    function u0() defining the initial condition, as well the function f()
    defining the right-hand side.

    DOLFIN provides two types of ODE solvers: a set of standard mono-
    adaptive solvers with equal adaptive time steps for all components as
    well as a set of multi-adaptive solvers with individual and adaptive
    time steps for the different components. The right-hand side f() is
    defined differently for the two sets of methods, with the multi-
    adaptive solvers requiring a component-wise evaluation of the right-
    hand side. Only one right-hand side function f() needs to be defined
    for use of any particular solver.

    It is also possible to solve implicit systems of the form

    M(u(t), t) u'(t) = f(u(t),t) on (0,T],

    u(0) = u0,

    by setting the option "implicit" to true and defining the function
    M().

    Two different solve() functions are provided, one to solve the ODE on
    the time interval [0, T], including the solution of a dual problem for
    error control:

    ode.solve();

    Alternatively, a time interval may be given in which case the solution
    will be computed in a single sweep over the given time interval
    without solution of dual problems:

    ode.solve(t0, t1);

    This mode allows the state to be specified and retrieved in between
    intervals by calling set_state() and get_state().

    C++ includes: ODE.h 
    
";

%feature("docstring")  dolfin::ODE::set_state "

        set_state(self, real u)

        Set state for ODE (only available during interval stepping). 
        
";

%feature("docstring")  dolfin::ODE::JT "

        JT(self, real dx, real dy, real u, real t)

        Compute product dy = tranpose(J) dx for Jacobian J (optional, for dual
        problem). 
        
";

%feature("docstring")  dolfin::ODE::timestep "

        timestep(self, real t, real k0) -> real
        timestep(self, real t, uint i, real k0) -> real

        Time step to use for a given component at a given time t (optional).

        
";

%feature("docstring")  dolfin::ODE::dfdu "

        dfdu(self, real u, real t, uint i, uint j) -> real

        Compute entry of Jacobian (optional). 
        
";

%feature("docstring")  dolfin::ODE::J "

        J(self, real dx, real dy, real u, real t)

        Compute product dy = J dx for Jacobian J (optional). 
        
";

%feature("docstring")  dolfin::ODE::M "

        M(self, real dx, real dy, real u, real t)

        Compute product dy = M dx for implicit system (optional). 
        
";

%feature("docstring")  dolfin::ODE::update "

        update(self, real u, real t, bool end) -> bool

        Update ODE, return false to stop (optional). 
        
";

%feature("docstring")  dolfin::ODE::get_state "

        get_state(self, real u)

        Get state for ODE (only available during interval stepping). 
        
";

%feature("docstring")  dolfin::ODE::endtime "

        endtime(self) -> real

        Return end time (final time T). 
        
";

%feature("docstring")  dolfin::ODE::solve_dual "

        solve_dual(self, ODESolution u)
        solve_dual(self, ODESolution u, ODESolution z)

        Solve dual and save soution in z. 
        
";

%feature("docstring")  dolfin::ODE::analyze_stability "

        analyze_stability(self, uint q, ODESolution u)

        Compute stability factors as function of T (including solving the dual
        problem). The stability factor is the integral of the norm of the q'th
        derivative of the dual. 
        
";

%feature("docstring")  dolfin::ODE::ODE "

        __init__(self, uint N, real T) -> ODE

        Create an ODE of size N with final time T. 
        
";

%feature("docstring")  dolfin::ODE::size "

        size(self) -> uint

        Return number of components N. 
        
";

%feature("docstring")  dolfin::ODE::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::ODE::f "

        f(self, real u, real t, real y)
        f(self, real u, real t, uint i) -> real

        Evaluate right-hand side f_i(u, t), multi-adaptive version (optional).

        
";

%feature("docstring")  dolfin::ODE::analyze_stability_computation "

        analyze_stability_computation(self, ODESolution u)

        Compute stability factors as function of T (including solving the dual
        problem). The stability factor accounts for stability wrt the round-
        off errors. 
        
";

%feature("docstring")  dolfin::ODE::analyze_stability_initial "

        analyze_stability_initial(self, ODESolution u)

        Compute stability factors as function of T (including solving the dual
        problem). The stability factor accounts for stability wrt errors in
        initial data. 
        
";

%feature("docstring")  dolfin::ODE::u0 "

        u0(self, real u)

        Set initial values. 
        
";

%feature("docstring")  dolfin::ODE::analyze_stability_discretization "

        analyze_stability_discretization(self, ODESolution u)

        Compute stability factors as function of T (including solving the dual
        problem). The stability factor accounts for stability wrt the
        discretization scheme. 
        
";

%feature("docstring")  dolfin::ODE::solve "

        solve(self)
        solve(self, real t0, real t1)
        solve(self, ODESolution u)
        solve(self, ODESolution u, real t0, real t1)

        Solve ODE on [t0, t1]. Save solution in u. 
        
";

%feature("docstring")  dolfin::ODE::sparse "

        sparse(self)

        Automatically detect sparsity (optional). 
        
";

%feature("docstring")  dolfin::ODE::time "

        time(self) -> real
        time(self, real t) -> real

        Return real time (might be flipped backwards for dual). 
        
";

%feature("docstring")  dolfin::ODE::save "

        save(self, Sample sample)

        Save sample (optional). 
        
";

%feature("docstring")  dolfin::GenericLUSolver "

    This a base class for LU solvers.

    C++ includes: GenericLUSolver.h 
    
";

%feature("docstring")  dolfin::GenericLUSolver::solve "

        solve(self, GenericVector x, GenericVector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint

        Solve linear system Ax = b. 
        
";

%feature("docstring")  dolfin::GenericLUSolver::GenericLUSolver "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::Interval "

    Interval mesh of the 1D line (a,b). Given the number of cells (nx) in
    the axial direction, the total number of intervals will be nx and the
    total number of vertices will be (nx + 1).

    C++ includes: Interval.h 
    
";

%feature("docstring")  dolfin::Interval::Interval "
__init__(self, uint nx, double a, double b) -> Interval
";

%feature("docstring")  dolfin::GaussQuadrature "

    Gauss (Gauss-Legendre) quadrature on the interval [-1,1]. The n
    quadrature points are given by the zeros of the n:th Legendre Pn(x).

    The quadrature points are computed using Newton's method, and the
    quadrature weights are computed by solving a linear system determined
    by the condition that Gauss quadrature with n points should be exact
    for polynomials of degree 2n-1.

    C++ includes: GaussQuadrature.h 
    
";

%feature("docstring")  dolfin::GaussQuadrature::GaussQuadrature "

        __init__(self, unsigned int n) -> GaussQuadrature

        Create Gauss quadrature with n points. 
        
";

%feature("docstring")  dolfin::Face "

    A Face is a MeshEntity of topological dimension 2.

    C++ includes: Face.h 
    
";

%feature("docstring")  dolfin::Face::Face "

        __init__(self, Mesh mesh, uint index) -> Face

        Constructor. 
        
";

%feature("docstring")  dolfin::Edge "

    An Edge is a MeshEntity of topological dimension 1.

    C++ includes: Edge.h 
    
";

%feature("docstring")  dolfin::Edge::length "

        length(self) -> double

        Compute Euclidean length of edge. 
        
";

%feature("docstring")  dolfin::Edge::Edge "

        __init__(self, Mesh mesh, uint index) -> Edge
        __init__(self, MeshEntity entity) -> Edge

        Create edge from mesh entity. 
        
";

%feature("docstring")  dolfin::PETScKrylovMatrix "
Proxy of C++ dolfin::PETScKrylovMatrix class
";

%feature("docstring")  dolfin::PETScKrylovMatrix::PETScKrylovMatrix "

        __init__(self) -> PETScKrylovMatrix
        __init__(self, uint m, uint n) -> PETScKrylovMatrix
        
";

%feature("docstring")  dolfin::PETScKrylovMatrix::mult "
mult(self, PETScVector x, PETScVector y)
";

%feature("docstring")  dolfin::GenericFunction "

    This is a common base class for functions. Functions can be evaluated
    at a given point and they can be restricted to a given cell in a
    finite element mesh. This functionality is implemented by sub-classes
    that implement the eval() and restrict() functions.

    DOLFIN provides two implementations of the GenericFunction interface
    in the form of the classes Function and Expression.

    Sub-classes may optionally implement the gather() function that will
    be called prior to restriction when running in parallel.

    C++ includes: GenericFunction.h 
    
";

%feature("docstring")  dolfin::GenericFunction::gather "

        gather(self)

        Collect off-process coefficients to prepare for interpolation. 
        
";

%feature("docstring")  dolfin::GenericFunction::eval_data "

        eval_data(self, DoubleArray values, Data data)

        Evaluate function for given data. 
        
";

%feature("docstring")  dolfin::GenericFunction::value_dimension "

        value_dimension(self, uint i) -> uint

        Return value dimension for given axis. 
        
";

%feature("docstring")  dolfin::GenericFunction::value_size "

        value_size(self) -> uint

        Return value size (product of value dimensions). 
        
";

%feature("docstring")  dolfin::GenericFunction::value_rank "

        value_rank(self) -> uint

        Return value rank. 
        
";

%feature("docstring")  dolfin::GenericFunction::restrict "

        restrict(self, double w, FiniteElement element, Cell dolfin_cell, 
            cell ufc_cell, int local_facet)
        restrict(self, double w, FiniteElement element, Cell dolfin_cell, 
            cell ufc_cell)

        Convenience function for restriction when facet is unknown. 
        
";

%feature("docstring")  dolfin::GenericFunction::compute_vertex_values "

        compute_vertex_values(self, DoubleArray vertex_values, Mesh mesh)

        Compute values at all mesh vertices. 
        
";

%feature("docstring")  dolfin::GenericFunction::str "
Return a string representation of it self
";

%feature("docstring")  dolfin::GenericFunction::GenericFunction "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::Expression "

    This class represents a user-defined expression. Expressions can be
    used as coefficients in variational forms or interpolated into finite
    element spaces.

    An expression is defined by overloading the eval() method. Users may
    choose to overload either a simple version of eval(), in the case of
    expressions only depending on the coordinate x, or an optional version
    for expressions depending on x and mesh data like cell indices or
    facet normals.

    The geometric dimension (the size of x) and the value rank and
    dimensions of an expression must supplied as arguments to the
    constructor.

    C++ includes: Expression.h 
    
";

%feature("docstring")  dolfin::Expression::eval "

        eval(self, DoubleArray values, DoubleArray x)

        Evaluate expression, must be overloaded by user (simple version). 
        
";

%feature("docstring")  dolfin::Expression::eval_data "

        eval_data(self, DoubleArray values, Data data)

        Evaluate expression, must be overloaded by user (simple version). 
        
";

%feature("docstring")  dolfin::Expression::Expression "

        __init__(self) -> Expression
        __init__(self, uint dim) -> Expression
        __init__(self, uint dim0, uint dim1) -> Expression
        __init__(self, std::vector<(dolfin::uint)> value_shape) -> Expression
        __init__(self, Expression expression) -> Expression

        Copy constructor. 
        
";

%feature("docstring")  dolfin::BlockMatrix "
Proxy of C++ dolfin::BlockMatrix class
";

%feature("docstring")  dolfin::BlockMatrix::set "

        set(self, uint i, uint j, Matrix m)

        Set block. 
        
";

%feature("docstring")  dolfin::BlockMatrix::get "

        get(self, uint i, uint j) -> Matrix
        get(self, uint i, uint j) -> Matrix

        Get block. 
        
";

%feature("docstring")  dolfin::BlockMatrix::zero "

        zero(self)

        Set all entries to zero and keep any sparse structure. 
        
";

%feature("docstring")  dolfin::BlockMatrix::BlockMatrix "

        __init__(self, uint n = 0, uint m = 0, bool owner = False) -> BlockMatrix
        __init__(self, uint n = 0, uint m = 0) -> BlockMatrix
        __init__(self, uint n = 0) -> BlockMatrix
        __init__(self) -> BlockMatrix
        
";

%feature("docstring")  dolfin::BlockMatrix::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::BlockMatrix::apply "

        apply(self, string mode)

        Finalize assembly of tensor. 
        
";

%feature("docstring")  dolfin::BlockMatrix::mult "

        mult(self, BlockVector x, BlockVector y, bool transposed = False)
        mult(self, BlockVector x, BlockVector y)

        Matrix-vector product, y = Ax. 
        
";

%feature("docstring")  dolfin::BlockMatrix::size "

        size(self, uint dim) -> uint

        Return size of given dimension. 
        
";

%feature("docstring")  dolfin::uBLASKrylovMatrix "

    This class provides an interface for matrices that define linear
    systems for the uBLASKrylovSolver. This interface is implemented by
    the classes uBLASSparseMatrix and DenseMatrix. Users may also overload
    the mult() function to specify a linear system only in terms of its
    action.

    C++ includes: uBLASKrylovMatrix.h 
    
";

%feature("docstring")  dolfin::uBLASKrylovMatrix::solve "

        solve(self, uBLASVector x, uBLASVector b)

        Solve linear system Ax = b for a Krylov matrix using uBLAS and dense
        matrices. 
        
";

%feature("docstring")  dolfin::uBLASKrylovMatrix::mult "

        mult(self, uBLASVector x, uBLASVector y)

        Compute product y = Ax. 
        
";

%feature("docstring")  dolfin::uBLASKrylovMatrix::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). 
        
";

%feature("docstring")  dolfin::uBLASKrylovMatrix::uBLASKrylovMatrix "

        __init__(self) -> uBLASKrylovMatrix

        Constructor. 
        
";

%feature("docstring")  dolfin::uBLASKrylovMatrix::size "

        size(self, uint dim) -> uint

        Return number of rows (dim = 0) or columns (dim = 1). 
        
";

%feature("docstring")  dolfin::Data "

    This class holds data for function evaluation, including the
    coordinates x, the time t, and auxiliary data that a function may
    depend on.

    C++ includes: Data.h 
    
";

%feature("docstring")  dolfin::Data::cell "

        cell(self) -> Cell

        Return current cell (if available). 
        
";

%feature("docstring")  dolfin::Data::set "

        set(self, Cell dolfin_cell, cell ufc_cell, int local_facet)
        set(self, cell ufc_cell, double x)
        set(self, uint gdim, double x)
        
";

%feature("docstring")  dolfin::Data::normal "

        normal(self) -> Point

        Return current facet normal (if available). 
        
";

%feature("docstring")  dolfin::Data::clear "

        clear(self)

        Clear all cell data. 
        
";

%feature("docstring")  dolfin::Data::on_facet "

        on_facet(self) -> bool

        Check if we are on a facet. 
        
";

%feature("docstring")  dolfin::Data::geometric_dimension "

        geometric_dimension(self) -> uint

        Return geometric dimension of cell. 
        
";

%feature("docstring")  dolfin::Data::ufc_cell "

        ufc_cell(self) -> cell

        Return current UFC cell (if available). 
        
";

%feature("docstring")  dolfin::Data::facet "

        facet(self) -> uint

        Return current facet (if available). 
        
";

%feature("docstring")  dolfin::Data::x "
x(self) -> PyObject
";

%feature("docstring")  dolfin::Data::Data "

        __init__(self) -> Data

        Constructor. 
        
";

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

%feature("docstring")  dolfin::PETScMatrix "
Proxy of C++ dolfin::PETScMatrix class
";

%feature("docstring")  dolfin::PETScMatrix::size "

        size(self, uint dim) -> uint

        Return size of given dimension. 
        
";

%feature("docstring")  dolfin::PETScMatrix::zero "

        zero(self)
        zero(self, uint m)

        Set given rows to zero. 
        
";

%feature("docstring")  dolfin::PETScMatrix::copy "
copy(self) -> PETScMatrix
";

%feature("docstring")  dolfin::PETScMatrix::assign "

        assign(self, GenericMatrix A) -> GenericMatrix
        assign(self, PETScMatrix A) -> PETScMatrix
        
";

%feature("docstring")  dolfin::PETScMatrix::PETScMatrix "

        __init__(self) -> PETScMatrix
        __init__(self, uint M, uint N) -> PETScMatrix
        __init__(self, PETScMatrix A) -> PETScMatrix
        __init__(self, boost::shared_ptr<(Mat)> A) -> PETScMatrix
        
";

%feature("docstring")  dolfin::MTL4Matrix "
Proxy of C++ dolfin::MTL4Matrix class
";

%feature("docstring")  dolfin::MTL4Matrix::zero "

        zero(self)
        zero(self, uint m)

        Set given rows to zero. 
        
";

%feature("docstring")  dolfin::MTL4Matrix::copy "
copy(self) -> MTL4Matrix
";

%feature("docstring")  dolfin::MTL4Matrix::assign "

        assign(self, GenericMatrix A) -> GenericMatrix
        assign(self, MTL4Matrix A) -> MTL4Matrix
        
";

%feature("docstring")  dolfin::MTL4Matrix::MTL4Matrix "

        __init__(self) -> MTL4Matrix
        __init__(self, uint M, uint N) -> MTL4Matrix
        __init__(self, MTL4Matrix A) -> MTL4Matrix
        __init__(self, uint M, uint N, uint nz) -> MTL4Matrix
        
";

%feature("docstring")  dolfin::GenericLinearSolver "

    This class provides a general solver for linear systems Ax = b.

    C++ includes: GenericLinearSolver.h 
    
";

%feature("docstring")  dolfin::GenericLinearSolver::set_operators "

        set_operators(self, GenericMatrix A, GenericMatrix P)

        Solve the operator (matrix) and preconditioner matrix. 
        
";

%feature("docstring")  dolfin::GenericLinearSolver::solve "

        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint
        solve(self, GenericVector x, GenericVector b) -> uint

        Solve linear system Ax = b. 
        
";

%feature("docstring")  dolfin::GenericLinearSolver::set_operator "

        set_operator(self, GenericMatrix A)

        Solve the operator (matrix). 
        
";

%feature("docstring")  dolfin::GenericLinearSolver::GenericLinearSolver "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::StabilityAnalysis "
Proxy of C++ dolfin::StabilityAnalysis class
";

%feature("docstring")  dolfin::StabilityAnalysis::analyze_endpoint "

        analyze_endpoint(self)

        Compute z(0) (the endpoint of the dual) as function of (primal)
        endtime T. 
        
";

%feature("docstring")  dolfin::StabilityAnalysis::analyze_integral "

        analyze_integral(self, uint q)

        Compute the integral of the q'th derivative of the dual as function of
        (primal) endtime T. 
        
";

%feature("docstring")  dolfin::StabilityAnalysis::StabilityAnalysis "

        __init__(self, ODE ode, ODESolution u) -> StabilityAnalysis

        Constructor. 
        
";

%feature("docstring")  dolfin::MeshEditor "

    A simple mesh editor for creating simplicial meshes in 1D, 2D and 3D.

    C++ includes: MeshEditor.h 
    
";

%feature("docstring")  dolfin::MeshEditor::add_higher_order_cell_data "

        add_higher_order_cell_data(self, uint c, uint v0, uint v1, uint v2, uint v3, uint v4, 
            uint v5)

        Add higher order cell data (assume P2 triangle for now). 
        
";

%feature("docstring")  dolfin::MeshEditor::init_higher_order_vertices "

        init_higher_order_vertices(self, uint num_higher_order_vertices)

        Specify number of vertices. 
        
";

%feature("docstring")  dolfin::MeshEditor::open "

        open(self, Mesh mesh, uint tdim, uint gdim)
        open(self, Mesh mesh, string type, uint tdim, uint gdim)

        Open mesh of given cell type, topological and geometrical dimension.

        
";

%feature("docstring")  dolfin::MeshEditor::init_vertices "

        init_vertices(self, uint num_vertices)

        Specify number of vertices. 
        
";

%feature("docstring")  dolfin::MeshEditor::set_affine_cell_indicator "

        set_affine_cell_indicator(self, uint c, string affine_str)

        Set boolean indicator inside MeshGeometry. 
        
";

%feature("docstring")  dolfin::MeshEditor::add_cell "

        add_cell(self, uint c, std::vector<(dolfin::uint)> v)
        add_cell(self, uint c, uint v0, uint v1)
        add_cell(self, uint c, uint v0, uint v1, uint v2)
        add_cell(self, uint c, uint v0, uint v1, uint v2, uint v3)

        Add cell (tetrahedron) with given vertices. 
        
";

%feature("docstring")  dolfin::MeshEditor::init_cells "

        init_cells(self, uint num_cells)

        Specify number of cells. 
        
";

%feature("docstring")  dolfin::MeshEditor::init_higher_order_cells "

        init_higher_order_cells(self, uint num_higher_order_cells, uint num_higher_order_cell_dof)

        Specify number of cells. 
        
";

%feature("docstring")  dolfin::MeshEditor::close "

        close(self, bool order = True)
        close(self)

        Close mesh, finish editing, and order entities locally. 
        
";

%feature("docstring")  dolfin::MeshEditor::add_higher_order_vertex "

        add_higher_order_vertex(self, uint v, Point p)
        add_higher_order_vertex(self, uint v, double x)
        add_higher_order_vertex(self, uint v, double x, double y)
        add_higher_order_vertex(self, uint v, double x, double y, double z)

        Add vertex v at given coordinate (x, y, z). 
        
";

%feature("docstring")  dolfin::MeshEditor::add_vertex "

        add_vertex(self, uint v, Point p)
        add_vertex(self, uint v, double x)
        add_vertex(self, uint v, double x, double y)
        add_vertex(self, uint v, double x, double y, double z)

        Add vertex v at given coordinate (x, y, z). 
        
";

%feature("docstring")  dolfin::MeshEditor::MeshEditor "

        __init__(self) -> MeshEditor

        Constructor. 
        
";

%feature("docstring")  dolfin::UIntArray "

    This class provides a simple wrapper for a pointer to an array. A
    purpose of this class is to enable the simple and safe exchange of
    data between C++ and Python.

    C++ includes: Array.h 
    
";

%feature("docstring")  dolfin::UIntArray::min "

        min(self) -> unsigned int

        Return minimum value of array. 
        
";

%feature("docstring")  dolfin::UIntArray::zero_eps "

        zero_eps(self, double eps = 3.0e-16)
        zero_eps(self)
        
";

%feature("docstring")  dolfin::UIntArray::update "

        update(self, uint N, unsigned int _x)

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::UIntArray::zero "

        zero(self)

        Zero array. 
        
";

%feature("docstring")  dolfin::UIntArray::resize "

        resize(self, uint N)

        Resize array to size N. If size changes, contents will be destroyed.

        
";

%feature("docstring")  dolfin::UIntArray::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). Note that the
        Array class is not a subclass of Variable (for efficiency) which means
        that one needs to call str() directly instead of using the info()
        function on Array objects. 
        
";

%feature("docstring")  dolfin::UIntArray::max "

        max(self) -> unsigned int

        Return maximum value of array. 
        
";

%feature("docstring")  dolfin::UIntArray::array "
array(self) -> PyObject
";

%feature("docstring")  dolfin::UIntArray::data "

        data(self) -> boost::shared_array<(unsigned int)>
        data(self) -> boost::shared_array<(unsigned int)>

        Return pointer to data (non-const version). 
        
";

%feature("docstring")  dolfin::UIntArray::UIntArray "

        __init__(self) -> UIntArray
        __init__(self, uint N) -> UIntArray
        __init__(self, uint N) -> UIntArray

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::UIntArray::size "

        size(self) -> uint

        Return size of array. 
        
";

%feature("docstring")  dolfin::SubMesh "

    A SubMesh is a mesh defined as a subset of a given mesh. It provides a
    convenient way to create matching meshes for multiphysics applications
    by creating meshes for subdomains as subsets of a single global mesh.

    C++ includes: SubMesh.h 
    
";

%feature("docstring")  dolfin::SubMesh::SubMesh "

        __init__(self, Mesh mesh, SubDomain sub_domain) -> SubMesh
        __init__(self, Mesh mesh, MeshFunctionUInt sub_domains, uint sub_domain) -> SubMesh

        Create subset of given mesh marked by mesh function. 
        
";

%feature("docstring")  dolfin::SubSpace "

    This class represents a subspace (component) of a function space.

    The subspace is specified by an array of indices. For example, the
    array [3, 0, 2] specifies subspace 2 of subspace 0 of subspace 3.

    A typical example is the function space W = V x P for Stokes. Here, V
    = W[0] is the subspace for the velocity component and P = W[1] is the
    subspace for the pressure component. Furthermore, W[0][0] = V[0] is
    the first component of the velocity space etc.

    C++ includes: SubSpace.h 
    
";

%feature("docstring")  dolfin::SubSpace::SubSpace "

        __init__(self, FunctionSpace V, uint component) -> SubSpace
        __init__(self, FunctionSpace V, uint component, uint sub_component) -> SubSpace
        __init__(self, FunctionSpace V, std::vector<(dolfin::uint)> component) -> SubSpace

        Create subspace for given component (n levels). 
        
";

%feature("docstring")  dolfin::ConstDoubleArray "

    This class provides a simple wrapper for a pointer to an array. A
    purpose of this class is to enable the simple and safe exchange of
    data between C++ and Python.

    C++ includes: Array.h 
    
";

%feature("docstring")  dolfin::ConstDoubleArray::min "

        min(self) -> double

        Return minimum value of array. 
        
";

%feature("docstring")  dolfin::ConstDoubleArray::zero_eps "

        zero_eps(self, double eps = 3.0e-16)
        zero_eps(self)
        
";

%feature("docstring")  dolfin::ConstDoubleArray::update "

        update(self, uint N, double _x)

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::ConstDoubleArray::str "

        str(self, bool verbose) -> string

        Return informal string representation (pretty-print). Note that the
        Array class is not a subclass of Variable (for efficiency) which means
        that one needs to call str() directly instead of using the info()
        function on Array objects. 
        
";

%feature("docstring")  dolfin::ConstDoubleArray::max "

        max(self) -> double

        Return maximum value of array. 
        
";

%feature("docstring")  dolfin::ConstDoubleArray::data "

        data(self) -> boost::shared_array<(q(const).double)>
        data(self) -> boost::shared_array<(q(const).double)>

        Return pointer to data (non-const version). 
        
";

%feature("docstring")  dolfin::ConstDoubleArray::ConstDoubleArray "

        __init__(self) -> ConstDoubleArray
        __init__(self, uint N) -> ConstDoubleArray

        Construct array from a pointer. Array will not take ownership. 
        
";

%feature("docstring")  dolfin::ConstDoubleArray::size "

        size(self) -> uint

        Return size of array. 
        
";

%feature("docstring")  dolfin::vertices "

    A VertexIterator is a MeshEntityIterator of topological dimension 0.

    C++ includes: Vertex.h 
    
";

%feature("docstring")  dolfin::vertices::_dereference "
_dereference(self) -> Vertex
";

%feature("docstring")  dolfin::vertices::vertices "

        VertexIterator(Mesh mesh) -> vertices
        __init__(self, MeshEntity entity) -> vertices
        
";

%feature("docstring")  dolfin::LobattoQuadrature "

    Lobatto (Gauss-Lobatto) quadrature on the interval [-1,1]. The n
    quadrature points are given by the end-points -1 and 1, and the zeros
    of P{n-1}'(x), where P{n-1}(x) is the (n-1):th Legendre polynomial.

    The quadrature points are computed using Newton's method, and the
    quadrature weights are computed by solving a linear system determined
    by the condition that Lobatto quadrature with n points should be exact
    for polynomials of degree 2n-3.

    C++ includes: LobattoQuadrature.h 
    
";

%feature("docstring")  dolfin::LobattoQuadrature::LobattoQuadrature "

        __init__(self, unsigned int n) -> LobattoQuadrature

        Create Lobatto quadrature with n points. 
        
";

%feature("docstring")  dolfin::MTL4Factory "
Proxy of C++ dolfin::MTL4Factory class
";

%feature("docstring")  dolfin::MTL4Factory::create_matrix "
create_matrix(self) -> MTL4Matrix
";

%feature("docstring")  dolfin::MTL4Factory::create_vector "
create_vector(self) -> MTL4Vector
";

%feature("docstring")  dolfin::MTL4Factory::create_krylov_solver "

        create_krylov_solver(self, string method, string pc) -> ITLKrylovSolver

        Create Krylov solver. 
        
";

%feature("docstring")  dolfin::MTL4Factory::instance "
instance() -> MTL4Factory
";

%feature("docstring")  dolfin::MTL4Factory::create_lu_solver "

        create_lu_solver(self) -> UmfpackLUSolver

        Create LU solver. 
        
";

%feature("docstring")  dolfin::MTL4Factory::create_local_vector "

        create_local_vector(self) -> MTL4Vector

        Create empty vector (local). 
        
";

%feature("docstring")  dolfin::MTL4Factory::MTL4Factory "
No constructor defined
";

%feature("docstring")  dolfin::DirichletBC "

    This class specifies the interface for setting (strong).

    This class specifies the interface for setting (strong) Dirichlet
    boundary conditions for partial differential equations,

    u = g on G,

    where u is the solution to be computed, g is a function and G is a sub
    domain of the mesh.

    A DirichletBC is specified by the function g, the function space
    (trial space) and boundary indicators on (a subset of) the mesh
    boundary.

    The boundary indicators may be specified in a number of different
    ways.

    The simplest approach is to specify the boundary by a SubDomain
    object, using the inside() function to specify on which facets the
    boundary conditions should be applied.

    Alternatively, the boundary may be specified by a MeshFunction
    labeling all mesh facets together with a number that specifies which
    facets should be included in the boundary.

    The third option is to attach the boundary information to the mesh.
    This is handled automatically when exporting a mesh from for example
    VMTK.

    The BCMethod variable may be used to specify the type of method used
    to identify degrees of freedom on the boundary. Available methods are:
    topological approach (default), geometric approach, and pointwise
    approach. The topological approach is faster, but will only identify
    degrees of freedom that are located on a facet that is entirely on the
    boundary. In particular, the topological approach will not identify
    degrees of freedom for discontinuous elements (which are all internal
    to the cell). A remedy for this is to use the geometric approach. To
    apply pointwise boundary conditions e.g. pointloads, one will have to
    use the pointwise approach which in turn is the slowest of the three
    possible methods. The three possibilties are "topological",
    "geometric" and "pointwise".

    C++ includes: DirichletBC.h 
    
";

%feature("docstring")  dolfin::DirichletBC::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::DirichletBC::value_ptr "

        value_ptr(self) -> __dummy_23__

        Return shared pointer to boundary value g Testing multiline comment 
        
";

%feature("docstring")  dolfin::DirichletBC::is_compatible "

        is_compatible(self, GenericFunction v) -> bool

        Check if given function is compatible with boundary condition
        (checking only vertex values). 
        
";

%feature("docstring")  dolfin::DirichletBC::markers "

        markers(self) -> std::vector<(std::pair<(dolfin::uint,dolfin::uint)>)>

        Return boundary markers (facets stored as pairs of cells and local
        facet numbers). 
        
";

%feature("docstring")  dolfin::DirichletBC::zero "

        zero(self, GenericMatrix A)

        Make row associated with boundary conditions zero, useful for non-
        diagonal matrices in a block matrix. 
        
";

%feature("docstring")  dolfin::DirichletBC::get_bc "

        get_bc(self, uint indicators, double values)

        Get Dirichlet values and indicators. 
        
";

%feature("docstring")  dolfin::DirichletBC::value "

        value(self) -> GenericFunction

        Return boundary value g. 
        
";

%feature("docstring")  dolfin::DirichletBC::apply "

        apply(self, GenericMatrix A)
        apply(self, GenericVector b)
        apply(self, GenericMatrix A, GenericVector b)
        apply(self, GenericVector b, GenericVector x)
        apply(self, GenericMatrix A, GenericVector b, GenericVector x)

        Apply boundary condition to a linear system for a nonlinear problem.

        
";

%feature("docstring")  dolfin::DirichletBC::DirichletBC "

        __init__(self, __dummy_19__ V, __dummy_23__ g, __dummy_59__ sub_domain, 
            string method = "topological") -> DirichletBC
        __init__(self, __dummy_19__ V, __dummy_23__ g, __dummy_59__ sub_domain) -> DirichletBC
        __init__(self, __dummy_19__ V, __dummy_23__ g, MeshFunctionUInt sub_domains, 
            uint sub_domain, string method = "topological") -> DirichletBC
        __init__(self, __dummy_19__ V, __dummy_23__ g, MeshFunctionUInt sub_domains, 
            uint sub_domain) -> DirichletBC
        __init__(self, __dummy_19__ V, __dummy_23__ g, uint sub_domain, string method = "topological") -> DirichletBC
        __init__(self, __dummy_19__ V, __dummy_23__ g, uint sub_domain) -> DirichletBC
        __init__(self, FunctionSpace V, GenericFunction g, std::vector<(std::pair<(dolfin::uint,dolfin::uint)>)> markers, 
            string method = "topological") -> DirichletBC
        __init__(self, FunctionSpace V, GenericFunction g, std::vector<(std::pair<(dolfin::uint,dolfin::uint)>)> markers) -> DirichletBC
        __init__(self, __dummy_19__ V, __dummy_23__ g, std::vector<(std::pair<(dolfin::uint,dolfin::uint)>)> markers, 
            string method = "topological") -> DirichletBC
        __init__(self, __dummy_19__ V, __dummy_23__ g, std::vector<(std::pair<(dolfin::uint,dolfin::uint)>)> markers) -> DirichletBC
        __init__(self, DirichletBC bc) -> DirichletBC

        Copy constructor. 
        
";

%feature("docstring")  dolfin::DirichletBC::set_value "

        set_value(self, GenericFunction g)
        set_value(self, __dummy_23__ g)

        Set value g for boundary condition, domain remains unchanged. 
        
";

%feature("docstring")  dolfin::FiniteElement "

    This is a wrapper for a UFC finite element (ufc::finite_element).

    C++ includes: FiniteElement.h 
    
";

%feature("docstring")  dolfin::FiniteElement::hash "

        hash(self) -> uint

        Return simple hash of the signature string. 
        
";

%feature("docstring")  dolfin::FiniteElement::ufc_element "

        ufc_element(self) -> __dummy_3__

        Return ufc::finite_element. 
        
";

%feature("docstring")  dolfin::FiniteElement::evaluate_basis_derivatives "

        evaluate_basis_derivatives(self, unsigned int i, unsigned int n, double values, double x, 
            cell cell)
        
";

%feature("docstring")  dolfin::FiniteElement::value_rank "
value_rank(self) -> uint
";

%feature("docstring")  dolfin::FiniteElement::extract_sub_element "

        extract_sub_element(self, std::vector<(dolfin::uint)> component) -> __dummy_17__

        Extract sub finite element for component. 
        
";

%feature("docstring")  dolfin::FiniteElement::evaluate_dof "
evaluate_dof(self, uint i, function function, cell cell) -> double
";

%feature("docstring")  dolfin::FiniteElement::FiniteElement "

        __init__(self, __dummy_3__ element) -> FiniteElement

        Create finite element from UFC finite element (data may be shared). 
        
";

%feature("docstring")  dolfin::FiniteElement::evaluate_basis "

        evaluate_basis(self, uint i, double values, double x, cell cell)
        evaluate_basis(self, uint i, double values, double x, Cell cell)
        
";

%feature("docstring")  dolfin::FiniteElement::value_dimension "
value_dimension(self, uint i) -> uint
";

%feature("docstring")  dolfin::FiniteElement::space_dimension "
space_dimension(self) -> uint
";

%feature("docstring")  dolfin::FiniteElement::create_sub_element "

        create_sub_element(self, uint i) -> __dummy_17__

        Create sub element. 
        
";

%feature("docstring")  dolfin::FiniteElement::signature "
signature(self) -> string
";

%feature("docstring")  dolfin::FiniteElement::interpolate_vertex_values "
interpolate_vertex_values(self, double vertex_values, double coefficients, cell cell)
";

%feature("docstring")  dolfin::FiniteElement::num_sub_elements "
num_sub_elements(self) -> uint
";

%feature("docstring")  dolfin::DomainBoundary "

    This class provides a SubDomain which picks out the boundary of a
    mesh, and provides a convenient way to specify boundary conditions on
    the entire boundary of a mesh.

    C++ includes: DomainBoundary.h 
    
";

%feature("docstring")  dolfin::DomainBoundary::DomainBoundary "

        __init__(self) -> DomainBoundary

        Constructor. 
        
";

%feature("docstring")  dolfin::StringParameter "

    Parameter with value type string.

    C++ includes: Parameter.h 
    
";

%feature("docstring")  dolfin::StringParameter::_assign "

        _assign(self, string value) -> StringParameter
        _assign(self, char value) -> StringParameter
        
";

%feature("docstring")  dolfin::StringParameter::StringParameter "

        __init__(self, string key, string value) -> StringParameter

        Create string-valued parameter. 
        
";

%feature("docstring")  dolfin::DofMap "

    This class handles the mapping of degrees of freedom. It builds a dof
    map based on a ufc::dof_map on a specific mesh. It will reorder the
    dofs when running in parallel.

    If ufc_offset != 0, then the dof map provides a view into a larger dof
    map. A dof map which is a view, can be 'collapsed' such that the dof
    indices are contiguous.

    C++ includes: DofMap.h 
    
";

%feature("docstring")  dolfin::DofMap::collapse "

        collapse(self, std::map<(dolfin::uint,dolfin::uint)> collapsed_map, 
            Mesh dolfin_mesh) -> DofMap

        "Collapse" a sub dofmap 
        
";

%feature("docstring")  dolfin::DofMap::extract_sub_dofmap "

        extract_sub_dofmap(std::vector<(dolfin::uint)> component, Mesh dolfin_mesh) -> DofMap
        extract_sub_dofmap(dof_map ufc_dof_map, uint offset, std::vector<(dolfin::uint)> component, 
            mesh ufc_mesh, Mesh dolfin_mesh) -> dof_map

        Extract sub dofmap component. 
        
";

%feature("docstring")  dolfin::DofMap::tabulate_coordinates "

        tabulate_coordinates(self, double coordinates, cell ufc_cell)
        tabulate_coordinates(self, double coordinates, Cell cell)

        Tabulate the coordinates of all dofs on a cell (DOLFIN cell version).

        
";

%feature("docstring")  dolfin::DofMap::tabulate_dofs "

        tabulate_dofs(self, uint dofs, cell ufc_cell, uint cell_index)
        tabulate_dofs(self, uint dofs, Cell cell)

        Tabulate the local-to-global mapping of dofs on a cell (DOLFIN cell
        version). 
        
";

%feature("docstring")  dolfin::DofMap::dofs "

        dofs(self, Mesh mesh, bool sort = False) -> dolfin::Set<(dolfin::uint)>
        dofs(self, Mesh mesh) -> dolfin::Set<(dolfin::uint)>

        Return the set of dof indices. 
        
";

%feature("docstring")  dolfin::DofMap::DofMap "

        __init__(self, __dummy_4__ ufc_dofmap, Mesh dolfin_mesh) -> DofMap

        Create dof map on mesh (const mesh version). 
        
";

%feature("docstring")  dolfin::STLMatrix "

    Simple implementation of a GenericMatrix for experimenting with new
    assembly. Not sure this will be used later but it might be useful.

    C++ includes: STLMatrix.h 
    
";

%feature("docstring")  dolfin::STLMatrix::resize "

        resize(self, uint M, uint N)
        resize(self, uint rank, uint dims, bool reset)

        Resize tensor of given rank and dimensions. 
        
";

%feature("docstring")  dolfin::STLMatrix::copy "

        copy(self) -> STLMatrix

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::STLMatrix::zero "

        zero(self)
        zero(self, uint m)

        Set given rows to zero. 
        
";

%feature("docstring")  dolfin::STLMatrix::STLMatrix "

        __init__(self) -> STLMatrix
        __init__(self, uint M, uint N) -> STLMatrix
        __init__(self, STLMatrix A) -> STLMatrix

        Copy constructor. 
        
";

%feature("docstring")  dolfin::BarycenterQuadrature "
Proxy of C++ dolfin::BarycenterQuadrature class
";

%feature("docstring")  dolfin::BarycenterQuadrature::points "
points(self) -> std::vector<(dolfin::Point)>
";

%feature("docstring")  dolfin::BarycenterQuadrature::weights "
weights(self) -> std::vector<(double)>
";

%feature("docstring")  dolfin::BarycenterQuadrature::BarycenterQuadrature "
__init__(self, Nef_polyhedron_3 polyhedron) -> BarycenterQuadrature
";

%feature("docstring")  dolfin::BarycenterQuadrature::size "
size(self) -> uint
";

%feature("docstring")  dolfin::FunctionPlotData "

    This class is used for communicating plot data for functions to and
    from (XML) files. It is used by DOLFIN for plotting Function objects.
    The data is stored as a mesh and a vector of interpolated vertex
    values.

    C++ includes: FunctionPlotData.h 
    
";

%feature("docstring")  dolfin::FunctionPlotData::FunctionPlotData "

        __init__(self, GenericFunction v, Mesh mesh) -> FunctionPlotData
        __init__(self) -> FunctionPlotData

        Create empty data to be read from file. 
        
";

%feature("docstring")  dolfin::FunctionPlotData::vertex_values "

        vertex_values(self) -> GenericVector

        Return vertex values. 
        
";

%feature("docstring")  dolfin::GenericTensor "

    This class defines a common interface for arbitrary rank tensors.

    C++ includes: GenericTensor.h 
    
";

%feature("docstring")  dolfin::GenericTensor::init "

        init(self, GenericSparsityPattern sparsity_pattern)

        Initialize zero tensor using sparsity pattern. 
        
";

%feature("docstring")  dolfin::GenericTensor::factory "

        factory(self) -> LinearAlgebraFactory

        Return linear algebra backend factory. 
        
";

%feature("docstring")  dolfin::GenericTensor::rank "

        rank(self) -> uint

        Return tensor rank (number of dimensions). 
        
";

%feature("docstring")  dolfin::GenericTensor::zero "

        zero(self)

        Set all entries to zero and keep any sparse structure. 
        
";

%feature("docstring")  dolfin::GenericTensor::resize "

        resize(self, uint rank, uint dims)

        Resize tensor with given dimensions. 
        
";

%feature("docstring")  dolfin::GenericTensor::apply "

        apply(self, string mode)

        Finalize assembly of tensor. 
        
";

%feature("docstring")  dolfin::GenericTensor::copy "

        copy(self) -> GenericTensor

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::GenericTensor::GenericTensor "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::GenericTensor::size "

        size(self, uint dim) -> uint

        Return size of given dimension. 
        
";

%feature("docstring")  dolfin::UnitCube "

    Tetrahedral mesh of the 3D unit cube (0,1) x (0,1) x (0,1). Given the
    number of cells (nx, ny, nz) in each direction, the total number of
    tetrahedra will be 6*nx*ny*nz and the total number of vertices will be
    (nx + 1)*(ny + 1)*(nz + 1).

    C++ includes: UnitCube.h 
    
";

%feature("docstring")  dolfin::UnitCube::UnitCube "
__init__(self, uint nx, uint ny, uint nz) -> UnitCube
";

%feature("docstring")  dolfin::Function "

    This class represents a function u_h in a finite element function
    space V_h, given by

    u_h = sum_i U_i phi_i

    where {phi_i}_i is a basis for V_h, and U is a vector of expansion
    coefficients for u_h.

    C++ includes: Function.h 
    
";

%feature("docstring")  dolfin::Function::function_space "
Return the FunctionSpace
";

%feature("docstring")  dolfin::Function::_function_space "

        _function_space(self) -> __dummy_19__

        Return shared pointer to function space. 
        
";

%feature("docstring")  dolfin::Function::geometric_dimension "

        geometric_dimension(self) -> uint

        Return geometric dimension. 
        
";

%feature("docstring")  dolfin::Function::interpolate "

        interpolate(self, GenericFunction v)

        Interpolate function (possibly non-matching meshes). 
        
";

%feature("docstring")  dolfin::Function::_in "

        _in(self, FunctionSpace V) -> bool

        Check if function is a member of the given function space. 
        
";

%feature("docstring")  dolfin::Function::extrapolate "

        extrapolate(self, Function v)

        Extrapolate function (from a possibly lower-degree function space). 
        
";

%feature("docstring")  dolfin::Function::vector "

        vector(self) -> GenericVector
        vector(self) -> GenericVector

        Return vector of expansion coefficients (const version). 
        
";

%feature("docstring")  dolfin::Function::eval "

        eval(self, DoubleArray values, DoubleArray x)
        eval(self, DoubleArray values, DoubleArray x, Cell dolfin_cell, 
            cell ufc_cell)

        Evaluate function for given data. 
        
";

%feature("docstring")  dolfin::Function::_sub "
_sub(self, uint i) -> Function
";

%feature("docstring")  dolfin::Function::assign "

        assign(self, Function v) -> Function
        assign(self, Expression v) -> Function
        
";

%feature("docstring")  dolfin::Function::Function "

        __init__(self, __dummy_19__ V) -> Function
        __init__(self, __dummy_19__ V, boost::shared_ptr<(dolfin::GenericVector)> x) -> Function
        __init__(self, __dummy_19__ V, GenericVector x) -> Function
        __init__(self, __dummy_19__ V, string filename) -> Function
        __init__(self, Function v) -> Function
        __init__(self, Function v, uint i) -> Function

        Sub-function constructor with shallow copy of vector (used in Python
        interface) 
        
";

%feature("docstring")  dolfin::UmfpackLUSolver "

    This class implements the direct solution (LU factorization) of linear
    systems of the form Ax = b using UMFPACK
    (http://www.cise.ufl.edu/research/sparse/umfpack/) if installed.

    C++ includes: UmfpackLUSolver.h 
    
";

%feature("docstring")  dolfin::UmfpackLUSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::UmfpackLUSolver::solve "

        solve(self, GenericVector x, GenericVector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint

        Solve linear system. 
        
";

%feature("docstring")  dolfin::UmfpackLUSolver::UmfpackLUSolver "

        __init__(self) -> UmfpackLUSolver
        __init__(self, GenericMatrix A) -> UmfpackLUSolver
        __init__(self, boost::shared_ptr<(q(const).dolfin::GenericMatrix)> A) -> UmfpackLUSolver

        Constructor. 
        
";

%feature("docstring")  dolfin::FunctionSpace "

    This class represents a finite element function space defined by a
    mesh, a finite element, and a local-to-global mapping of the degrees
    of freedom (dofmap).

    C++ includes: FunctionSpace.h 
    
";

%feature("docstring")  dolfin::FunctionSpace::str "
Return a string representation of it self
";

%feature("docstring")  dolfin::FunctionSpace::extract_sub_space "

        extract_sub_space(self, std::vector<(dolfin::uint)> component) -> __dummy_18__

        Extract sub space for component. 
        
";

%feature("docstring")  dolfin::FunctionSpace::collapse_sub_space "

        collapse_sub_space(self, __dummy_12__ dofmap) -> __dummy_18__

        Return function space with a new dof map. 
        
";

%feature("docstring")  dolfin::FunctionSpace::component "

        component(self) -> dolfin::Array<(dolfin::uint)>

        Return component (relative to super space). 
        
";

%feature("docstring")  dolfin::FunctionSpace::interpolate "

        interpolate(self, GenericVector expansion_coefficients, GenericFunction v)

        Interpolate function v into function space, returning the vector of
        expansion coefficients 
        
";

%feature("docstring")  dolfin::FunctionSpace::mesh "

        mesh(self) -> Mesh

        Return mesh. 
        
";

%feature("docstring")  dolfin::FunctionSpace::FunctionSpace "

        __init__(self, __dummy_36__ mesh, __dummy_17__ element, __dummy_13__ dofmap) -> FunctionSpace
        __init__(self, __dummy_37__ mesh, __dummy_17__ element, __dummy_13__ dofmap) -> FunctionSpace
        __init__(self, FunctionSpace V) -> FunctionSpace

        Copy constructor. 
        
";

%feature("docstring")  dolfin::FunctionSpace::dim "

        dim(self) -> uint

        Return dimension of function space. 
        
";

%feature("docstring")  dolfin::FunctionSpace::sub "
sub(self, uint i) -> __dummy_18__
";

%feature("docstring")  dolfin::FunctionSpace::element "

        element(self) -> FiniteElement

        Return finite element. 
        
";

%feature("docstring")  dolfin::FunctionSpace::dofmap "

        dofmap(self) -> GenericDofMap

        Return dofmap. 
        
";

%feature("docstring")  dolfin::FunctionSpace::has_element "

        has_element(self, FiniteElement element) -> bool

        Check if function space has given element. 
        
";

%feature("docstring")  dolfin::FunctionSpace::print_dofmap "

        print_dofmap(self)

        Print dofmap (useful for debugging). 
        
";

%feature("docstring")  dolfin::FunctionSpace::assign "
assign(self, FunctionSpace V) -> FunctionSpace
";

%feature("docstring")  dolfin::FunctionSpace::has_cell "

        has_cell(self, Cell cell) -> bool

        Check if function space has given cell. 
        
";

%feature("docstring")  dolfin::faces "

    A FaceIterator is a MeshEntityIterator of topological dimension 2.

    C++ includes: Face.h 
    
";

%feature("docstring")  dolfin::faces::_dereference "
_dereference(self) -> Face
";

%feature("docstring")  dolfin::faces::faces "

        FaceIterator(Mesh mesh) -> faces
        __init__(self, MeshEntity entity) -> faces
        
";

%feature("docstring")  dolfin::RealParameter "

    Parameter with value type double.

    C++ includes: Parameter.h 
    
";

%feature("docstring")  dolfin::RealParameter::_assign "

        _assign(self, double value) -> RealParameter
        _assign(self, real value) -> RealParameter
        
";

%feature("docstring")  dolfin::RealParameter::RealParameter "

        __init__(self, string key, real value) -> RealParameter

        Create double-valued parameter. 
        
";

%feature("docstring")  dolfin::PETScKrylovSolver "
Proxy of C++ dolfin::PETScKrylovSolver class
";

%feature("docstring")  dolfin::PETScKrylovSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::PETScKrylovSolver::set_operators "

        set_operators(self, GenericMatrix A, GenericMatrix P)
        set_operators(self, PETScBaseMatrix A, PETScBaseMatrix P)

        Solve the operator (matrix) and preconditioner matrix. 
        
";

%feature("docstring")  dolfin::PETScKrylovSolver::ksp "
ksp(self) -> boost::shared_ptr<(KSP)>
";

%feature("docstring")  dolfin::PETScKrylovSolver::set_operator "

        set_operator(self, GenericMatrix A)
        set_operator(self, PETScBaseMatrix A)

        Solve the operator (matrix). 
        
";

%feature("docstring")  dolfin::PETScKrylovSolver::solve "

        solve(self, GenericVector x, GenericVector b) -> uint
        solve(self, PETScVector x, PETScVector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint
        solve(self, PETScBaseMatrix A, PETScVector x, PETScVector b) -> uint

        Solve linear system Ax = b. 
        
";

%feature("docstring")  dolfin::PETScKrylovSolver::PETScKrylovSolver "

        __init__(self, string method = "default", string pc_type = "default") -> PETScKrylovSolver
        __init__(self, string method = "default") -> PETScKrylovSolver
        __init__(self) -> PETScKrylovSolver
        __init__(self, string method, PETScPreconditioner preconditioner) -> PETScKrylovSolver
        __init__(self, string method, PETScUserPreconditioner preconditioner) -> PETScKrylovSolver
        __init__(self, boost::shared_ptr<(KSP)> ksp) -> PETScKrylovSolver
        
";

%feature("docstring")  dolfin::ODESolution "
Proxy of C++ dolfin::ODESolution class
";

%feature("docstring")  dolfin::ODESolution::begin "
begin(self) -> iterator
";

%feature("docstring")  dolfin::ODESolution::add_timeslab "
add_timeslab(self, real a, real b, real values)
";

%feature("docstring")  dolfin::ODESolution::nsize "
nsize(self) -> uint
";

%feature("docstring")  dolfin::ODESolution::get_weights "

        get_weights(self) -> real

        Get pointer to weights. 
        
";

%feature("docstring")  dolfin::ODESolution::eval "

        eval(self, real t, real y)

        Evaluate (interpolate) value of solution at given time. 
        
";

%feature("docstring")  dolfin::ODESolution::flush "

        flush(self)

        Make object ready for evaluating, set to read mode. 
        
";

%feature("docstring")  dolfin::ODESolution::end "
end(self) -> iterator
";

%feature("docstring")  dolfin::ODESolution::endtime "
endtime(self) -> real
";

%feature("docstring")  dolfin::ODESolution::ODESolution "

        __init__(self) -> ODESolution
        __init__(self, string filename, uint number_of_files = 1) -> ODESolution
        __init__(self, string filename) -> ODESolution
        
";

%feature("docstring")  dolfin::ODESolution::size "
size(self) -> uint
";

%feature("docstring")  dolfin::ODESolution::save_to_file "
save_to_file(self)
";

%feature("docstring")  dolfin::ODESolution::get_timeslab "

        get_timeslab(self, uint index) -> ODESolutionData

        Get timeslab (used when iterating). 
        
";

%feature("docstring")  dolfin::ODESolution::set_filename "
set_filename(self, string filename)
";

%feature("docstring")  dolfin::ODESolution::init "
init(self, uint N, Lagrange trial, real quad_weights)
";

%feature("docstring")  dolfin::ODESolution::str "
str(self, bool verbose) -> string
";

%feature("docstring")  dolfin::ParameterValue "

    Base class for parameters.

    C++ includes: Parameter.h 
    
";

%feature("docstring")  dolfin::ParameterValue::description "

        description(self) -> string

        Return parameter description. 
        
";

%feature("docstring")  dolfin::ParameterValue::_assign "

        _assign(self, int value) -> ParameterValue
        _assign(self, double value) -> ParameterValue
        _assign(self, real value) -> ParameterValue
        _assign(self, string value) -> ParameterValue
        _assign(self, char value) -> ParameterValue
        
";

%feature("docstring")  dolfin::ParameterValue::_get_real_range "

        _get_real_range(self)

        Get range for string-valued parameter. 
        
";

%feature("docstring")  dolfin::ParameterValue::key "

        key(self) -> string

        Return parameter key. 
        
";

%feature("docstring")  dolfin::ParameterValue::range_str "

        range_str(self) -> string

        Return range string. 
        
";

%feature("docstring")  dolfin::ParameterValue::value_str "

        value_str(self) -> string

        Return value string. 
        
";

%feature("docstring")  dolfin::ParameterValue::type_str "

        type_str(self) -> string

        Return value type string. 
        
";

%feature("docstring")  dolfin::ParameterValue::_get_string_range "

        _get_string_range(self)

        Get range for string-valued parameter. 
        
";

%feature("docstring")  dolfin::ParameterValue::data "
Missing docstring
";

%feature("docstring")  dolfin::ParameterValue::ParameterValue "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::ParameterValue::change_count "

        change_count(self) -> uint

        Return change count (number of times parameter has been changed). 
        
";

%feature("docstring")  dolfin::ParameterValue::get_real "

        get_real(self) -> real

        Get real value of parameter with (possibly) extended precision. 
        
";

%feature("docstring")  dolfin::ParameterValue::access_count "

        access_count(self) -> uint

        Return access count (number of times parameter has been accessed). 
        
";

%feature("docstring")  dolfin::ParameterValue::set_range "

        set_range(self, int min_value, int max_value)
        set_range(self, real min_value, real max_value)
        set_range(self, std::set<(std::string)> range)

        Set range for string-valued parameter. 
        
";

%feature("docstring")  dolfin::ParameterValue::warn_once "
Missing docstring
";

%feature("docstring")  dolfin::ParameterValue::get_range "
Missing docstring
";

%feature("docstring")  dolfin::ParameterValue::value "
Missing docstring
";

%feature("docstring")  dolfin::ParameterValue::check_key "
check_key(string key)
";

%feature("docstring")  dolfin::ParameterValue::_assign_bool "
_assign_bool(self, bool value) -> ParameterValue
";

%feature("docstring")  dolfin::ParameterValue::_get_int_range "

        _get_int_range(self)

        Get range for string-valued parameter. 
        
";

%feature("docstring")  dolfin::ParameterValue::str "

        str(self) -> string

        Return short string description. 
        
";

%feature("docstring")  dolfin::edges "

    An EdgeIterator is a MeshEntityIterator of topological dimension 1.

    C++ includes: Edge.h 
    
";

%feature("docstring")  dolfin::edges::_dereference "
_dereference(self) -> Edge
";

%feature("docstring")  dolfin::edges::edges "

        EdgeIterator(Mesh mesh) -> edges
        __init__(self, MeshEntity entity) -> edges
        
";

%feature("docstring")  dolfin::TimeSeries "

    This class stores a time series of objects to file(s) in a binary
    format which is efficient for reading and writing.

    When objects are retrieved, the object stored at the time closest to
    the given time will be used.

    A new time series will check if values have been stored to file before
    (for a series with the same name) and in that case reuse those values.
    If new values are stored, old values will be cleared.

    C++ includes: TimeSeries.h 
    
";

%feature("docstring")  dolfin::TimeSeries::store "

        store(self, GenericVector vector, double t)
        store(self, Mesh mesh, double t)

        Store mesh at given time. 
        
";

%feature("docstring")  dolfin::TimeSeries::retrieve "

        retrieve(self, GenericVector vector, double t)
        retrieve(self, Mesh mesh, double t)

        Retrieve mesh at given time. 
        
";

%feature("docstring")  dolfin::TimeSeries::vector_times "

        vector_times(self) -> DoubleArray

        Return array of sample times for vectors. 
        
";

%feature("docstring")  dolfin::TimeSeries::filename_data "
filename_data(string series_name, string type_name, uint index) -> string
";

%feature("docstring")  dolfin::TimeSeries::filename_times "
filename_times(string series_name, string type_name) -> string
";

%feature("docstring")  dolfin::TimeSeries::TimeSeries "

        __init__(self, string name) -> TimeSeries

        Create empty time series. 
        
";

%feature("docstring")  dolfin::TimeSeries::clear "

        clear(self)

        Clear time series. 
        
";

%feature("docstring")  dolfin::TimeSeries::mesh_times "

        mesh_times(self) -> DoubleArray

        Return array of sample times for meshes. 
        
";

%feature("docstring")  dolfin::GaussianQuadrature "

    Gaussian-type quadrature rule on the double line, including Gauss,
    Radau, and Lobatto quadrature.

    Points and weights are computed to be exact within a tolerance of
    DOLFIN_EPS. Comparing with known exact values for n <= 3 shows that we
    obtain full precision (16 digits, error less than 2e-16).

    C++ includes: GaussianQuadrature.h 
    
";

%feature("docstring")  dolfin::GaussianQuadrature::GaussianQuadrature "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::Legendre "

    Legendre polynomial of given degree n on the interval [-1,1].

    P0(x) = 1 P1(x) = x P2(x) = (3x^2 - 1) / 2 ...

    The function values and derivatives are computed using three-term
    recurrence formulas.

    C++ includes: Legendre.h 
    
";

%feature("docstring")  dolfin::Legendre::ddx "

        ddx(self, real x) -> real
        ddx(self, uint n, real x) -> real
        
";

%feature("docstring")  dolfin::Legendre::eval "

        eval(self, uint nn, real x) -> real

        Evaluation of arbitrary order, nn <= n (useful ie in RadauQuadrature).

        
";

%feature("docstring")  dolfin::Legendre::d2dx "

        d2dx(self, real x) -> real
        d2dx(self, uint n, real x) -> real
        
";

%feature("docstring")  dolfin::Legendre::Legendre "
__init__(self, uint n) -> Legendre
";

%feature("docstring")  dolfin::ITLKrylovSolver "
Proxy of C++ dolfin::ITLKrylovSolver class
";

%feature("docstring")  dolfin::ITLKrylovSolver::default_parameters "
default_parameters() -> Parameters
";

%feature("docstring")  dolfin::ITLKrylovSolver::solve "

        solve(self, GenericVector x, GenericVector b) -> uint
        solve(self, MTL4Vector x, MTL4Vector b) -> uint
        solve(self, GenericMatrix A, GenericVector x, GenericVector b) -> uint

        Solve linear system Ax = b. 
        
";

%feature("docstring")  dolfin::ITLKrylovSolver::ITLKrylovSolver "

        __init__(self, string method = "default", string pc_type = "default") -> ITLKrylovSolver
        __init__(self, string method = "default") -> ITLKrylovSolver
        __init__(self) -> ITLKrylovSolver
        
";

%feature("docstring")  dolfin::SubDomain "

    This class defines the interface for definition of sub domains.
    Alternatively, sub domains may be defined by a Mesh and a
    MeshFunction<uint> over the mesh.

    C++ includes: SubDomain.h 
    
";

%feature("docstring")  dolfin::SubDomain::map "

        map(self, DoubleArray x, DoubleArray arg0)

        Map coordinate x in domain H to coordinate y in domain G (used for
        periodic boundary conditions) 
        
";

%feature("docstring")  dolfin::SubDomain::inside "

        inside(self, DoubleArray x, bool on_boundary) -> bool

        Return true for points inside the subdomain. 
        
";

%feature("docstring")  dolfin::SubDomain::geometric_dimension "

        geometric_dimension(self) -> uint

        Return geometric dimension. 
        
";

%feature("docstring")  dolfin::SubDomain::mark "

        mark(self, MeshFunctionUInt sub_domains, uint sub_domain)

        Set sub domain markers for given subdomain. 
        
";

%feature("docstring")  dolfin::SubDomain::snap "

        snap(self, DoubleArray x)

        Snap coordinate to boundary of sub domain. 
        
";

%feature("docstring")  dolfin::SubDomain::SubDomain "

        __init__(self) -> SubDomain

        Constructor. 
        
";

%feature("docstring")  dolfin::ODESolutionData "
Proxy of C++ dolfin::ODESolutionData class
";

%feature("docstring")  dolfin::ODESolutionData::b "
b(self) -> real
";

%feature("docstring")  dolfin::ODESolutionData::eval_a "
eval_a(self, real u)
";

%feature("docstring")  dolfin::ODESolutionData::ODESolutionData "

        __init__(self, real a, real k, uint nodal_size, uint N, real values) -> ODESolutionData
        __init__(self, ODESolutionData cp) -> ODESolutionData
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix "

    This class provides a simple matrix class based on uBLAS. It is a
    simple wrapper for a uBLAS matrix implementing the GenericMatrix
    interface.

    The interface is intentionally simple. For advanced usage, access the
    underlying uBLAS matrix and use the standard uBLAS interface which is
    documented athttp://www.boost.org/libs/numeric/ublas/doc/index.htm.

    Developer note: specialised member functions must be inlined to avoid
    link errors.

    C++ includes: uBLASMatrix.h 
    
";

%feature("docstring")  dolfin::uBLASDenseMatrix::mat "

        mat(self) -> dolfin::ublas::matrix<(double)>
        mat(self) -> dolfin::ublas::matrix<(double)>

        Return reference to uBLAS matrix (non-const version). 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::invert "

        invert(self)

        Compute inverse of matrix. 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::compress "

        compress(self)

        Compress matrix (eliminate all non-zeros from a sparse matrix). 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::lump "

        lump(self, uBLASVector m)

        Lump matrix into vector m. 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::solveInPlace "

        solveInPlace(self, uBLASVector x, uBLASVector b)

        Solve Ax = b in-place using uBLAS(A is destroyed). 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::zero "

        zero(self)
        zero(self, uint m)

        Set given rows to zero. 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::solve "

        solve(self, uBLASVector x, uBLASVector b)

        Solve Ax = b out-of-place using uBLAS (A is not destroyed). 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::copy "

        copy(self) -> uBLASDenseMatrix

        Return copy of tensor. 
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::assign "

        assign(self, GenericMatrix A) -> GenericMatrix
        assign(self, uBLASDenseMatrix A) -> uBLASDenseMatrix
        
";

%feature("docstring")  dolfin::uBLASDenseMatrix::uBLASDenseMatrix "

        __init__(self) -> uBLASDenseMatrix
        __init__(self, uint M, uint N) -> uBLASDenseMatrix
        __init__(self, uBLASDenseMatrix A) -> uBLASDenseMatrix

        Create matrix from given uBLAS matrix expression. 
        
";

%feature("docstring")  dolfin::ODESolutionIterator "
Proxy of C++ dolfin::ODESolutionIterator class
";

%feature("docstring")  dolfin::ODESolutionIterator::get_index "
get_index(self) -> uint
";

%feature("docstring")  dolfin::ODESolutionIterator::get_ODESolution "
get_ODESolution(self) -> ODESolution
";

%feature("docstring")  dolfin::ODESolutionIterator::ODESolutionIterator "

        __init__(self, ODESolution u) -> ODESolutionIterator
        __init__(self, ODESolution u, int index) -> ODESolutionIterator
        __init__(self, ODESolutionIterator it) -> ODESolutionIterator
        
";

%feature("docstring")  dolfin::Assembler "

    This class provides automated assembly of linear systems, or more
    generally, assembly of a sparse tensor from a given variational form.

    The MeshFunction arguments can be used to specify assembly over
    subdomains of the mesh cells, exterior facets or interior facets.
    Either a null pointer or an empty MeshFunction may be used to specify
    that the tensor should be assembled over the entire set of cells or
    facets.

    C++ includes: Assembler.h 
    
";

%feature("docstring")  dolfin::Assembler::assemble "

        assemble(GenericTensor A, Form a, bool reset_sparsity = True, 
            bool add_values = False)
        assemble(GenericTensor A, Form a, bool reset_sparsity = True)
        assemble(GenericTensor A, Form a)
        assemble(GenericTensor A, Form a, SubDomain sub_domain, bool reset_sparsity = True, 
            bool add_values = False)
        assemble(GenericTensor A, Form a, SubDomain sub_domain, bool reset_sparsity = True)
        assemble(GenericTensor A, Form a, SubDomain sub_domain)
        assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains, bool reset_sparsity = True, 
            bool add_values = False)
        assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains, bool reset_sparsity = True)
        assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
            MeshFunctionUInt exterior_facet_domains, 
            MeshFunctionUInt interior_facet_domains)
        
";

%feature("docstring")  dolfin::Assembler::Assembler "
__init__(self) -> Assembler
";

%feature("docstring")  dolfin::CellSize "

    This Function represents the local cell size on a given mesh.

    C++ includes: SpecialFunctions.h 
    
";

%feature("docstring")  dolfin::CellSize::CellSize "

        __init__(self, Mesh mesh) -> CellSize

        Constructor. 
        
";

%feature("docstring")  dolfin::Sample "

    Sample of solution values at a given point.

    C++ includes: Sample.h 
    
";

%feature("docstring")  dolfin::Sample::k "

        k(self, uint index) -> real

        Return time step for component with given index. 
        
";

%feature("docstring")  dolfin::Sample::r "

        r(self, uint index) -> real

        Return residual for component with given index. 
        
";

%feature("docstring")  dolfin::Sample::u "

        u(self, uint index) -> real

        Return value of component with given index. 
        
";

%feature("docstring")  dolfin::Sample::t "

        t(self) -> real

        Return time t. 
        
";

%feature("docstring")  dolfin::Sample::Sample "

        __init__(self, TimeSlab timeslab, real t, string name, string label) -> Sample

        Constructor. 
        
";

%feature("docstring")  dolfin::Sample::size "

        size(self) -> uint

        Return number of components. 
        
";

%feature("docstring")  dolfin::File "

    A File represents a data file for reading and writing objects. Unless
    specified explicitly, the format is determined by the file name
    suffix. A list of objects that can be read/written to file can be
    found in GenericFile.h

    C++ includes: File.h 
    
";

%feature("docstring")  dolfin::File::exists "
exists(string filename) -> bool
";

%feature("docstring")  dolfin::File::File "

        __init__(self, string filename, string encoding = "ascii") -> File
        __init__(self, string filename) -> File
        __init__(self, string filename, Type type, string encoding = "ascii") -> File
        __init__(self, string filename, Type type) -> File
        __init__(self, std::ostream outstream) -> File

        Create a outfile object writing to stream. 
        
";

%feature("docstring")  dolfin::Table "

    This class provides storage and pretty-printing for tables. Example
    usage:

    Table table("Timings");

    table("uBLAS", "Assemble") = 0.010; table("uBLAS", "Solve") =
    0.020; table("PETSc", "Assemble") = 0.011; table("PETSc",
    "Solve") = 0.019; table("Epetra", "Assemble") = 0.012;
    table("Epetra", "Solve") = 0.018;

    info(table);

    C++ includes: Table.h 
    
";

%feature("docstring")  dolfin::Table::set "

        set(self, string row, string col, int value)
        set(self, string row, string col, double value)
        set(self, string row, string col, string value)

        Set value of table entry. 
        
";

%feature("docstring")  dolfin::Table::get "

        get(self, string row, string col) -> string

        Get value of table entry. 
        
";

%feature("docstring")  dolfin::Table::title "

        title(self) -> string

        Return table title. 
        
";

%feature("docstring")  dolfin::Table::str_latex "

        str_latex(self) -> string

        Return informal string representation for LaTeX. 
        
";

%feature("docstring")  dolfin::Table::Table "

        __init__(self, string title = "") -> Table
        __init__(self) -> Table

        Create empty table. 
        
";

%feature("docstring")  dolfin::MeshData "

    The class MeshData is a container for auxiliary mesh data, represented
    either as MeshFunctions over topological mesh entities, arrays or
    maps. Each dataset is identified by a unique user-specified string.
    Only uint- valued data are currently supported.

    The following named mesh data are recognized by DOLFIN:

    Boundary indicators

    "boundary facet cells" - Array<uint> of size num_facets "boundary
    facet numbers" - Array<uint> of size num_facets "boundary
    indicators" - Array<uint> of size num_facets "material indicators"
    - MeshFunction<uint> of dimension D

    Boundary indicators (alternative)

    "exterior facet domains" - MeshFunction<uint> of dimension D - 1

    Facet orientation (used for assembly over interior facets)

    "facet orientation" - MeshFunction<uint> of dimension D - 1

    Boundary extraction

    "vertex map" - MeshFunction<uint> of dimension 0 "cell map" -
    MeshFunction<uint> of dimension D

    Mesh partitioning

    "global entity indices %d" - MeshFunction<uint> of dimension 0, 1,
    ..., D "exterior facets" - MeshFunction<uint> of dimension D - 1
    "num global entities" - Array<uint> of size D + 1 "overlap" -
    vector mapping

    Sub meshes

    "global vertex indices" - MeshFunction<uint> of dimension 0

    C++ includes: MeshData.h 
    
";

%feature("docstring")  dolfin::MeshData::vector_mapping "

        vector_mapping(self, string name) -> std::map<(dolfin::uint,std::vector<(dolfin::uint)>)>

        Return vector mapping with given name (returning zero if data is not
        available). 
        
";

%feature("docstring")  dolfin::MeshData::erase_mesh_function "

        erase_mesh_function(self, string name)

        Erase MeshFunction with given name. 
        
";

%feature("docstring")  dolfin::MeshData::erase_array "

        erase_array(self, string name)

        Erase array with given name. 
        
";

%feature("docstring")  dolfin::MeshData::mapping "

        mapping(self, string name) -> std::map<(dolfin::uint,dolfin::uint)>

        Return mapping with given name (returning zero if data is not
        available). 
        
";

%feature("docstring")  dolfin::MeshData::create_mapping "

        create_mapping(self, string name) -> std::map<(dolfin::uint,dolfin::uint)>

        Create mapping from uint to uint with given name. 
        
";

%feature("docstring")  dolfin::MeshData::create_mesh_function "

        create_mesh_function(self, string name) -> MeshFunctionUInt
        create_mesh_function(self, string name, uint dim) -> MeshFunctionUInt

        Create MeshFunction with given name and dimension. 
        
";

%feature("docstring")  dolfin::MeshData::create_vector_mapping "

        create_vector_mapping(self, string name) -> std::map<(dolfin::uint,std::vector<(dolfin::uint)>)>

        Create mapping from uint to vector of uint with given name. 
        
";

%feature("docstring")  dolfin::MeshData::erase_mapping "

        erase_mapping(self, string name)

        Erase mapping with given name. 
        
";

%feature("docstring")  dolfin::MeshData::array "

        array(self, string name) -> std::vector<(dolfin::uint)>

        Return array with given name (returning zero if data is not
        available). 
        
";

%feature("docstring")  dolfin::MeshData::MeshData "

        __init__(self, Mesh mesh) -> MeshData

        Constructor. 
        
";

%feature("docstring")  dolfin::MeshData::erase_vector_mapping "

        erase_vector_mapping(self, string name)

        Erase vector mapping with given name. 
        
";

%feature("docstring")  dolfin::MeshData::create_array "

        create_array(self, string name, uint size) -> std::vector<(dolfin::uint)>

        Create array (vector) with given name and size. 
        
";

%feature("docstring")  dolfin::MeshData::clear "

        clear(self)

        Clear all data. 
        
";

%feature("docstring")  dolfin::MeshData::mesh_function "

        mesh_function(self, string name) -> MeshFunctionUInt

        Return MeshFunction with given name (returning zero if data is not
        available). 
        
";

%feature("docstring")  dolfin::uBLASPreconditioner "

    This class specifies the interface for preconditioners for the uBLAS
    Krylov solver.

    C++ includes: uBLASPreconditioner.h 
    
";

%feature("docstring")  dolfin::uBLASPreconditioner::init "

        init(self, uBLASSparseMatrix P)
        init(self, uBLASDenseMatrix P)
        init(self, uBLASKrylovMatrix P)

        Initialise preconditioner (virtual matrix). 
        
";

%feature("docstring")  dolfin::uBLASPreconditioner::solve "

        solve(self, uBLASVector x, uBLASVector b)

        Solve linear system (M^-1)Ax = y. 
        
";

%feature("docstring")  dolfin::uBLASPreconditioner::uBLASPreconditioner "
No constructor defined - class is abstract
";

%feature("docstring")  dolfin::MeshPartitioning_partition "

    partition(Mesh mesh)
    MeshPartitioning_partition(Mesh mesh, LocalMeshData data)
    
";

%feature("docstring")  dolfin::down_cast_PETScMatrix "
down_cast_PETScMatrix(GenericTensor tensor) -> PETScMatrix
";

%feature("docstring")  dolfin::has_type_uBLASVector "
has_type_uBLASVector(GenericTensor tensor) -> bool
";

%feature("docstring")  dolfin::ParameterValue_check_key "
ParameterValue_check_key(string key)
";

%feature("docstring")  dolfin::PETScLUSolver_default_parameters "
PETScLUSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::info_underline "

    info_underline(string msg, v(...) *args)

    Print underlined message. 
    
";

%feature("docstring")  dolfin::has_mpi "
has_mpi() -> bool
";

%feature("docstring")  dolfin::down_cast "
Cast tensor to the given subclass, passing the wrong class is an error.
";

%feature("docstring")  dolfin::has_scotch "
has_scotch() -> bool
";

%feature("docstring")  dolfin::has_slepc "
has_slepc() -> bool
";

%feature("docstring")  dolfin::ITLKrylovSolver_default_parameters "
ITLKrylovSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::isnormal "
isnormal(real x) -> int
";

%feature("docstring")  dolfin::has_la_backend "
has_la_backend(string backend) -> bool
";

%feature("docstring")  dolfin::real_max "
real_max(real x, real y) -> real
";

%feature("docstring")  dolfin::get_log_level "

    get_log_level() -> int

    Get log level. 
    
";

%feature("docstring")  dolfin::_get_matrix_single_item "
_get_matrix_single_item(GenericMatrix self, int m, int n) -> double
";

%feature("docstring")  dolfin::SystemAssembler_assemble "

    assemble(GenericMatrix A, GenericVector b, Form a, Form L, bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, bool reset_sparsity = True)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc, 
        bool reset_sparsity = True, bool add_values = True)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc, 
        bool reset_sparsity = True)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        bool reset_sparsity = True)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        GenericVector x0, 
        bool reset_sparsity = True, bool add_values = False)
    assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        GenericVector x0, 
        bool reset_sparsity = True)
    SystemAssembler_assemble(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        GenericVector x0)
    
";

%feature("docstring")  dolfin::real_mat_pow "
real_mat_pow(uint n, real A, real B, uint q)
";

%feature("docstring")  dolfin::uBLASDenseFactory_instance "
uBLASDenseFactory_instance() -> uBLASDenseFactory
";

%feature("docstring")  dolfin::PrimitiveIntersector_do_intersect "

    do_intersect(MeshEntity entity_1, MeshEntity entity_2) -> bool
    PrimitiveIntersector_do_intersect(MeshEntity entity_1, Point point) -> bool
    
";

%feature("docstring")  dolfin::has_type "
Return wether tensor is of the given subclass.
";

%feature("docstring")  dolfin::MPI_is_broadcaster "
MPI_is_broadcaster() -> bool
";

%feature("docstring")  dolfin::UnitInterval_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::CholmodCholeskySolver_default_parameters "
CholmodCholeskySolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::CellSize_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::TimeSeries_filename_times "
TimeSeries_filename_times(string series_name, string type_name) -> string
";

%feature("docstring")  dolfin::uBLASSparseFactory_instance "
uBLASSparseFactory_instance() -> uBLASSparseFactory
";

%feature("docstring")  dolfin::MPI_local_range "

    local_range(uint N) -> std::pair<(dolfin::uint,dolfin::uint)>
    MPI_local_range(uint process, uint N) -> std::pair<(dolfin::uint,dolfin::uint)>
    
";

%feature("docstring")  dolfin::dolfin_swigversion "
dolfin_swigversion() -> int
";

%feature("docstring")  dolfin::BoundaryMesh_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::get_tensor_type "
Return the concrete subclass of tensor.
";

%feature("docstring")  dolfin::down_cast_MTL4Vector "
down_cast_MTL4Vector(GenericTensor tensor) -> MTL4Vector
";

%feature("docstring")  dolfin::has_type_uBLASDenseMatrix "
has_type_uBLASDenseMatrix(GenericTensor tensor) -> bool
";

%feature("docstring")  dolfin::real_mat_prod_inplace "
real_mat_prod_inplace(uint n, real A, real B)
";

%feature("docstring")  dolfin::timing "

    timing(string task, bool reset = False) -> double
    timing(string task) -> double

    Return timing (average) for given task, optionally clearing timing for
    task. 
    
";

%feature("docstring")  dolfin::real_min "
real_min(real x, real y) -> real
";

%feature("docstring")  dolfin::to_double "
to_double(real x) -> double
";

%feature("docstring")  dolfin::has_cgal "
has_cgal() -> bool
";

%feature("docstring")  dolfin::real_log "

    real_log(real x) -> real

    Logarithmic function (note: not full precision!). 
    
";

%feature("docstring")  dolfin::KrylovSolver_default_parameters "
KrylovSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::MeshCoordinates_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::has_type_uBLASSparseMatrix "
has_type_uBLASSparseMatrix(GenericTensor tensor) -> bool
";

%feature("docstring")  dolfin::summary "

    summary(bool reset = False)
    summary()

    Print summary of timings and tasks, optionally clearing stored
    timings. 
    
";

%feature("docstring")  dolfin::has_type_PETScMatrix "
has_type_PETScMatrix(GenericTensor tensor) -> bool
";

%feature("docstring")  dolfin::has_type_PETScVector "
has_type_PETScVector(GenericTensor tensor) -> bool
";

%feature("docstring")  dolfin::debug "
Missing docstring
";

%feature("docstring")  dolfin::CellType_string2type "
CellType_string2type(string type) -> Type
";

%feature("docstring")  dolfin::_get_vector_values "
_get_vector_values(GenericVector self) -> DoubleArray
";

%feature("docstring")  dolfin::_set_vector_items_value "
_set_vector_items_value(GenericVector self, PyObject op, double value)
";

%feature("docstring")  dolfin::down_cast_uBLASVector "
down_cast_uBLASVector(GenericTensor tensor) -> uBLASVector
";

%feature("docstring")  dolfin::_set_matrix_items_matrix "
_set_matrix_items_matrix(GenericMatrix self, GenericMatrix arg1)
";

%feature("docstring")  dolfin::GlobalParameters_default_parameters "
GlobalParameters_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::MPI_index_owner "
MPI_index_owner(uint index, uint N) -> uint
";

%feature("docstring")  dolfin::Box_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::Constant_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::assemble_system "

    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, bool reset_sparsitys = True, 
        bool add_values = False)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, bool reset_sparsitys = True)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc, 
        bool reset_sparsitys = True, bool add_values = False)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc, 
        bool reset_sparsitys = True)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, DirichletBC bc)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        bool reset_sparsitys = True, 
        bool add_values = False)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        bool reset_sparsitys = True)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        GenericVector x0, 
        bool reset_sparsitys = True, bool add_values = False)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        GenericVector x0, 
        bool reset_sparsitys = True)
    assemble_system(GenericMatrix A, GenericVector b, Form a, Form L, std::vector<(p.q(const).dolfin::DirichletBC)> bcs, 
        MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        GenericVector x0)

    Assemble system (A, b) on sub domains and apply Dirichlet boundary
    conditions. 
    
";

%feature("docstring")  dolfin::MPI_is_receiver "
MPI_is_receiver() -> bool
";

%feature("docstring")  dolfin::LUSolver_default_parameters "
LUSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::_get_vector_single_item "
_get_vector_single_item(GenericVector self, int index) -> double
";

%feature("docstring")  dolfin::TimeSeries_filename_data "
TimeSeries_filename_data(string series_name, string type_name, uint index) -> string
";

%feature("docstring")  dolfin::has_gmp "
has_gmp() -> bool
";

%feature("docstring")  dolfin::DofMap_extract_sub_dofmap "

    extract_sub_dofmap(std::vector<(dolfin::uint)> component, Mesh dolfin_mesh) -> DofMap
    DofMap_extract_sub_dofmap(dof_map ufc_dof_map, uint offset, std::vector<(dolfin::uint)> component, 
        mesh ufc_mesh, Mesh dolfin_mesh) -> dof_map

    Extract sub dofmap component. 
    
";

%feature("docstring")  dolfin::check_equal "

    check_equal(uint value, uint valid_value, string task, string value_name)

    Check value and print an informative error message if invalid. 
    
";

%feature("docstring")  dolfin::MPI_global_maximum "
MPI_global_maximum(uint size) -> uint
";

%feature("docstring")  dolfin::real_mat_prod "
real_mat_prod(uint n, real res, real A, real B)
";

%feature("docstring")  dolfin::has_zlib "
has_zlib() -> bool
";

%feature("docstring")  dolfin::SubMesh_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::STLFactory_instance "
STLFactory_instance() -> STLFactory
";

%feature("docstring")  dolfin::residual "

    residual(GenericMatrix A, GenericVector x, GenericVector b) -> double

    Compute residual ||Ax - b||. 
    
";

%feature("docstring")  dolfin::Interval_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::UnitCircle_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::assemble "

    assemble(GenericTensor A, Form a, bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericTensor A, Form a, bool reset_sparsity = True)
    assemble(GenericTensor A, Form a)
    assemble(GenericTensor A, Form a, SubDomain sub_domain, bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericTensor A, Form a, SubDomain sub_domain, bool reset_sparsity = True)
    assemble(GenericTensor A, Form a, SubDomain sub_domain)
    assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, bool reset_sparsity = True)
    assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains)
    assemble(Form a, bool reset_sparsity = True, bool add_values = False) -> double
    assemble(Form a, bool reset_sparsity = True) -> double
    assemble(Form a) -> double
    assemble(Form a, SubDomain sub_domain, bool reset_sparsity = True, 
        bool add_values = False) -> double
    assemble(Form a, SubDomain sub_domain, bool reset_sparsity = True) -> double
    assemble(Form a, SubDomain sub_domain) -> double
    assemble(Form a, MeshFunctionUInt cell_domains, MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        bool reset_sparsity = True, 
        bool add_values = False) -> double
    assemble(Form a, MeshFunctionUInt cell_domains, MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, 
        bool reset_sparsity = True) -> double
    assemble(Form a, MeshFunctionUInt cell_domains, MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains) -> double

    Assemble scalar on sub domains. 
    
";

%feature("docstring")  dolfin::_set_matrix_items_array_of_float "
_set_matrix_items_array_of_float(GenericMatrix self, PyObject op, PyObject other)
";

%feature("docstring")  dolfin::info_stream "

    info_stream(std::ostream out, string msg)

    Print message to stream. 
    
";

%feature("docstring")  dolfin::real_pow "

    real_pow(real x, uint y) -> real
    real_pow(real x, real y) -> real
    
";

%feature("docstring")  dolfin::logging "

    logging(bool active = True)
    logging()

    Turn logging on or off. 
    
";

%feature("docstring")  dolfin::UmfpackLUSolver_default_parameters "
UmfpackLUSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::DofMap_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::PrimitiveIntersector_do_intersect_exact "

    do_intersect_exact(MeshEntity entity_1, MeshEntity entity_2) -> bool
    PrimitiveIntersector_do_intersect_exact(MeshEntity entity_1, Point point) -> bool
    
";

%feature("docstring")  dolfin::_info "

    _info(string msg, v(...) *args)
    _info(int debug_level, string msg, v(...) *args)

    Print variable (using output of str() method). 
    
";

%feature("docstring")  dolfin::MeshPartitioning_number_entities "
MeshPartitioning_number_entities(Mesh mesh, uint d)
";

%feature("docstring")  dolfin::SLEPcEigenSolver_default_parameters "
SLEPcEigenSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::rand "

    rand() -> double

    Return a random number, uniformly distributed between [0.0, 1.0). 
    
";

%feature("docstring")  dolfin::DirichletBC_default_parameters "
DirichletBC_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::pow "
pow(real x, real y) -> real
";

%feature("docstring")  dolfin::UnitSphere_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::_get_matrix_sub_vector "
_get_matrix_sub_vector(GenericMatrix self, uint single, PyObject op, bool row) -> GenericVector
";

%feature("docstring")  dolfin::File_exists "
File_exists(string filename) -> bool
";

%feature("docstring")  dolfin::MPI_scatter "

    scatter(std::vector<(dolfin::uint)> values, uint sending_process = 0)
    scatter(std::vector<(dolfin::uint)> values)
    scatter(std::vector<(std::vector<(dolfin::uint)>)> values, 
        uint sending_process = 0)
    scatter(std::vector<(std::vector<(dolfin::uint)>)> values)
    scatter(std::vector<(std::vector<(double)>)> values, uint sending_process = 0)
    MPI_scatter(std::vector<(std::vector<(double)>)> values)
    
";

%feature("docstring")  dolfin::to_real "
to_real(double x) -> real
";

%feature("docstring")  dolfin::Expression_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::real_decimal_prec "
real_decimal_prec() -> int
";

%feature("docstring")  dolfin::_contains "
_contains(GenericVector self, double value) -> bool
";

%feature("docstring")  dolfin::_refine "

    _refine(Mesh mesh) -> Mesh
    _refine(Mesh refined_mesh, Mesh mesh)
    _refine(Mesh mesh, MeshFunctionBool cell_markers) -> Mesh
    _refine(Mesh refined_mesh, Mesh mesh, MeshFunctionBool cell_markers)

    Create locally refined mesh. 
    
";

%feature("docstring")  dolfin::real_pi "

    real_pi() -> real

    Compute pi. 
    
";

%feature("docstring")  dolfin::toc "

    toc() -> double

    Return elapsed CPU time (should not be used internally in DOLFIN!). 
    
";

%feature("docstring")  dolfin::VariationalProblem_default_parameters "
VariationalProblem_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::has_cholmod "
has_cholmod() -> bool
";

%feature("docstring")  dolfin::UnitSquare_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::Rectangle_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::ipow "

    ipow(uint a, uint n) -> uint

    Return a to the power n. 
    
";

%feature("docstring")  dolfin::down_cast_uBLASDenseMatrix "
down_cast_uBLASDenseMatrix(GenericTensor tensor) -> uBLASDenseMatrix
";

%feature("docstring")  dolfin::info "
Missing docstring
";

%feature("docstring")  dolfin::ALE_move "

    move(Mesh mesh, BoundaryMesh new_boundary, ALEType method = lagrange)
    move(Mesh mesh, BoundaryMesh new_boundary)
    move(Mesh mesh0, Mesh mesh1, ALEType method = lagrange)
    move(Mesh mesh0, Mesh mesh1)
    ALE_move(Mesh mesh, Function displacement)
    
";

%feature("docstring")  dolfin::ODE_default_parameters "
ODE_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::UnitCube_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::_compare_vector_with_vector "
_compare_vector_with_vector(GenericVector self, GenericVector other, DolfinCompareType cmp_type) -> PyObject
";

%feature("docstring")  dolfin::seed "

    seed(unsigned int s)

    Seed random number generator. 
    
";

%feature("docstring")  dolfin::MPI_send_recv "

    send_recv(uint send_buffer, uint send_size, uint dest, uint recv_buffer, 
        uint recv_size, uint source) -> uint
    MPI_send_recv(double send_buffer, uint send_size, uint dest, double recv_buffer, 
        uint recv_size, uint source) -> uint
    
";

%feature("docstring")  dolfin::error "

    error(string msg, v(...) *args)

    Print error message and throw an exception. 
    
";

%feature("docstring")  dolfin::MPI_distribute "

    distribute(std::vector<(dolfin::uint)> values, std::vector<(dolfin::uint)> partition)
    MPI_distribute(std::vector<(dolfin::uint)> partition)
    
";

%feature("docstring")  dolfin::not_working_in_parallel "

    not_working_in_parallel(string what)

    Report that functionality has not (yet) been implemented to work in
    parallel. 
    
";

%feature("docstring")  dolfin::down_cast_uBLASSparseMatrix "
down_cast_uBLASSparseMatrix(GenericTensor tensor) -> uBLASSparseMatrix
";

%feature("docstring")  dolfin::Assembler_assemble "

    assemble(GenericTensor A, Form a, bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericTensor A, Form a, bool reset_sparsity = True)
    assemble(GenericTensor A, Form a)
    assemble(GenericTensor A, Form a, SubDomain sub_domain, bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericTensor A, Form a, SubDomain sub_domain, bool reset_sparsity = True)
    assemble(GenericTensor A, Form a, SubDomain sub_domain)
    assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, bool reset_sparsity = True, 
        bool add_values = False)
    assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains, bool reset_sparsity = True)
    Assembler_assemble(GenericTensor A, Form a, MeshFunctionUInt cell_domains, 
        MeshFunctionUInt exterior_facet_domains, 
        MeshFunctionUInt interior_facet_domains)
    
";

%feature("docstring")  dolfin::dolfin_set_precision "
dolfin_set_precision(uint prec)
";

%feature("docstring")  dolfin::real_abs "
real_abs(real x) -> real
";

%feature("docstring")  dolfin::PETScPreconditioner_default_parameters "
PETScPreconditioner_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::PETScFactory_instance "
PETScFactory_instance() -> PETScFactory
";

%feature("docstring")  dolfin::MTL4Factory_instance "
MTL4Factory_instance() -> MTL4Factory
";

%feature("docstring")  dolfin::FacetArea_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::warning "

    warning(string msg, v(...) *args)

    Print warning. 
    
";

%feature("docstring")  dolfin::PETScKrylovSolver_default_parameters "
PETScKrylovSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::NewtonSolver_default_parameters "
NewtonSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::dolfin_init "

    dolfin_init(int argc, char argv)

    Initialize DOLFIN (and PETSc) with command-line arguments. This should
    not be needed in most cases since the initialization is otherwise
    handled automatically. 
    
";

%feature("docstring")  dolfin::MPI_barrier "
MPI_barrier()
";

%feature("docstring")  dolfin::normalize "

    normalize(GenericVector x, string normalization_type = "average") -> double
    normalize(GenericVector x) -> double

    Normalize vector according to given normalization type. 
    
";

%feature("docstring")  dolfin::CellType_type2string "
CellType_type2string(Type type) -> string
";

%feature("docstring")  dolfin::end "

    end()

    End task (decrease indentation level). 
    
";

%feature("docstring")  dolfin::LinearSolver_default_parameters "
LinearSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::has_parmetis "
has_parmetis() -> bool
";

%feature("docstring")  dolfin::set_log_level "

    set_log_level(int level)

    Set log level. 
    
";

%feature("docstring")  dolfin::DomainBoundary_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::PETScUserPreconditioner_setup "
PETScUserPreconditioner_setup(KSP ksp, PETScUserPreconditioner pc)
";

%feature("docstring")  dolfin::CellType_create "

    create(Type type) -> CellType
    CellType_create(string type) -> CellType
    
";

%feature("docstring")  dolfin::real_frexp "
real_frexp(int exp, real x) -> double
";

%feature("docstring")  dolfin::real_mat_exp "

    real_mat_exp(uint n, real res, real A, uint p = 6)
    real_mat_exp(uint n, real res, real A)

    Compute matrix exponential using Pade approximation og degree p. 
    
";

%feature("docstring")  dolfin::MPI_gather "

    gather(uint value) -> std::vector<(dolfin::uint)>
    gather(std::vector<(dolfin::uint)> values)
    MPI_gather()
    
";

%feature("docstring")  dolfin::tic "

    tic()

    Start timing (should not be used internally in DOLFIN!).

    Timing functions measure CPU time as determined by clock(), the
    precision of which seems to be 0.01 seconds. 
    
";

%feature("docstring")  dolfin::real_sqrt "

    real_sqrt(real a) -> real

    Square root. 
    
";

%feature("docstring")  dolfin::begin "

    begin(string msg, v(...) *args)
    begin(int debug_level, string msg, v(...) *args)

    Begin task (increase indentation level). 
    
";

%feature("docstring")  dolfin::_set_matrix_single_item "
_set_matrix_single_item(GenericMatrix self, int m, int n, double value)
";

%feature("docstring")  dolfin::uBLASKrylovSolver_default_parameters "
uBLASKrylovSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::_set_vector_items_vector "
_set_vector_items_vector(GenericVector self, PyObject op, GenericVector other)
";

%feature("docstring")  dolfin::sqr "

    sqr(double x) -> double

    Return the square of x. 
    
";

%feature("docstring")  dolfin::has_type_MTL4Matrix "
has_type_MTL4Matrix(GenericTensor tensor) -> bool
";

%feature("docstring")  dolfin::has_type_MTL4Vector "
has_type_MTL4Vector(GenericTensor tensor) -> bool
";

%feature("docstring")  dolfin::has_umfpack "
has_umfpack() -> bool
";

%feature("docstring")  dolfin::MPI_process_number "
MPI_process_number() -> uint
";

%feature("docstring")  dolfin::_compare_vector_with_value "
_compare_vector_with_value(GenericVector self, double value, DolfinCompareType cmp_type) -> PyObject
";

%feature("docstring")  dolfin::_get_vector_sub_vector "
_get_vector_sub_vector(GenericVector self, PyObject op) -> GenericVector
";

%feature("docstring")  dolfin::MPI_global_offset "
MPI_global_offset(uint range, bool exclusive) -> uint
";

%feature("docstring")  dolfin::MPI_sum "

    sum(double value) -> double
    MPI_sum(uint value) -> uint
    
";

%feature("docstring")  dolfin::real_epsilon "
real_epsilon() -> real
";

%feature("docstring")  dolfin::_set_vector_items_array_of_float "
_set_vector_items_array_of_float(GenericVector self, PyObject op, PyObject other)
";

%feature("docstring")  dolfin::down_cast_PETScVector "
down_cast_PETScVector(GenericTensor tensor) -> PETScVector
";

%feature("docstring")  dolfin::real_mat_vector_prod "
real_mat_vector_prod(uint n, real y, real A, real x)
";

%feature("docstring")  dolfin::_set_matrix_items_vector "
_set_matrix_items_vector(GenericMatrix self, PyObject op, GenericVector other)
";

%feature("docstring")  dolfin::down_cast_MTL4Matrix "
down_cast_MTL4Matrix(GenericTensor tensor) -> MTL4Matrix
";

%feature("docstring")  dolfin::SubSpace_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::solve "

    solve(GenericMatrix A, GenericVector x, GenericVector b, 
        string solver_type = "lu", string pc_type = "default")
    solve(GenericMatrix A, GenericVector x, GenericVector b, 
        string solver_type = "lu")
    solve(GenericMatrix A, GenericVector x, GenericVector b)

    Solve linear system Ax = b. 
    
";

%feature("docstring")  dolfin::Function_SWIGSharedPtrUpcast "
Missing docstring
";

%feature("docstring")  dolfin::time "

    time() -> double

    Return current CPU time used by process. 
    
";

%feature("docstring")  dolfin::SingularSolver_default_parameters "
SingularSolver_default_parameters() -> Parameters
";

%feature("docstring")  dolfin::MPI_num_processes "
MPI_num_processes() -> uint
";

%feature("docstring")  dolfin::real_exp "

    real_exp(real x) -> real

    Exponential function (note: not full precision!). 
    
";

%feature("docstring")  dolfin::_get_matrix_sub_matrix "
_get_matrix_sub_matrix(GenericMatrix self, PyObject row_op, PyObject col_op) -> GenericMatrix
";

