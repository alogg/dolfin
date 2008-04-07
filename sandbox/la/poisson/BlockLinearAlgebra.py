#
# objdoc: Block Linear Algebra module
# Ola Skavhaug
# hacks to port to uBLAS by Kent-Andre Mardal (not completed yet)

"""
Python block matrices and block vectors. This module provides three classes for block linear algebra systems in
Python. The main class, C{ BlockMatrix} implements a dense block matrix, 
C{DiagBlockMatrix} implements a block diagonal matrix (used e.g. for
preconditioning, and C{BlockVector} is a simple block vector consiting of
numpy.array blocks.
"""

import operator, numpy, math

class ZeroBlock(object):
    """
    This class is a zeros block matrix used to fill upp empty blocks
    in BlockMatrix.
    """

    def __init__(self, *args):
        """Input arguments are two integer values
        specifying the size of the block"""
        if isinstance(args[0], int) and isinstance(args[1], int):
            n = args[0]
            m = args[1]
            self.n = n
            self.m = m
        else:
            print "Wrong input value in ZeroBlock"


    def shape(self):
        """Return the number of block rows and columns."""
        return (self.n, self.m)

    shape = property(shape)

    def __mul__(self, other):
        """Matrix vector product."""
        return 0.0*other


class BlockMatrix(object):
    """A simple block matrix implementation in Python. The blocks are typically
    matrices implemented in compile languages and later exposed to Python
    either manually or using some sort of wrapper generator software."""

    def __init__(self, *args):
        """Input arguments are either two integer values
        specifying the number of block rows and columns, or a nested Python
        list containing the block elements directly."""
        if isinstance(args[0], int):
            n,m = args
            self.n = n
            self.m = m
            self.data = []
            self.rdim = self.cdim = None
            for i in range(n):
                self.data.append([])
                for j in range(m):
                    self.data[i].append(0)
        elif isinstance(args[0], (list, tuple)):
            self.n = len(args)
            self.m = len(args[0])
            self.data = list(args)
            self.check()
        else:
            print "ERROR", args

    def check(self):
        m,n = self.shape
        rdims = numpy.zeros((m,n), dtype='l')
        cdims = numpy.zeros((m,n), dtype='l')
        for i in xrange(m):
            for j in xrange(n):
                mat = self.data[i][j]
                if not isinstance(mat, int):
                    _m, _n = mat.size(0), mat.size(1)
                    rdims[i,j] = _m
                    cdims[i,j] = _n

        rdim = numpy.zeros(m, dtype="l")
        cdim = numpy.zeros(n, dtype="l")
        for i in xrange(m):
            row = rdims[i]
            row = numpy.where(row==0, max(row), row)
            if min(row) == max(row):
                rdim[i] = max(row)
            else:
                raise RuntimeError, "Wrong row dimensions in block matrix"
        for j in xrange(n):
            col = cdims[:,j]
            col = numpy.where(col==0, max(col), col)
            if min(col) == max(col):
                cdim[j] = max(col)
            else:
                raise RuntimeError, "Wrong col dimensions in block matrix"
        self.rdim = rdim
        self.cdim = cdim
        for i in xrange(m):
            for j in xrange(n):
                if isinstance(self.data[i][j], int):
                    import Misc
                    self.data[i][j] = Misc.genZeroBlock(self.rdim[i])
                    self.data[i][j].n = self.cdim[j]

            
    def row_dim(self, i):
        return self.rdim[i]

    def col_dim(self, j):
        return self.cdim[j]

    def shape(self):
        """Return the number of block rows and columns."""
        return (self.n, self.m)

    shape = property(shape)

    def __str__(self):
        """Pretty print (ugly)."""
        return str(self.data)

    def __repr__(self):
        """Return string capable of copying the object when calling
        C{eval(repr(self))}"""
        return "BlockMatrix(%d, %d)" % (self.n, self.m)

    def __setitem__(self, idx, entry):
        """Set a block matrix entry."""
        if (self._check_tuple(idx)):
            i = idx[0]; j = idx[1]
            self.data[i][j] = entry

    def __getitem__(self, idx):
        """Get a block matrix entry."""
        if (self._check_tuple(idx)):
            return self.data[idx[0]][idx[1]]

    def __add__(self, other):
        """Add two block matrices."""
        if (not isinstance(other, BlockMatrix)):
            raise TypeError, "Can not add a BlockMatrix and a %s" %(str(other))
        if (not self.shape == other.shape):
            raise ValueError, "The BlockMatrices have different shapes: %s, %s" % (str(self.shape), str(other.shape))
        b = Blockmatrix(self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                b[i,j] = self[i]+other[j]
        return b

    def __mul__(self, other):
        """Matrix vector product."""
        if (not isinstance(other, BlockVector)):
            raise TypeError, "Can not multiply BlockMatrix with %s" % (str(type(other)))
        if (not other.n == self.m):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the BlockMatrix (%d, %d)" % (other.n, self.n, self.m)
        res = BlockVector()
        res.n = self.n
        if self.rdim == None:
            self.check()
        for i in range(self.n):
            n = self.rdim[i]
# OLD CODE:
#            res.data.append(numpy.zeros(n, dtype='d'))
# NEW CODE:
            res.data.append(other.data[i].copy())
            for j in range(self.m):
# OLD CODE: 
#                tmp = self.data[i][j]*other.data[j]
#                res[i] += self.data[i][j]*other.data[j]
# NEW CODE: 
                tmp = res[i].copy()
                self.data[i][j].mult(other.data[j], tmp)
                res[i].add(tmp)
        return res

    def prod(self, x, b,transposed=False): #Compute b = A*x
        """Bug in prod!!"""
        if (not isinstance(x, BlockVector)):
            raise TypeError, "Can not multiply BlockMatrix with %s" % (str(type(other)))
        if (not x.n == self.m):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the BlockMatrix (%d, %d)" % (x.n, self.n, self.m)
        if not transposed:
            for i in range(self.n):
                b.data[i] = numpy.multiply(b.data[i], 0.0)
                for j in range(self.m):
                    self.data[i][j].mxv(x.data[j], b.data[i])
        else:
            for i in xrange(self.n):
                b.data[i] = numpy.multiply(b.data[i], 0.0)
                for j in xrange(self.m):
                    self.data[j][i].mxv(x.data[j], b.data[i])

    def _check_tuple(self, idx):
        """Assure that idx is a tuple containing two-elements inside the
        range (self.n-1, self.m-1)"""
        if (not isinstance(idx, tuple)):
            raise TypeError, "wrong type: %s as index" % (str(type(idx)))
        if (not len(idx) == 2):
            raise ValueError, "length of index tuple must be 2, got %d" %(len(idx))
        if (idx[0] > self.n-1):
            raise ValueError, "index out of range, %d" %(idx[0])
        if (idx[1] > self.m-1):
            raise ValueError, "index out of range, %d" %(idx[1])
        return True

    def collect(self):
        """Convert a block matrix to crs"""
        from MatSparse import MapMatSparse, CRS
        B = MapMatSparse()
        for i in xrange(self.shape[0]):
            if i==0: row_offset=0
            else: row_offset = self.rdim[i-1]
            for j in xrange(self.shape[1]):
                if j==0: col_offset=0
                else: col_offset = self.cdim[j-1]
                B.add(self[i,j], False, row_offset, col_offset)
        return CRS(B)

class LowerTriagPrecBlockMatrix(object):
    """A simple block lower triangular matrix implementation in Python.
    The blocks are typically
    matrices implemented in compile languages and later exposed to Python
    either manually or using some sort of wrapper generator software."""

    def __init__(self, *args):
        """Input arguments are either two integer values
        specifying the number of block rows and columns, or a nested Python
        list containing the block elements directly."""
        if isinstance(args[0], int):
            n,m = args
            self.n = n
            self.m = m
            self.data = []
            for i in range(n):
                self.data.append([])
                for j in range(m):
                    self.data[i].append(0)
        elif isinstance(args[0], (list, tuple)):
            self.n = len(args)
            self.m = len(args[0])
            self.data = list(args)


    def shape(self):
        """Return the number of block rows and columns."""
        return (self.n,self.m)

    shape = property(shape)

    def __str__(self):
        """Pretty print (ugly)."""
        return str(self.data)

    def __repr__(self):
        """Return string capable of copying the object when calling
        C{eval(repr(self))}"""
        return "BlockMatrix(%d, %d)" % (self.n, self.m)

    def __setitem__(self, idx, entry):
        """Set a block matrix entry."""
        if (self._check_tuple(idx)):
            i = idx[0]; j = idx[1]
            self.data[i][j] = entry

    def __getitem__(self, idx):
        """Get a block matrix entry."""
        if (self._check_tuple(idx)):
            return self.data[idx[0]][idx[1]]

    def __add__(self, other):
        """Add two block matrices."""
        if (not isinstance(other, BlockMatrix)):
            raise TypeError, "Can not add a BlockMatrix and a %s" %(str(other))
        if (not self.shape == other.shape):
            raise ValueError, "The BlockMatrices have different shapes: %s, %s" % (str(self.shape), str(other.shape))
        b = Blockmatrix(self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                b[i,j] = self[i]+other[j]
        return b

    def __mul__(self, other):
        """Matrix vector product."""
        if (not isinstance(other, BlockVector)):
            raise TypeError, "Can not multiply BlockMatrix with %s" % (str(type(other)))
        if (not other.n == self.m):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the BlockMatrix (%d, %d)" % (other.n, self.n, self.m)
        res = BlockVector()
        res.n = self.n
        for i in range(self.n):
            n = self.data[i][0].shape[0]
            res.data.append(numpy.zeros(n, dtype='d'))
            tmp = other.data[i].copy()#numpy.zeros(n, dtype='d')
            for j in range(i):
                tmp -= self.data[i][j]*res[j]
            res[i]=self.data[i][i]*tmp
        return res

    def prod(self, x, b,transposed=False): #Compute b = A*x
        """Bug in prod!!"""
        print "LowerTriagPrecBlockMatrix::prod"
        if (not isinstance(x, BlockVector)):
            raise TypeError, "Can not multiply BlockMatrix with %s" % (str(type(other)))
        if (not x.n == self.m):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the BlockMatrix (%d, %d)" % (x.n, self.n, self.m)
        if not transposed:
            for i in range(self.n):
                b.data[i] = numpy.multiply(b.data[i], 0.0)
                for j in range(i):
                    self.data[i][j].mxv(x.data[j], b.data[i])
        else:
            for i in range(self.n):
                b.data[i] = numpy.multiply(b.data[i], 0.0)
                for j in range(self.m):
                    self.data[j][i].transposed = True
                    self.data[j][i].mxv(x.data[j], b.data[i])
                    self.data[j][i].transposed = False
            
#                res[i] += self.data[i][j]*other.data[j]
#        return res

    def _check_tuple(self, idx):
        """Assure that idx is a tuple containing two-elements inside the
        range (self.n-1, self.m-1)"""
        if (not isinstance(idx, tuple)):
            raise TypeError, "wrong type: %s as index" % (str(type(idx)))
        if (not len(idx) == 2):
            raise ValueError, "length of index tuple must be 2, got %d" %(len(idx))
        if (idx[0] > self.n-1):
            raise ValueError, "index out of range, %d" %(idx[0])
        if (idx[1] > self.m-1):
            raise ValueError, "index out of range, %d" %(idx[1])
        return True


class SSORPrecBlockMatrix(object):
    """A simple Symmetric Gauss Seidel block  matrix implementation in Python.
    The blocks are typically
    matrices implemented in compile languages and later exposed to Python
    either manually or using some sort of wrapper generator software."""

    def __init__(self, *args):
        """Input arguments are either two integer values
        specifying the number of block rows and columns, or a nested Python
        list containing the block elements directly."""
        if isinstance(args[0], int):
            n,m = args
            self.n = n
            self.m = m
            self.data = []
            self.A = None
            self.iter = 1
            self.w    = 1 # SGS
            for i in range(n):
                self.data.append([])
                for j in range(m):
                    self.data[i].append(0)
        elif isinstance(args[0], (list, tuple)):
            self.n = len(args)
            self.m = len(args[0])
            self.data = list(args)

    def setA(self,A_):
        self.A=A_

    def setIter(self,it_):
        self.iter = it_

    def setRelaxation(self,w_):
        self.w = w_

    def shape(self):
        """Return the number of block rows and columns."""
        return (self.n,self.m)

    shape = property(shape)

    def __str__(self):
        """Pretty print (ugly)."""
        return str(self.data)

    def __repr__(self):
        """Return string capable of copying the object when calling
        C{eval(repr(self))}"""
        return "BlockMatrix(%d, %d)" % (self.n, self.m)

    def __setitem__(self, idx, entry):
        """Set a block matrix entry."""
        if (self._check_tuple(idx)):
            i = idx[0]; j = idx[1]
            self.data[i][j] = entry

    def __getitem__(self, idx):
        """Get a block matrix entry."""
        if (self._check_tuple(idx)):
            return self.data[idx[0]][idx[1]]

    def __add__(self, other):
        """Add two block matrices."""
        if (not isinstance(other, BlockMatrix)):
            raise TypeError, "Can not add a BlockMatrix and a %s" %(str(other))
        if (not self.shape == other.shape):
            raise ValueError, "The BlockMatrices have different shapes: %s, %s" % (str(self.shape), str(other.shape))
        b = Blockmatrix(self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                b[i,j] = self[i]+other[j]
        return b

    def __mul__(self, other):
        """Matrix vector product."""
        if (not isinstance(other, BlockVector)):
            raise TypeError, "Can not multiply BlockMatrix with %s" % (str(type(other)))
        if (not other.n == self.m):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the BlockMatrix (%d, %d)" % (other.n, self.n, self.m)
        res = BlockVector()
        res.n = self.n
        
        # We append memory to the return BlockVector
        for i in xrange(self.n):
            n = self.data[i][0].shape[0]
            res.data.append(numpy.zeros(n, dtype='d'))


        if self.iter==1:
            # Forward SOR
            for i in xrange(self.n):
                tmp = other.data[i].copy()
                for j in xrange(i):
                    tmp -= self.data[i][j]*res[j]
                #for j in xrange(i+1,self.n):
                #    tmp -= self.data[i][j]*res[j]
                res[i]+=self.w*(self.data[i][i]*tmp-res[i])

            # Backward SOR
            for i in xrange(self.n-1,-1,-1):
                tmp = other.data[i].copy()
                for j in xrange(i):
                    tmp -= self.data[i][j]*res[j]
                for j in xrange(i+1,self.n):
                    tmp -= self.data[i][j]*res[j]
                res[i]+=self.w*(self.data[i][i]*tmp-res[i])

        else:
            for k in xrange(self.iter):
                # Forward SOR
                for i in xrange(self.n):
                    tmp = other.data[i].copy()
                    for j in xrange(i):
                        tmp -= self.data[i][j]*res[j]
                    for j in xrange(i+1,self.n):
                        tmp -= self.data[i][j]*res[j]
                    res[i]+=self.w*(self.data[i][i]*tmp-res[i])

                # Backward SOR
                for i in xrange(self.n-1,-1,-1):
                    tmp = other.data[i].copy()
                    for j in xrange(i):
                        tmp -= self.data[i][j]*res[j]
                    for j in xrange(i+1,self.n):
                        tmp -= self.data[i][j]*res[j]
                    res[i]+=self.w*(self.data[i][i]*tmp-res[i])
        return res


    def prod(self, x, b,transposed=False): #Compute b = A*x
        """Bug in prod!!"""
        print "SSORPrecBlockMatrix::prod"
        if (not isinstance(x, BlockVector)):
            raise TypeError, "Can not multiply BlockMatrix with %s" % (str(type(other)))
        if (not x.n == self.m):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the BlockMatrix (%d, %d)" % (x.n, self.n, self.m)
        if not transposed:
            for i in range(self.n):
                b.data[i] = numpy.multiply(b.data[i], 0.0)
                for j in range(i):
                    self.data[i][j].mxv(x.data[j], b.data[i])
        else:
            for i in range(self.n):
                b.data[i] = numpy.multiply(b.data[i], 0.0)
                for j in range(self.m):
                    self.data[j][i].transposed = True
                    self.data[j][i].mxv(x.data[j], b.data[i])
                    self.data[j][i].transposed = False
            
#                res[i] += self.data[i][j]*other.data[j]
#        return res

    def _check_tuple(self, idx):
        """Assure that idx is a tuple containing two-elements inside the
        range (self.n-1, self.m-1)"""
        if (not isinstance(idx, tuple)):
            raise TypeError, "wrong type: %s as index" % (str(type(idx)))
        if (not len(idx) == 2):
            raise ValueError, "length of index tuple must be 2, got %d" %(len(idx))
        if (idx[0] > self.n-1):
            raise ValueError, "index out of range, %d" %(idx[0])
        if (idx[1] > self.m-1):
            raise ValueError, "index out of range, %d" %(idx[1])
        return True





class DiagBlockMatrix(object):
    """A diagonal block matrix implementation in Python. This class is
    typically used to store block diagonal preconditioners"""

    def __init__(self, input):
        """Constructor. 
        
        Input arguments are either an integer value specifying the number of
        block rows, or a Python list containing the diagonal block elements
        directly."""
        if isinstance(input,  int):
            n = input
            self.n = n
            self.data = []
            for i in range(n):
                self.data.append([])
        elif isinstance(input, (list, tuple)):
            self.n = len(input)
            self.data = input[:]


    def shape(self):
        """Return the number of block rows."""
        return (self.n,)

    shape = property(shape)

    def __str__(self):
        """Pretty print (ugly)."""
        return str(self.data)

    def __repr__(self):
        """Return string capable of copying the object when calling
        eval(repr(self))"""
        return "BlockMatrix(%d)" % (self.n)

    def __setitem__(self, idx, entry):
        """Set a diagonal block matrix entry."""
        if (self._check_tuple(idx)):
            if not idx[0] == idx[1] :
                raise ValueError, "diagonal block matrix can not assign to indices (%d, %d)" % (idx[0], idx[1])
            self.data[idx[0]] = entry

    def __getitem__(self, idx):
        """Get a diagonal block matrix entry."""
        if (self._check_tuple(idx)):
            if not idx[0] == idx[1] :
                raise ValueError, "diagonal block matrix dows not provide indices (%d, %d)" % (idx[0], idx[1])
            return self.data[idx[0]]

    def __add__(self, other):
        """Add two diagonal block matrices."""
        if (not isinstance(other, DiagBlockMatrix)):
            raise TypeError, "Can not add a DiagBlockMatrix and a %s" %(str(other))
        if (not self.shape == other.shape):
            raise ValueError, "The BlockMatrices have different shapes: %s, %s" % (str(self.shape), str(other.shape))
        b = DiagBlockmatrix(self.n)
        for i in range(self.n):
            b.data[i] = self[i]+other[i]
        return b

    def __mul__(self, other):
        """Matrix vector product."""
        if (not isinstance(other, BlockVector)):
            raise TypeError, "Can not multiply DiagBlockMatrix with %s" % (str(type(other)))
        if (not other.n == self.n):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the DiagBlockMatrix (%d)" % (other.n, self.n)
        res = BlockVector()
        res.n = self.n
        for i in range(self.n):
            res.data.append(self.data[i]*other.data[i])
        return res


    def prod(self, x, b,transposed=False): #Compute b = A*x
        """Bug in prod!!"""
        if (not isinstance(x, BlockVector)):
            raise TypeError, "Can not multiply DiagBlockMatrix with %s" % (str(type(other)))
        if (not x.n == self.n):
            raise ValueError, "The length of the BlockVector (%d) does not match the dimention of the DiagBlockMatrix (%d, %d)" % (x.n, self.n, self.n)
        if not transposed:
            for i in range(self.n):
                #b.data[i] = numpy.multiply(b.data[i], 0.0)
                #self.data[i].mxv(x.data[i], b.data[i])
                b.data[i]=self.data[i]*x.data[i]
        else:
            for i in range(self.n):
                #b.data[i] = numpy.multiply(b.data[i], 0.0)
                self.data[i].transposed = True
                b.data[i]=self.data[i]*x.data[i]
                self.data[i].transposed = False
            
#                res[i] += self.data[i][j]*other.data[j]
#        return res



    def _check_tuple(self, idx):
        """Assure that idx is a tuple containing two-elements inside the
        range (self.n-1, self.m-1)"""
        if (not isinstance(idx, tuple)):
            raise TypeError, "wrong type: %s as index" % (str(type(idx)))
        if (not len(idx) == 2):
            raise ValueError, "length of index tuple must be 2, got %d" %(len(idx))
        if (idx[0] > self.n-1):
            raise ValueError, "index out of range, %d" %(idx[0])
        if (idx[1] > self.n-1):
            raise ValueError, "index out of range, %d" %(idx[1])
        return True



class BlockVector(object):
    """A block vector implementation in Python. 
    
    The blocks in \c BlockVector are numpy arrays. Instances of this class
    can be matrix vector multiplied with both BlockVectors and
    DiagBlockVectors."""

    def __init__(self, *x):
        if (x):
            self.data = list(x)
            self.n = len(x)
        else:
            self.data = []
            self.n = 0

    def shape(self):
        """Return the number of blocks."""
        return self.n
    
    shape = property(shape)

    def norm(self):
        return reduce(operator.add, map(lambda x: numpy.dot(x,x), self.data))

    def inner(self, other):
#   OLD CODE: 
#        return reduce(operator.add, map(lambda x,y: numpy.dot(x,y), self.data, other.data))
#   NEW CODE: 
        return reduce(operator.add, map(lambda x,y: x.inner(y), self.data, other.data))

        return reduce(operator.add, map(lambda x,y: numpy.dot(x,y), self.data, other.data))

    def __add__(self, other):
        res = BlockVector()
        res.n = self.n
        for i in range(self.n):
#   OLD CODE: 
#            res.data.append(self.data[i]+other.data[i])
#   NEW CODE: 
            tmp = self.data[i].copy()
            tmp.add(other.data[i], 1.0)
            res.data.append(tmp)
        return res

    def __sub__(self, other):
        res = BlockVector()
        res.n = self.n
        for i in range(self.n):
# OLD CODE: 
#            res.data.append(self.data[i]-other.data[i])
# NEW CODE: 
            tmp = self.data[i].copy()
            tmp.add(other.data[i], -1.0)
            res.data.append(tmp)
        return res


    def __mul__(self, other):
        res = BlockVector()
        res.n = self.n
        if isinstance(other,BlockVector):
            for i in range(self.n):
                res.data.append(self.data[i]*other[i])
        else:
            for i in range(self.n):
# OLD CODE: 
#                res.data.append(self.data[i]*other)
# NEW CODE
                tmp = self.data[i].copy()
                tmp.mult(other)
                res.data.append(tmp)
        return res

    def __rmul__(self, other):
        return self.__mul__(other)


    def __imul__(self, other):
        for i in range(self.n):
            self.data[i] *= other
        return self

    def __iadd__(self, other):
        for i in range(self.n):
# OLD CODE:
#            self.data[i] += other.data[i]
# NEW CODE: 
            self.data[i].add(other.data[i], 1.0)
        return self

    def __isub__(self, other):
        for i in range(self.n):
            self.data[i] -= other.data[i]
        return self



    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, x):
        if self.n > i:
            self.data[i] = x
        elif self.n == i:
            self.data.append(x)
            self.n += 1
        else:
            raise "BlockVector error","Unable to attach vector in block %d" % (i)

    def copy(self):
        res = BlockVector()
        res.n = self.n
        for i in range(self.n):
            res.data.append(self.data[i].copy())
        return res

if __name__ == '__main__':
    """Test suite"""

    A = BlockMatrix(2,2)
    u = numpy.arange(10,dtype='d')
    v = numpy.arange(10,dtype='d')
    x = BlockVector(u,v)
    print x.norm()
    print (x+x).norm()
    print (2*x).norm()
    x *=2
    print type(x)
    print x.norm()
    print (x+2*x).norm()

    B = DiagBlockMatrix(2)
