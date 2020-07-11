with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_VecVecs;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;

package Standard_Coefficient_Circuits is

-- DESCRIPTION :
--   Separating real parts from imaginary parts, inlining the algorithms
--   for complex multiplication, fusing loops, results in faster code.

-- DATA STRUCTURES :
--   A circuit stores the exponents and coefficients and hold work space to
--   apply the reverse mode of algorithmic differentiation and evaluation.

  type Circuit ( nbr : integer32 ) is record
    dim : integer32; -- the number of variables in the circuit
    pdg : integer32; -- the polynomial degree of the circuit
    xps : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponents
    idx : Standard_Integer_VecVecs.VecVec(1..nbr); -- indices of exponents
    fac : Standard_Integer_VecVecs.VecVec(1..nbr); -- factor indices
   -- the coefficients are stored as vectors of real and imaginary parts
    rcf : Standard_Floating_Vectors.Vector(1..nbr); -- real parts
    icf : Standard_Floating_Vectors.Vector(1..nbr); -- imaginary parts
   -- the constant is stored as its real and imaginary parts
    rcst : double_float; -- real part of the constant 
    icst : double_float; -- imaginary part of the constant 
   -- the work space for the forward products
    rfwd : Standard_Floating_Vectors.Link_to_Vector; -- real parts
    ifwd : Standard_Floating_Vectors.Link_to_Vector; -- imaginary parts
   -- the work space for the backward products
    rbck : Standard_Floating_Vectors.Link_to_Vector; -- real parts
    ibck : Standard_Floating_Vectors.Link_to_Vector; -- imaginary parts
   -- the work space for the cross products
    rcrs : Standard_Floating_Vectors.Link_to_Vector; -- real parts
    icrs : Standard_Floating_Vectors.Link_to_Vector; -- imaginary parts
  end record;

  type Link_to_Circuit is access Circuit;

  type Circuits is array ( integer32 range <> ) of Link_to_Circuit;

-- A system stores the sequence of circuits for each polynomial,
-- along with work space and the final outcomes.

  type System ( neq,dim : integer32 ) is record
    crc : Circuits(1..neq);                          -- polynomials
    mxe : Standard_Integer_Vectors.Vector(1..dim);   -- exponent maxima
    rpwt : Standard_Floating_VecVecs.VecVec(1..dim); -- power table real part
    ipwt : Standard_Floating_VecVecs.VecVec(1..dim); -- power table imag part
    ryd : Standard_Floating_Vectors.Link_to_Vector;  -- real part of gradient
    iyd : Standard_Floating_Vectors.Link_to_Vector;  -- imag part of gradient
    fx : Standard_Complex_Vectors.Vector(1..neq);    -- function value
    jm : Standard_Complex_Matrices.Matrix(1..neq,1..dim); -- Jacobian matrix
    jrc,jic : Standard_Floating_VecVecs.Link_to_VecVec;
   -- real and imaginary parts of the columns of the Jacobian matrix
    hrp,hip : Standard_Floating_VecVecs.Link_to_VecVec;
   -- work space for the rows of Hessians, real and imaginary parts
  end record;

  type Link_to_System is access System;

  type System_Array is array ( integer32 range<> ) of Link_to_System;

-- CONSTRUCTORS :

  function Allocate ( nbr,dim : integer32 ) return Circuit;

  -- DESCRIPTION :
  --   Returns a circuit for a polynomial with nbr monomials,
  --   with dim variables, and with allocated work space vectors.

  function Exponent_Maxima
             ( c : Circuits; dim : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the maximal exponents of the dim variables in the circuits,
  --   to allocate the power table in a system of circuits.

  function Create ( c : Circuits; dim : integer32 ) return System;

  -- DESCRIPTION :
  --   Given well defined circuits for dimension dim,
  --   computes mxe and allocates space for a system,
  --   including the work space columns for the Jacobian matrix,
  --   as they are defined also by EvalDiff.
  --   The hrp and hip work spaces are not allocated,
  --   as Hessians are optional, only used in EvalDiff2.

  procedure Allocate_Jacobian_Space
              ( neq,dim : in integer32;
                jrc,jic : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Allocates in jrc and jic space for the real and imaginary parts
  --   of the complex numbers in the columns of the Jacobian matrix.
  --   Every vector jrc(k) and jic(k) will have range 1..neq,
  --   where neq is the number of equations = the number of rows
  --   in the Jacobian matrix.
  --   The number of variables or the number of columns in the
  --   Jacobian matrix is given in the value of dim.

  -- REQUIRED : jrc'range = jic'range = 1..dim.

  procedure Allocate_Jacobian_Space
              ( neq,dim : in integer32;
                jrc,jic : out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns in jrc and jic two vectors or range 1..dim.
  --   Every vector in jrc and jic has range 1..neq
  --   and is initialized to zero.

  procedure Allocate_Jacobian_Space ( s : in out System );
  procedure Allocate_Jacobian_Space ( s : in Link_to_System );

  -- DESCRIPTION :
  --   Allocates the space for all columns of the real and imaginary parts
  --   of the Jacobian matrix, in s.jrc and s.jic.

  procedure Allocate_Hessian_Space
              ( dim : in integer32;
                hrp,hip : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Allocates in hrp and hip space for vectors of range 1..dim
  --   and initializes all entries to zero.

  -- REQUIRED : hrp'range = hip'range = 1..dim.

  procedure Allocate_Hessian_Space
              ( dim : in integer32;
                hrp,hip : out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns in hrp and hip two vectors of vectors of range 1..dim.
  --   Each vector in hrp and hip has range 1..dim and is initialized to zero.

  procedure Allocate_Hessian_Space ( s : in out System );
  procedure Allocate_Hessian_Space ( s : in Link_to_System );

  -- DESCRIPTION :
  --   Allocates the space for all rows of the real and imaginary parts
  --   of the Hessian matrices, in s.hrp and s.hip.

  function Merge ( hrp,hip : Standard_Floating_VecVecs.VecVec )
                 return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Merges the real and imaginary parts of the rows in hrp and hip
  --   into one complex matrix.

  -- REQUIRED : hrp'range = hip'range,
  --   and for all k in hrp'range: hrp(k)'range = hip(k)'range.

  function Copy ( c : Circuit ) return Circuit;
  function Copy ( c : Link_to_Circuit ) return Link_to_Circuit;
  function Copy ( c : Circuits ) return Circuits;

  -- DESCRIPTION :
  --   Returns a deep copy of the circuit (or circuits) c.
  --   The work space for the forward, backward, and cross products
  --   is allocated but the values are not copied from c.

  function Copy ( s : System ) return System;
  function Copy ( s : Link_to_System ) return Link_to_System;

  -- DESCRIPTION :
  --   Returns a deep copy of the circuits in s
  --   and allocates all work space data structures,
  --   including the data structures for the Hessians.

-- RADIUS COEFFICIENTS :

  procedure AbsVal ( c : in out Circuit );
  procedure AbsVal ( c : in Link_to_Circuit );
  procedure AbsVal ( c : in Circuits );
  procedure AbsVal ( s : in System );
  procedure AbsVal ( s : in Link_to_System );

  -- DESCRIPTION :
  --   Replaces all coefficients in c by their radius.

-- EVALUATION OF CIRCUITS :

  function Eval ( c : in Circuit;
                  xr : in Standard_Floating_Vectors.Link_to_Vector;
                  rpwt : in Standard_Floating_VecVecs.VecVec )
                return double_float;

  -- DESCRIPTION :
  --   Returns the value of the circuit at a vector with zero
  --   imaginary parts, with the aid of a power table.

  function Eval ( c : in Circuit;
                  xr : in Standard_Floating_Vectors.Link_to_Vector;
                  xi : in Standard_Floating_Vectors.Link_to_Vector;
                  rpwt : in Standard_Floating_VecVecs.VecVec;
                  ipwt : in Standard_Floating_VecVecs.VecVec ) 
                return Standard_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the circuit c at the vector with
  --   real parts in xr and imaginary parts in xi,
  --   with the aid of the power table, with real parts in rpwt
  --   and imaginary parts in ipwt.

  procedure Eval ( s : in out System;
                   xr : in Standard_Floating_Vectors.Link_to_Vector;
                   xi : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Eval ( s : in Link_to_System;
                   xr : in Standard_Floating_Vectors.Link_to_Vector;
                   xi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Updates the power table with the values of xr and xi
  --   and returns the value of the system in s.fx.

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUITS :

  procedure EvalDiff
              ( s : in out System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector );
  procedure EvalDiff
              ( s : in Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in s at x.

  -- REQUIRED :
  --   All space for the power table and yd has been allocated.

  -- ON ENTRY :
  --   s        properly defined and allocated system of circuits;
  --   xr       real parts of values for the variables in the system;
  --   xi       imaginary parts of values for the variables in the system.

  -- ON RETURN :
  --   s.pwt    power table updated for the values in x;
  --   s.fx     function value of the circuits at x;
  --   s.jm     the Jacobian matrix evaluated at x.

  procedure EvalDiff2
              ( s : in out System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat );
  procedure EvalDiff2
              ( s : in Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in s at x,
  --   and computes the values of all Hessian matrices at x.

  -- REQUIRED :
  --   All space for the power table and yd has been allocated.
  --   Also s.hrp and s.hip are allocated.

  -- ON ENTRY :
  --   s        properly defined and allocated system of circuits;
  --   xr       real parts of values for the variables in the system;
  --   xi       imaginary parts of values for the variables in the system.
  --   vh       space allocated for s.neq matrices,
  --            all matrices have 1..dim for range(1) and range(2).

  -- ON RETURN :
  --   s.pwt    power table updated for the values in x;
  --   s.fx     function value of the circuits at x;
  --   s.jm     the Jacobian matrix evaluated at x;
  --   vh       vector of evaluated Hessian matrices.

  procedure EvalDiff
              ( c : in Circuits;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                jrc : in Standard_Floating_VecVecs.VecVec;
                jic : in Standard_Floating_VecVecs.VecVec;
                fx : out Standard_Complex_Vectors.Vector;
                jm : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in c at x.

  -- ON ENTRY :
  --   c        a sequence of circuits, properly defined and allocated;
  --   xr       real parts of the values for the variables;
  --   xi       imaginary parts of the values for the variables;
  --   ryd      work space for the real part of function value and gradient,
  --            of range 0..dim, where dim = x'last;
  --   ryd      work space for the imag part of function value and gradient,
  --            of range 0..dim, where dim = x'last;
  --   rpwt     real part of power table defined and computed for x;
  --   ipwt     imag part of power table defined and computed for x;
  --   jrc      space allocated for the real parts of the complex numbers
  --            in the columns of the Jacobian matrix;
  --   jic      space allocated for the imaginary parts of the complex numbers
  --            in the columns of the Jacobian matrix.

  -- ON RETURN :
  --   jrc      real parts of the columns of the Jacobian matrix;
  --   jic      imaginary parts of the columns of the Jacobian matrix;
  --   fx       vector of function values of the circuits at x;
  --   jm       matrix of partial derivatives.

  procedure EvalDiff2
              ( c : in Circuits;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                jrc : in Standard_Floating_VecVecs.VecVec;
                jic : in Standard_Floating_VecVecs.VecVec;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec;
                fx : out Standard_Complex_Vectors.Vector;
                jm : out Standard_Complex_Matrices.Matrix;
                vh : in Standard_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in c at x,
  --   returns the evaluated Jacobian and vector of Hessians.

  -- ON ENTRY :
  --   c        a sequence of circuits, properly defined and allocated;
  --   xr       real parts of the values for the variables;
  --   xi       imaginary parts of the values for the variables;
  --   ryd      work space for the real part of function value and gradient,
  --            of range 0..dim, where dim = x'last;
  --   ryd      work space for the imag part of function value and gradient,
  --            of range 0..dim, where dim = x'last;
  --   rpwt     real part of power table defined and computed for x;
  --   ipwt     imag part of power table defined and computed for x.
  --   jrc      space allocated for the real parts of the complex numbers
  --            in the columns of the Jacobian matrix;
  --   jic      space allocated for the imaginary parts of the complex numbers
  --            in the columns of the Jacobian matrix.
  --   hrp      works space for the real parts of Hessian matrices;
  --   hip      works space for the imaginary parts of Hessian matrices;
  --   vh       space allocated for dim matrices,
  --            all matrices have 1..dim for range(1) and range(2).

  -- ON RETURN :
  --   jrc      real parts of the columns of the Jacobian matrix;
  --   jic      imaginary parts of the columns of the Jacobian matrix;
  --   fx       vector of function values of the circuits at x;
  --   jm       matrix of partial derivatives;
  --   vh       vector of evaluated Hessian matrices.

-- SINGULAR VALUE DECOMPOSITIONS :

  procedure Singular_Values
              ( s : in out System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                svls : in Standard_Complex_VecVecs.VecVec );
  procedure Singular_Values
              ( s : in Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                svls : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in s at x,
  --   computes the values of all Hessian matrices at x,
  --   computes all singular values of the Jacobian matrix
  --   and of all Hessian matrices.

  -- REQUIRED :
  --   All space for the power table and yd has been allocated.

  -- ON ENTRY :
  --   s        properly defined and allocated system of circuits;
  --   xr       real parts of the values for the variables in the system;
  --   xi       imaginary parts of the values for the variables in the system;
  --   vh       space allocated for s.neq matrices,
  --            all matrices have 1..dim for range(1) and range(2);
  --   svls     vector of range 0..s.neq, with allocated space
  --            for all vectors of singular values,
  --            the range of s(k) is 1..s.dim+1.

  -- ON RETURN :
  --   s.pwt    power table updated for the values in x;
  --   s.fx     function value of the circuits at x;
  --   s.jm     the Jacobian matrix evaluated at x,
  --            but destroyed by the SVD computation;
  --   vh       vector of evaluated Hessian matrices,
  --            by destroyed by the SVD computation;
  --   svls     svls(0) contains the singular values of s.jm, and
  --            svls(k) contains the singular values of vh(k),
  --            for k in vh'range.

  procedure Singular_Values
              ( c : in Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Vectors.Vector );
  procedure Singular_Values
              ( c : in Link_to_Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates the circuit c at x, computes the Hessian matrix A,
  --   and computes the singular value decomposition of A.
  --   In the special case of c.pdg = 1, all values in s are set to zero
  --   and nothing else is computed.

  -- REQUIRED :
  --   x'range = 1..c.dim, yd'range = 0..c.dim,
  --   pwt'range = x'range, pwt(k)'range extends to the maximal exponent
  --   of the k-th variable in the circuit,
  --   A'range(1) = A'range(2) = x'range,
  --   U'range(1) = U'range(2) = x'range, V'range(1) = V'range(2) = x'range,
  --   e'range = 1..c.dim, s'range = 1..c.dim+1.

  -- ON ENTRY :
  --   c        a circuit, properly defined and allocated;
  --   xr       real parts of the values for the variables;
  --   xi       imaginary parts of the values for the variables;
  --   ryd      work space for the real parts of the function value
  --            and gradient, of range 0..dim, where dim = xr'last;
  --   iyd      work space for the imaginary parts of the function value
  --            and gradient, of range 0..dim, where dim = xi'last;
  --   rpwt     real parts of the power table, for xr and xi;
  --   ipwt     imaginary parts of the power table, for xr and xi;
  --   hrp      works space for the real parts of Hessian matrices;
  --   hip      works space for the imaginary parts of Hessian matrices.

  -- ON RETURN :
  --   A        the Hessian matrix of c at x, 
  --            however, its values are destroyed by the SVD;
  --   U        the U matrix in the SVD of A;
  --   V        the V matrix in the SVD of A;
  --   e        contains error information on the SVD computation;
  --   s        the first c.dim entries are the singular values of A.

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUIT :
--   The Indexed_Speel procedures are for circuits where the exponents
--   are either zero or one.  There are an important subclass to deal
--   with monomials that have no higher powers.
--   The general case is handled by the Speel procedures,
--   with wrappers working on circuits.
--   Both indexed and general speel procedure compute the gradient
--   and optionally, the evaluated Hessian matrix.

  procedure Indexed_Speel 
              ( c : in Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates and differentiates the circuit c at x
  --   and stores the result in yd.
  --   Wraps the next Speel procedure, using the c.xps as indices.
  --   The results are correct if all monomials are products of variables,
  --   with no exponent higher than one.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   hrp      vector of range 0..c.dim, with space allocated for
  --            the rows of the real parts of the Hessian;
  --   hip      vector of range 0..c.dim, with space allocated for
  --            the rows of the imaginary parts of the Hessian.      

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x;
  --   hrp      rows of the real parts of the Hessian at x;
  --   hip      rows of the imaginary parts of the Hessian at x.

  procedure Indexed_Speel
              ( idx : in Standard_Integer_Vectors.Vector;
                rcff,icff : in double_float;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                rbck : in Standard_Floating_Vectors.Link_to_Vector;
                ibck : in Standard_Floating_Vectors.Link_to_Vector;
                rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                icrs : in Standard_Floating_Vectors.Link_to_Vector;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates an indexed product, multiplied with a coefficient,
  --   computes its gradient and updates the Hessian matrix.
  --   Is called frequently by the next Indexed_Speel procedure.

  -- REQUIRED : idx'last >= 2.
  --   xr'range = xi'range = 1..dim and ryd'range = iyd'range = 0..dim,
  --   rfwd'range = ifwd'range = 1..dim-1,
  --   rbck'range = ibck'range = 1..dim-2 = rcrs'range = icrs'range,
  --   hrp'range = hip'range = 1..dim, and all vectors in hrp and hip
  --   have range 1..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in the monomial;
  --   rcff     real part of the coefficient of the monomial;
  --   icff     imaginary part of the coefficient of the monomial;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   rfwd     work space vector of range 1..dim-1,
  --            for the real parts of the forward products;
  --   ifwd     work space vector of range 1..dim-1,
  --            for the imaginary parts of the forward products;
  --   rbck     work space vector of range 1..dim-2,
  --            for the real parts of the backward products;
  --   ibck     work space vector of range 1..dim-2,
  --            for the imaginary parts of the backward products;
  --   rcrs     work space vector of range 1..dim-2,
  --            for the real parts of the cross products.
  --   icrs     work space vector of range 1..dim-2,
  --            for the imaginary parts of the cross products;
  --   hrp      vector of range 0..c.dim, with space allocated for
  --            the rows of the real parts of the Hessian;
  --   hip      vector of range 0..c.dim, with space allocated for
  --            the rows of the imaginary parts of the Hessian.      

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x;
  --   hrp      rows of the real parts of the Hessian at x,
  --            with updated upper triangular part;
  --   hip      rows of the imaginary parts of the Hessian at x,
  --            with updated upper triangular part.

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                rcf : in Standard_Floating_Vectors.Vector;
                icf : in Standard_Floating_Vectors.Vector;
                rcst,icst : in double_float;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                rbck : in Standard_Floating_Vectors.Link_to_Vector;
                ibck : in Standard_Floating_Vectors.Link_to_Vector;
                rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                icrs : in Standard_Floating_Vectors.Link_to_Vector;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Implements the Speel procedure on the circuit c.

  -- REQUIRED :
  --   idx'range = rcf'range = icf'range and all vectors in idx have 
  --   values in range 1..dim, where dim is the number of variables,
  --   xr'range = xi'range = 1..dim and ryd'range = iyd'range = 0..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in each monomial;
  --   rcf      real parts of the coefficients of the monomials;
  --   icf      imaginary parts of the coefficients of the monomials;
  --   rcst     real part of the constant coefficient of the circuit;
  --   icst     imaginary part of the constant coefficient of the circuit;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   rfwd     work space vector of range 1..dim-1,
  --            for the real parts of the forward products;
  --   ifwd     work space vector of range 1..dim-1,
  --            for the imaginary parts of the forward products;
  --   rbck     work space vector of range 1..dim-2,
  --            for the real parts of the backward products;
  --   ibck     work space vector of range 1..dim-2,
  --            for the imaginary parts of the backward products;
  --   rcrs     work space vector of range 1..dim-2,
  --            for the real parts of the cross products.
  --   icrs     work space vector of range 1..dim-2,
  --            for the imaginary parts of the cross products;
  --   hrp      vector of range 0..c.dim, with space allocated for
  --            the rows of the real parts of the Hessian;
  --   hip      vector of range 0..c.dim, with space allocated for
  --            the rows of the imaginary parts of the Hessian.      

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x;
  --   hrp      rows of the real parts of the Hessian at x;
  --   hip      rows of the imaginary parts of the Hessian at x.

  procedure Indexed_Speel 
              ( c : in Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluates and differentiates the circuit c at x
  --   and stores the result in yd.
  --   Wraps the next Speel procedure, using the c.xps as indices.
  --   The results are correct if all monomials are products of variables,
  --   with no exponent higher than one.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result.

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x.

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                rcf : in Standard_Floating_Vectors.Vector;
                icf : in Standard_Floating_Vectors.Vector;
                rcst,icst : in double_float;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                rbck : in Standard_Floating_Vectors.Link_to_Vector;
                ibck : in Standard_Floating_Vectors.Link_to_Vector;
                rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                icrs : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Implements the Speel procedure on the circuit c.

  -- REQUIRED :
  --   idx'range = rcf'range = icf'range and all vectors in idx have 
  --   values in range 1..dim, where dim is the number of variables,
  --   xr'range = xi'range = 1..dim and ryd'range = iyd'range = 0..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in each monomial;
  --   rcf      real parts of the coefficients of the monomials;
  --   icf      imaginary parts of the coefficients of the monomials;
  --   rcst     real part of the constant coefficient of the circuit;
  --   icst     imaginary part of the constant coefficient of the circuit;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   rfwd     work space vector of range 1..dim-1,
  --            for the real parts of the forward products;
  --   ifwd     work space vector of range 1..dim-1,
  --            for the imaginary parts of the forward products;
  --   rbck     work space vector of range 1..dim-2,
  --            for the real parts of the backward products;
  --   ibck     work space vector of range 1..dim-2,
  --            for the imaginary parts of the backward products;
  --   rcrs     work space vector of range 1..dim-2,
  --            for the real parts of the cross products.
  --   icrs     work space vector of range 1..dim-2,
  --            for the imaginary parts of the cross products.

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x.

  procedure Speel ( c : in Circuit;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   indexed sequence of products.
  --   Wraps the next Speel procedure, using the data in c.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   rpwt     real parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1);
  --   rpwt     imaginary parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1).

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x.

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcf : in Standard_Floating_Vectors.Vector;
                    icf : in Standard_Floating_Vectors.Vector;
                    rcst,icst : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   indexed sequence of products, with higher powers factored out.

  -- REQUIRED :
  --   idx'range = rcf'range = icf'range and all vectors in idx have 
  --   values in range 1..dim, where dim is the number of variables,
  --   xr'range = xi'range = 1..dim and ryd'range = iyd'range = 0..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in each monomial;
  --   rcf      real parts of the coefficients of the monomials;
  --   icf      imaginary parts of the coefficients of the monomials;
  --   rcst     real part of the constant coefficient of the circuit;
  --   icst     imaginary part of the constant coefficient of the circuit;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   rfwd     work space vector of range 1..dim-1,
  --            for the real parts of the forward products;
  --   ifwd     work space vector of range 1..dim-1,
  --            for the imaginary parts of the forward products;
  --   rbck     work space vector of range 1..dim-2,
  --            for the real parts of the backward products;
  --   ibck     work space vector of range 1..dim-2,
  --            for the imaginary parts of the backward products;
  --   rcrs     work space vector of range 1..dim-2,
  --            for the real parts of the cross products.
  --   icrs     work space vector of range 1..dim-2,
  --            for the imaginary parts of the cross products;
  --   rpwt     real parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1);
  --   rpwt     imaginary parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1).

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x.

  procedure Speel ( c : in Circuit;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec; 
                    hrp : in Standard_Floating_VecVecs.VecVec;
                    hip : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   indexed sequence of products.
  --   Wraps the next Speel procedure, using the data in c.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   rpwt     real parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1);
  --   rpwt     imaginary parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1);
  --   hrp      vector of range 0..c.dim, with space allocated for
  --            the rows of the real parts of the Hessian;
  --   hip      vector of range 0..c.dim, with space allocated for
  --            the rows of the imaginary parts of the Hessian.      

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x;
  --   hrp      real part of the Hessian matrix evaluated at x;
  --   hip      imaginary part of the Hessian matrix evaluated at x.

  procedure Speel ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                    idx : in Standard_Integer_Vectors.Vector;
                    rcff,icff : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Updates the function value and gradient in yd,
  --   for a general term in the circuit.
  --   This is a helper procedure in the next Speel procedure.

  -- REQUIRED : fac /= null.

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcf : in Standard_Floating_Vectors.Vector;
                    icf : in Standard_Floating_Vectors.Vector;
                    rcst,icst : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec;
                    hrp : in Standard_Floating_VecVecs.VecVec;
                    hip : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   indexed sequence of products, with higher powers factored out.

  -- REQUIRED :
  --   idx'range = rcf'range = icf'range and all vectors in idx have 
  --   values in range 1..dim, where dim is the number of variables,
  --   xr'range = xi'range = 1..dim and ryd'range = iyd'range = 0..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in each monomial;
  --   rcf      real parts of the coefficients of the monomials;
  --   icf      imaginary parts of the coefficients of the monomials;
  --   rcst     real part of the constant coefficient of the circuit;
  --   icst     imaginary part of the constant coefficient of the circuit;
  --   xr       vector of range 1..c.dim,
  --            with values for the real parts of x;
  --   xr       vector of range 1..c.dim,
  --            with values for the imaginary parts of x;
  --   ryd      vector of range 0..c.dim,
  --            allocated for the real parts of the result;
  --   iyd      vector of range 0..c.dim,
  --            allocated for the imaginary parts of the result;
  --   rfwd     work space vector of range 1..dim-1,
  --            for the real parts of the forward products;
  --   ifwd     work space vector of range 1..dim-1,
  --            for the imaginary parts of the forward products;
  --   rbck     work space vector of range 1..dim-2,
  --            for the real parts of the backward products;
  --   ibck     work space vector of range 1..dim-2,
  --            for the imaginary parts of the backward products;
  --   rcrs     work space vector of range 1..dim-2,
  --            for the real parts of the cross products.
  --   icrs     work space vector of range 1..dim-2,
  --            for the imaginary parts of the cross products;
  --   rpwt     real parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1);
  --   rpwt     imaginary parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1);
  --   hrp      vector of range 0..c.dim, with space allocated for
  --            the rows of the real parts of the Hessian;
  --   hip      vector of range 0..c.dim, with space allocated for
  --            the rows of the imaginary parts of the Hessian.      

  -- ON RETURN :
  --   ryd(0)   real part of the value of the circuit at x;
  --   iyd(0)   imaginary part of the value of the circuit at x;
  --   ryd(k)   real part of the k-th derivative of the circuit at x;
  --   iyd(k)   imaginary part of the k-th derivative of the circuit at x;
  --   hrp      real part of the Hessian matrix evaluated at x;
  --   hip      imaginary part of the Hessian matrix evaluated at x.

-- AUXILIARY PROCEDURES :

  procedure Forward ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                      xi : in Standard_Floating_Vectors.Link_to_Vector;
                      fr : in Standard_Floating_Vectors.Link_to_Vector;
                      fi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.

  -- REQUIRED :
  --   xr'range = xi'range, fr'first = xr'first = 1,
  --   and fi'last >= xi'last-1.

  procedure Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.

  procedure Fused_Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.
  --   Applies loop fusion.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.

  procedure Forward_Backward_Cross
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.
  --   Computes all cross products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in cr and the imaginary parts in ci.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.
  --   cr'first = cr'first = 1, ci'last >= ci'last-2.

  -- ON ENTRY :
  --   xr       real parts of the values for the variables;
  --   xi       imaginary parts of the values for the variables;
  --   fr       space allocated for real parts of forward products;
  --   fi       space allocated for imaginary parts of forward products;
  --   br       space allocated for real parts of backward products;
  --   br       space allocated for imaginary parts of backward products;
  --   cr       space allocated for real parts of cross products;
  --   cr       space allocated for imaginary parts of cross products.

  -- ON RETURN : let n be x'last-1,
  --   fr(n)    the real part of the product of all variables in x;
  --   fi(n)    the imaginary part of the product of all variables in x;
  --   fr(n-1)  is the real part of the partial derivative of the product
  --            with respect to the last variable;
  --   fi(n-1)  is the imaginary part of the partial derivative
  --            of the product with respect to the last variable;
  --   br(n-2)  is the real part of the partial derivative
  --            of the product with respect to the first variable;
  --   bi(n-2)  is the imaginary part of the partial derivative
  --            of the product with respect to the first variable;
  --   cr(k)    is the real part of the partial derivative 
  --            of the product with respect to the (k+1)-th variable;
  --   ci(k)    is the imaginary part of the partial derivative 
  --            of the product with respect to the (k+1)-th variable.

  procedure Fused_Forward_Backward_Cross
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Fused_Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.
  --   Computes all cross products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in cr and the imaginary parts in ci.
  --   Applies loop fusion.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1, or >= idx'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2, or >= idx'last-2,
  --   cr'first = cr'first = 1, ci'last >= ci'last-2, or >= idx'last-2.

  -- ON ENTRY :
  --   idx      if provided, then only those values of x
  --            as indexed by the entries in idx will be used;
  --   xr       real parts of the values for the variables;
  --   xi       imaginary parts of the values for the variables;
  --   fr       space allocated for real parts of forward products;
  --   fi       space allocated for imaginary parts of forward products;
  --   br       space allocated for real parts of backward products;
  --   br       space allocated for imaginary parts of backward products;
  --   cr       space allocated for real parts of cross products;
  --   cr       space allocated for imaginary parts of cross products.

  -- ON RETURN : let n be x'last-1, or idx'last-1,
  --   fr(n)    the real part of the product of all variables in x,
  --            or only those values indexed by idx;
  --   fi(n)    the imaginary part of the product of all variables in x,
  --            or only those values indexed by idx;
  --   fr(n-1)  is the real part of the partial derivative of the product
  --            with respect to the last variable (or as indexed by idx);
  --   fi(n-1)  is the imaginary part of the partial derivative of the product
  --            with respect to the last variable (or as indexed by idx);
  --   br(n-2)  is the real part of the partial derivative of the product
  --            with respect to the first variable (or idx(1));
  --   bi(n-2)  is the imaginary part of the partial derivative of the product
  --            with respect to the first variable (or idx(1));
  --   cr(k)    is the real part of the partial derivative of the product
  --            with respect to the (k+1)-th variable (or idx(k+1));
  --   ci(k)    is the imaginary part of the partial derivative of the product
  --            with respect to the (k+1)-th variable (or idx(k+1)).

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of range mxe'range with space for
  --   floating-point vectors of range 1..mxe(k)-1, for k in mxe'range.
  --            for i in range 1..mxe(k)-1.

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the power table for the values of the variables,
  --   with real parts in xr and imaginary parts in xi.

  -- REQUIRED :
  --   mxe'range = xr'range = xi'range,
  --   rpwt and ipwt are allocated according to mxe,
  --   rpwt'range = xr'range, rpwt(k)'range = 1..mxe(k)-1,
  --   ipwt'range = xr'range, ipwt(k)'range = 1..mxe(k)-1.

  -- ON ENTRY :
  --   mxe      highest exponents of the variables,
  --            mxe(k) is the highest exponent of the k-th variable
  --   xr       real parts of the values for all variables;
  --   xi       imaginary parts of the values for all variables;
  --   rpwt     allocated memory for the real parts of all powers
  --            of the values of the variables;
  --   ipwt     allocated memory for the imaginary parts of all powers
  --            of the values of the variables.

  -- ON RETURN :
  --   rpwt     real part of the power table,
  --            rpwt(k)(i) equals the real part of x(k)**(i+1),
  --            where x(k) is the complex value of the k-th variable,
  --            for i in range 1..mxe(k)-1;
  --   rpwt     imaginary part of the power table,
  --            rpwt(k)(i) equals the imaginary part of x(k)**(i+1),
  --            where x(k) is the complex value of the k-th variable,
  --            for i in range 1..mxe(k)-1.

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                rcf : in double_float;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                rpf : out double_float );

  -- DESCRIPTION :
  --   Computes the common factor, for higher powers of variables,
  --   in case all imaginary parts are zero.

  -- ON ENTRY :
  --   xps      values of the exponents for the powers of x;
  --   fac      factor indices;
  --   xr       real parts for the variables used for low powers;
  --   rcf      real part of the coefficient of the monomial;
  --   rpwt     real parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1).

  -- ON RETURN :
  --   rpf      the evaluated common factor.

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                rcf,icf : in double_float;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                rpf,ipf : out double_float );

  -- DESCRIPTION :
  --   Computes the common factor, for higher powers of variables.

  -- ON ENTRY :
  --   xps      values of the exponents for the powers of x;
  --   fac      factor indices;
  --   xr       real parts for the variables used for low powers;
  --   xi       imaginary parts for the variables used for low powers;
  --   rcf      real part of the coefficient of the monomial;
  --   icf      imaginary part of the coefficient of the monomial;
  --   rpwt     real parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1);
  --   rpwt     imaginary parts of the power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1).

  -- ON RETURN :
  --   rpf      real part of the common factor;
  --   ipf      imaginary part of the common factor.

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit );
  procedure Clear ( c : in out Link_to_Circuit );
  procedure Clear ( c : in out Circuits );

  -- DESCRIPION :
  --   Deallocates the space occupied by the circuit.

  procedure Clear ( s : in out System );
  procedure Clear ( s : in out Link_to_System );
  procedure Clear ( s : in out System_Array );

  -- DESCRIPTION :
  --   Deallocates the space occupied by s.

end Standard_Coefficient_Circuits;
