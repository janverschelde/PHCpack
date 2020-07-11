with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with DoblDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;

package DoblDobl_Complex_Circuits is

-- DESCRIPTION :
--   This package contains development code for standard_coefficient_circuits,
--   to test the correctness of the better performing algorithms for
--   algorithmic differentiation and evaluation.
--   Computations happen in double double precision.

-- DATA STRUCTURES :
--   A circuit stores the exponents and coefficients and hold work space to
--   apply the reverse mode of algorithmic differentiation and evaluation.

  type Circuit ( nbr : integer32 ) is record
    dim : integer32;                               -- dimension
    pdg : integer32;                               -- polynomial degree
    xps : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponents
    idx : Standard_Integer_VecVecs.VecVec(1..nbr); -- indices of exponents
    fac : Standard_Integer_VecVecs.VecVec(1..nbr); -- factor indices
    cff : DoblDobl_Complex_Vectors.Vector(1..nbr); -- coefficients
    cst : DoblDobl_Complex_Numbers.Complex_Number; -- constant
    fwd : DoblDobl_Complex_Vectors.Link_to_Vector; -- forward products
    bck : DoblDobl_Complex_Vectors.Link_to_Vector; -- backward products
    crs : DoblDobl_Complex_Vectors.Link_to_Vector; -- cross products
  end record;

  type Link_to_Circuit is access Circuit;

  type Circuits is array ( integer32 range <> ) of Link_to_Circuit;

-- A system stores the sequence of circuits for each polynomial,
-- along with work space and the final outcomes.

  type System ( neq,dim : integer32 ) is record
    crc : Circuits(1..neq);                        -- polynomials
    mxe : Standard_Integer_Vectors.Vector(1..dim); -- exponent maxima
    pwt : DoblDobl_Complex_VecVecs.VecVec(1..dim); -- power table
    yd : DoblDobl_Complex_Vectors.Link_to_Vector;  -- work space for a gradient
    fx : DoblDobl_Complex_Vectors.Vector(1..neq);  -- function value
    jm : DoblDobl_Complex_Matrices.Matrix(1..neq,1..dim); -- Jacobian matrix
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
  --   computes mxe and allocates space for a system.

  function Allocate ( neq,dim : integer32 )
                    return DoblDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Returns a vector of square matrices of dimension dim,
  --   the returned vector has range 1..neq.

-- ALGORITHMIC DIFFERENTIATION AND EVALUATION OF CIRCUITS :

  procedure EvalDiff
              ( s : in out System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure EvalDiff
              ( s : in Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in s at x.

  -- REQUIRED :
  --   All space for the power table and yd has been allocated.

  -- ON ENTRY :
  --   s        properly defined and allocated system of circuits;
  --   x        values for the variables in the system.

  -- ON RETURN :
  --   s.pwt    power table updated for the values in x;
  --   s.fx     function value of the circuits at x;
  --   s.jm     the Jacobian matrix evaluated at x.

  procedure EvalDiff2
              ( s : in out System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vh : in DoblDobl_Complex_VecMats.VecMat );
  procedure EvalDiff2
              ( s : in Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vh : in DoblDobl_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in s at x.

  -- REQUIRED :
  --   All space for the power table and yd has been allocated.

  -- ON ENTRY :
  --   s        properly defined and allocated system of circuits;
  --   x        values for the variables in the system;
  --   vh       space allocated for dim matrices,
  --            all matrices have 1..dim for range(1) and range(2).

  -- ON RETURN :
  --   s.pwt    power table updated for the values in x;
  --   s.fx     function value of the circuits at x;
  --   s.jm     the Jacobian matrix evaluated at x;
  --   vh       vector of evaluated Hessian matrices.

  procedure EvalDiff
              ( c : in Circuits;
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                pwt : in DoblDobl_Complex_VecVecs.VecVec;
                fx : out DoblDobl_Complex_Vectors.Vector;
                jm : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in c at x.

  -- ON ENTRY :
  --   c        a sequence of circuits, properly defined and allocated;
  --   x        a vector of values for the variables;
  --   yd       work space for the function value and gradient,
  --            of range 0..dim, where dim = x'last;
  --   pwt      power table defined and computed for x.

  -- ON RETURN :
  --   fx       vector of function values of the circuits at x;
  --   jm       matrix of partial derivatives.

  procedure EvalDiff2
              ( c : in Circuits;
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                pwt : in DoblDobl_Complex_VecVecs.VecVec;
                fx : out DoblDobl_Complex_Vectors.Vector;
                jm : out DoblDobl_Complex_Matrices.Matrix;
                vh : in DoblDobl_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in c at x,
  --   returns the evaluated Jacobian and vector of Hessians.

  -- ON ENTRY :
  --   c        a sequence of circuits, properly defined and allocated;
  --   x        a vector of values for the variables;
  --   yd       work space for the function value and gradient,
  --            of range 0..dim, where dim = x'last;
  --   pwt      power table defined and computed for x;
  --   vh       space allocated for dim matrices,
  --            all matrices have 1..dim for range(1) and range(2).

  -- ON RETURN :
  --   fx       vector of function values of the circuits at x;
  --   jm       matrix of partial derivatives;
  --   vh       vector of evaluated Hessian matrices.

-- SINGULAR VALUE DECOMPOSITIONS :

  procedure Singular_Values
              ( s : in out System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vh : in DoblDobl_Complex_VecMats.VecMat;
                U : out DoblDobl_Complex_Matrices.Matrix;
                V : out DoblDobl_Complex_Matrices.Matrix;
                e : out DoblDobl_Complex_Vectors.Vector;
                svls : in DoblDobl_Complex_VecVecs.VecVec );
  procedure Singular_Values
              ( s : in Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vh : in DoblDobl_Complex_VecMats.VecMat;
                U : out DoblDobl_Complex_Matrices.Matrix;
                V : out DoblDobl_Complex_Matrices.Matrix;
                e : out DoblDobl_Complex_Vectors.Vector;
                svls : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates and differentiations the circuits in s at x,
  --   computes the values of all Hessian matrices at x,
  --   computes all singular values of the Jacobian matrix
  --   and of all Hessian matrices.

  -- REQUIRED :
  --   All space for the power table and yd has been allocated.

  -- ON ENTRY :
  --   s        properly defined and allocated system of circuits;
  --   x        values for the variables in the system;
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
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                pwt : in DoblDobl_Complex_VecVecs.VecVec;
                A : out DoblDobl_Complex_Matrices.Matrix;
                U : out DoblDobl_Complex_Matrices.Matrix;
                V : out DoblDobl_Complex_Matrices.Matrix;
                e : out DoblDobl_Complex_Vectors.Vector;
                s : out DoblDobl_Complex_Vectors.Vector );
  procedure Singular_Values
              ( c : in Link_to_Circuit;
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                pwt : in DoblDobl_Complex_VecVecs.VecVec;
                A : out DoblDobl_Complex_Matrices.Matrix;
                U : out DoblDobl_Complex_Matrices.Matrix;
                V : out DoblDobl_Complex_Matrices.Matrix;
                e : out DoblDobl_Complex_Vectors.Vector;
                s : out DoblDobl_Complex_Vectors.Vector );

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
  --   x        a vector of values for the variables;
  --   yd       work space for the function value and gradient,
  --            of range 0..dim, where dim = x'last;
  --   pwt      power table defined and computed for x.

  -- ON RETURN :
  --   A        the Hessian matrix of c at x, 
  --            however, its values are destroyed by the SVD;
  --   U        the U matrix in the SVD of A;
  --   V        the V matrix in the SVD of A;
  --   e        contains error information on the SVD computation;
  --   s        the first c.dim entries are the singular values of A.

-- ALGORITHMIC DIFFERENTIATION AND EVALUATION OF ONE CIRCUIT :
--   The Indexed_Speel procedures are for circuits where the exponents
--   are either zero or one.  There are an important subclass to deal
--   with monomials that have no higher powers.
--   The general case is handled by the Speel procedures,
--   with wrappers working on circuits.
--   Both indexed and general speel procedure compute the gradient
--   and optionally, the evaluated Hessian matrix.

  procedure Indexed_Speel
               ( c : in Circuit;
                 x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                 h : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentatiates the circuit c at x,
  --   stores the function value at yd(0), the gradient at yd(x'range),
  --   and the Hessian at the matrix h.
  --   Wraps an Indexed_Speel procedure, using the c.xps as indices.
  --   This procedure is for monomials that are products of variables,
  --   with no exponent higher than one.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   x        vector of range 1..c.dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x;
  --   h        the Hessian matrix at x.

  procedure Indexed_Speel
              ( idx : in Standard_Integer_Vectors.Vector;
                cff : in DoblDobl_Complex_Numbers.Complex_Number;
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                fwd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                bck : in DoblDobl_Complex_Vectors.Link_to_Vector;
                crs : in DoblDobl_Complex_Vectors.Link_to_Vector;
                h : in out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates an indexed product, multiplied with a coefficient,
  --   computes its gradient and updates the Hessian matrix.
  --   Is called frequently by the next Indexed_Speel procedure.

  -- REQUIRED : idx'last >= 2.
  --   x'range = 1..dim and yd'range = 0..dim,
  --   fwd'range = 1..dim-1, bck'range = 1..dim-2 = crs'range,
  --   h'range(1) = h'range(2) = 1..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in the monomial;
  --   cff      coefficient of the monomial;
  --   x        vector of range 1..dim, with values for x;
  --   yd       vector of range 0..dim, allocated for the result.
  --   fwd      work space vector of range 1..dim-1,
  --            for the forward products;
  --   bck      work space vector of range 1..dim-2,
  --            for the backward products;
  --   crs      work space vector of range 1..dim-2,
  --            for the cross products;
  --   h        current values for the Hessian matrix,
  --            or initialized to zero.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x;
  --   h        the updated upper triangular Hessian matrix at x.

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                cff : in DoblDobl_Complex_Vectors.Vector;
                cst : in DoblDobl_Complex_Numbers.Complex_Number;
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                fwd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                bck : in DoblDobl_Complex_Vectors.Link_to_Vector;
                crs : in DoblDobl_Complex_Vectors.Link_to_Vector;
                h : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   sequence of indexed products of variables.

  -- REQUIRED :
  --   idx'range = cff'range and all vectors in idx have values
  --   in range 1..dim, where dim is the number of variables,
  --   x'range = 1..dim and yd'range = 0..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in each monomial;
  --   cff      coefficients of the monomials;
  --   cst      constant coefficient of the circuit;
  --   x        vector of range 1..dim, with values for x;
  --   yd       vector of range 0..dim, allocated for the result.
  --   fwd      work space vector of range 1..dim-1,
  --            for the forward products;
  --   bck      work space vector of range 1..dim-2,
  --            for the backward products;
  --   crs      work space vector of range 1..dim-2,
  --            for the cross products.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x;
  --   h        the Hessian matrix at x.

  procedure Indexed_Speel
              ( c : in Circuit;
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluates and differentiates the circuit c at x
  --   and stores the result in yd.
  --   Wraps the next Indexed_Speel procedure, using the c.xps as indices.
  --   This procedure is for monomials that are products of variables,
  --   with no exponent higher than one.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   x        vector of range 1..c.dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x.

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                cff : in DoblDobl_Complex_Vectors.Vector;
                cst : in DoblDobl_Complex_Numbers.Complex_Number;
                x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                fwd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                bck : in DoblDobl_Complex_Vectors.Link_to_Vector;
                crs : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   indexed sequence of products.

  -- REQUIRED :
  --   idx'range = cff'range and all vectors in idx have values
  --   in range 1..dim, where dim is the number of variables,
  --   x'range = 1..dim and yd'range = 0..dim.

  -- ON ENTRY :
  --   idx      indices to participating variables in each monomial;
  --   cff      coefficients of the monomials;
  --   cst      constant coefficient of the circuit;
  --   x        vector of range 1..dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result.
  --   fwd      work space vector of range 1..dim-1,
  --            for the forward products;
  --   bck      work space vector of range 1..dim-2,
  --            for the backward products;
  --   crs      work space vector of range 1..dim-2,
  --            for the cross products.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x.

  procedure Speel ( c : in Circuit;
                    x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates and differentiates the circuit c at x
  --   and stores the result in yd.
  --   Wraps the next Speel procedure, using the data in c.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   x        vector of range 1..c.dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result;
  --   pwt      power table to compute the common factors.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x.

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    cff : in DoblDobl_Complex_Vectors.Vector;
                    cst : in DoblDobl_Complex_Numbers.Complex_Number;
                    x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    fwd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    bck : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    crs : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   indexed sequence of products, with higher powers factored out.

  -- REQUIRED :
  --   idx'range = cff'range and all vectors in idx have values
  --   in range 1..dim, where dim is the number of variables,
  --   x'range = 1..dim and yd'range = 0..dim.

  -- ON ENTRY :
  --   xps      exponent vectors of the monomials in the circuit;
  --   idx      indices to participating variables in each monomial;
  --   fac      factor indices of the exponents;
  --   cff      coefficients of the monomials;
  --   cst      constant coefficient of the circuit;
  --   x        vector of range 1..dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result.
  --   fwd      work space vector of range 1..dim-1,
  --            for the forward products;
  --   bck      work space vector of range 1..dim-2,
  --            for the backward products;
  --   crs      work space vector of range 1..dim-2,
  --            for the cross products;
  --   pwt      power table to compute the common factors.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x.

  procedure Speel ( c : in Circuit;
                    x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in DoblDobl_Complex_VecVecs.VecVec;
                    h : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the circuit c at x,
  --   stores the result in yd, and computes the gradient.
  --   Wraps one of the next Speel procedures, using the data in c.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   x        vector of range 1..c.dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result;
  --   pwt      power table to compute the common factors.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x;
  --   h        the Hessian matrix evaluated at x.

  procedure Speel ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                    idx : in Standard_Integer_Vectors.Vector;
                    cff : in DoblDobl_Complex_Numbers.Complex_Number;
                    x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    fwd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    bck : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    crs : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Updates the function value and gradient in yd,
  --   for a general term in the circuit.
  --   This is a helper procedure in the next Speel procedure.

  -- REQUIRED : fac /= null.

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    cff : in DoblDobl_Complex_Vectors.Vector;
                    cst : in DoblDobl_Complex_Numbers.Complex_Number;
                    x,yd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    fwd : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    bck : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    crs : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in DoblDobl_Complex_VecVecs.VecVec;
                    h : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Runs the reverse mode of algorithmic differentiation on an
  --   indexed sequence of products, with higher powers factored out.

  -- REQUIRED :
  --   idx'range = cff'range and all vectors in idx have values
  --   in range 1..dim, where dim is the number of variables,
  --   x'range = 1..dim and yd'range = 0..dim.

  -- ON ENTRY :
  --   xps      exponent vectors of the monomials in the circuit;
  --   idx      indices to participating variables in each monomial;
  --   fac      factor indices of the exponents;
  --   cff      coefficients of the monomials;
  --   cst      constant coefficient of the circuit;
  --   x        vector of range 1..dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result.
  --   fwd      work space vector of range 1..dim-1,
  --            for the forward products;
  --   bck      work space vector of range 1..dim-2,
  --            for the backward products;
  --   crs      work space vector of range 1..dim-2,
  --            for the cross products;
  --   pwt      power table to compute the common factors.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x;
  --   h        the Hessian matrix evaluated at x.

-- AUXILIARY PROCEDURES :

  procedure Forward ( x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                      f : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f, in complex arithmetic.
  --   This procedure is for testing the next Forward procedure.

  -- REQUIRED : f'first = x'first = 1 and f'last >= x'last-1.

  procedure Forward_Backward
              ( x,f,b : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.

  -- REQUIRED :
  --    f'first = x'first = 1 and f'last >= x'last-1,
  --    b'first = b'first = 1 and b'last >= x'last-2.

  procedure Fused_Forward_Backward
              ( x,f,b : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.
  --   The two loops are fused, resulting in better performance.

  -- REQUIRED :
  --    f'first = x'first = 1 and f'last >= x'last-1.
  --    b'first = b'first = 1 and b'last >= x'last-2.

  procedure Forward_Backward_Cross
              ( x,f,b,c : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                x,f,b,c : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.
  --   Computes all cross products of the values in x
  --   and stores the products in b.

  -- REQUIRED : x'last > 2,
  --   f'first = x'first = 1 and f'last >= x'last-1, or >= idx'last-1,
  --   b'first = b'first = 1 and b'last >= x'last-2, or >= idx'last-2,
  --   c'first = c'first = 1 and c'last >= x'last-2, or >= idx'last-2.

  -- ON ENTRY :
  --   idx      if provided, then only those values of x
  --            as indexed by the entries in idx will be used;
  --   x        values for the variables;
  --   f        space allocated for forward products;
  --   b        space allocated for backward products;
  --   c        space allocated for cross products.

  -- ON RETURN : let n be x'last-1, or idx'last-1
  --   f(n)     holds the product of all (or those in idx) variables in x;
  --   f(n-1)   is the partial derivative of the product
  --            with respect to the last variable (in idx);
  --   b(n-2)   is the partial derivative of the product
  --            with respect to the first variable (or idx(1));
  --   c(k)     is the partial derivative of the product
  --            with respect to the (k+1)-th variable (or idx(k+1)).

  procedure Fused_Forward_Backward_Cross
              ( x,f,b,c : in DoblDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.
  --   Computes all cross products of the values in x
  --   and stores the products in b.
  --   Applies loop fusion.

  -- REQUIRED : x'last > 2, 
  --   f'first = x'first = 1 and f'last >= x'last-1,
  --   b'first = b'first = 1 and b'last >= x'last-2,
  --   c'first = c'first = 1 and c'last >= x'last-2.

  -- ON ENTRY :
  --   x        values for the variables;
  --   f        space allocated for forward products;
  --   b        space allocated for backward products;
  --   c        space allocated for cross products.

  -- ON RETURN : let n be x'last-1,
  --   f(n)     holds the product of all variables in x;
  --   f(n-1)   is the partial derivative of the product
  --            with respect to the last variable;
  --   b(n-2)   is the partial derivative of the product
  --            with respect to the first variable;
  --   c(k)     is the partial derivative of the product
  --            with respect to the (k+1)-th variable.

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of range mxe'range with space for
  --   complex vectors of range 1..mxe(k)-1, for k in mxe'range.

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                pwt : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the power table for the values of the variables in x,
  --   in complex arithmetic.
  --   This procedure is for testing the next Power_Table procedure.

  -- REQUIRED :
  --   mxe'range = x'range, pwt is allocated according to mxe,
  --   pwt'range = x'range and pwt(k)'range = 1..mxe(k)-1.

  -- ON ENTRY :
  --   mxe      highest exponents of the variables,
  --            mxe(k) is the highest exponent of the k-th variable
  --   x        values for all variables;
  --   pwt      allocated memory for all powers of the values in x.

  -- ON RETURN :
  --   pwt      power table, pwt(k)(i) equals x(k)**(i+1),
  --            for i in range 1..mxe(k)-1.

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                cff : in DoblDobl_Complex_Numbers.Complex_Number;
                pwt : in DoblDobl_Complex_VecVecs.VecVec;
                res : out DoblDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Given exponents in xps and the factor indices in fac,
  --   computes the value of the common factor at x,
  --   multiplied with the coefficient cff,
  --   and with the aid of the power table pwt.

  -- ON ENTRY :
  --   xps      values of the exponents for the powers of x;
  --   fac      factor indices;
  --   x        values for the variables used for low powers;
  --   cff      coefficient of the monomial;
  --   pwt      power table for the higher powers of x,
  --            pwt(k)(i) stores x(k)**(i+1).

  -- ON RETURN :
  --   res      the coefficient multiplied with the common factor.

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit );
  procedure Clear ( c : in out Link_to_Circuit );
  procedure Clear ( c : in out Circuits );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the circuit.

  procedure Clear ( s : in out System );
  procedure Clear ( s : in out Link_to_System );
  procedure Clear ( s : in out System_Array );

  -- DESCRIPTION :
  --   Deallocates the space occupied by s.

end DoblDobl_Complex_Circuits;
