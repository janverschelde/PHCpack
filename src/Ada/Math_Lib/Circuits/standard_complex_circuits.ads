with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Standard_Complex_Circuits is

-- DESCRIPTION :
--   This package contains development code for standard_coefficient_circuits,
--   to test the correctness of the better performing algorithms for
--   algorithmic differentiation and evaluation.

-- DATA STRUCTURES :
--   A circuit stores the exponents and coefficients and hold work space to
--   apply the reverse mode of algorithmic differentiation and evaluation.

  type Circuit ( nbr : integer32 ) is record
    dim : integer32;                               -- dimension
    xps : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponents
    idx : Standard_Integer_VecVecs.VecVec(1..nbr); -- indices of exponents
    fac : Standard_Integer_VecVecs.VecVec(1..nbr); -- factor indices
    cff : Standard_Complex_Vectors.Vector(1..nbr); -- coefficients
    cst : Standard_Complex_Numbers.Complex_Number; -- constant
    fwd : Standard_Complex_Vectors.Link_to_Vector; -- forward products
    bck : Standard_Complex_Vectors.Link_to_Vector; -- backward products
    crs : Standard_Complex_Vectors.Link_to_Vector; -- cross products
  end record;

  type Link_to_Circuit is access Circuit;

  type Circuits is array ( integer range <> ) of Link_to_Circuit;

  function Allocate ( nbr,dim : integer32 ) return Circuit;

  -- DESCRIPTION :
  --   Returns a circuit for a polynomial with nbr monomials,
  --   with dim variables, and with allocated work space vectors.

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF ONE CIRCUIT :

  procedure Speel ( c : in Circuit;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluates and differentiates the circuit c at x
  --   and stores the result in yd.
  --   Wraps the next Speel procedure, using the c.xps as indices.
  --   The results are correct if all monomials are products of variables,
  --   with no exponent higher than one.

  -- ON ENTRY :
  --   c        circuit properly defined and with allocated workspace;
  --   x        vector of range 1..c.dim, with values for x;
  --   yd       vector of range 0..c.dim, allocated for the result.

  -- ON RETURN :
  --   yd(0)    the value of the circuit at x;
  --   yd(k)    the k-th derivative of the circuit at x.

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    cff : in Standard_Complex_Vectors.Vector;
                    cst : in Standard_Complex_Numbers.Complex_Number;
                    x,yd : in Standard_Complex_Vectors.Link_to_Vector;
                    fwd : in Standard_Complex_Vectors.Link_to_Vector;
                    bck : in Standard_Complex_Vectors.Link_to_Vector;
                    crs : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Implements the Speel on the circuit c.

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

-- AUXILIARY PROCEDURES :

  procedure Forward ( x : in Standard_Complex_Vectors.Link_to_Vector;
                      f : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f, in complex arithmetic.
  --   This procedure is for testing the next Forward procedure.

  -- REQUIRED : f'first = x'first = 1 and f'last >= x'last-1.

  procedure Forward_Backward
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.

  -- REQUIRED :
  --    f'first = x'first = 1 and f'last >= x'last-1,
  --    b'first = b'first = 1 and b'last >= x'last-2.

  procedure Fused_Forward_Backward
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector );

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
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector );
  procedure Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector );

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
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector );

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
             return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of range mxe'range with space for
  --   complex vectors of range 1..mxe(k)-1, for k in mxe'range.

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                pwt : in Standard_Complex_VecVecs.VecVec );

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
                x : in Standard_Complex_Vectors.Link_to_Vector;
                cff : in Standard_Complex_Numbers.Complex_Number;
                pwt : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Numbers.Complex_Number );

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

  -- DESCRIPION :
  --   Deallocates the space occupied by the circuit.

end Standard_Complex_Circuits;
