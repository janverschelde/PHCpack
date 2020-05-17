with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Integer_VecVecs;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;

package Standard_Coefficient_Circuits is

-- DESCRIPTION :
--   Separating real parts from imaginary parts, inlining the algorithms
--   for complex multiplication, fusing loops, results in faster code.

-- DATA STRUCTURES :
--   A circuit stores the exponents and coefficients and hold work space to
--   apply the reverse mode of algorithmic differentiation and evaluation.

  type Circuit ( nbr : integer32 ) is record
    xps : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponents
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

  type Circuits is array ( integer range <> ) of Link_to_Circuit;

  function Allocate ( nbr,dim : integer32 ) return Circuit;

  -- DESCRIPTION :
  --   Returns a circuits for a polynomial with nbr monomials,
  --   with dim variables, and with allocated work space vectors.

-- PROCEDURES :

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

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit );
  procedure Clear ( c : in out Link_to_Circuit );
  procedure Clear ( c : in out Circuits );

  -- DESCRIPION :
  --   Deallocates the space occupied by the circuit.

end Standard_Coefficient_Circuits;
