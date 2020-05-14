with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Standard_Coefficient_Circuits is

-- DESCRIPTION :
--   Separating real parts from imaginary parts, inlining the algorithms
--   for complex multiplication, fusing loops, results in faster code.

  procedure Forward ( x : in Standard_Complex_Vectors.Link_to_Vector;
                      f : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f, in complex arithmetic.
  --   This procedure is for testing the next Forward procedure.

  -- REQUIRED : f'first = x'first = 1 and f'last >= x'last-1.

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
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.
  --   Computes all cross products of the values in x
  --   and stores the products in b.

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

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of range mxe'range with space for
  --   complex vectors of range 1..mxe(k)-1, for k in mxe'range.

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of range mxe'range with space for
  --   floating-point vectors of range 1..mxe(k)-1, for k in mxe'range.

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

end Standard_Coefficient_Circuits;
