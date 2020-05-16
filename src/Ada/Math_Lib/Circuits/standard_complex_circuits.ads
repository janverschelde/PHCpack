with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Standard_Complex_Circuits is

-- DESCRIPTION :
--   This package contains development code for standard_coefficient_circuits,
--   to test the correctness of the better performing algorithms for
--   algorithmic differentiation and evaluation.

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

end Standard_Complex_Circuits;
