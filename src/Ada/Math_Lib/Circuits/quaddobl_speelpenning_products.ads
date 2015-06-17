with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;

package QuadDobl_Speelpenning_Products is

-- DESCRIPTION :
--   Elaboration of the example of Speelpenning to evaluate a product
--   of n variables along with all its derivatives, using vectors of
--   complex quad doubles.

  function Straight_Speel
             ( x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..x'last, if y = Speel(x),
  --   then y(0) contains the value x(1)*x(2)*..*x(x'last)
  --   and y(k) contains the k-th derivative of y(0) at x.
  --   This is the straightforward and wasteful evaluation,
  --   mimicking what happens if the monomial and all its derivatives
  --   are evaluated independently from each other.

  -- REQUIRED : n = x'last > 1.

  -- COST : n^2 - 2*n + 2*n - 1 = n^2 - 1 multiplications

  function Straight_Speel
             ( e : Standard_Natural_Vectors.Vector;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..x'last, if y = Speel(x),
  --   then y(0) contains the product of all x(k) for which e(k) /= 0
  --   and y(k) contains the k-th derivative of y(0) at x.
  --   This is the straightforward and wasteful evaluation,
  --   mimicking what happens if the monomial and all its derivatives
  --   are evaluated independently from each other.

  -- REQUIRED : n = e'last = x'last > 1.

  function Reverse_Speel
             ( x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..x'last, if y = Speel(x),
  --   then y(0) contains the value x(1)*x(2)*..*x(x'last)
  --   and y(k) contains the k-th derivative of y(0) at x.
  --   This is the more efficient evaluation in reverse mode.

  -- REQUIRED : n > 1.

  -- COST : n-1 + n-2 + n-2  = 3*n - 5 multiplications,
  --   with n-1 additional storage for intermediate results.

  procedure Nonzeroes
             ( e : in Standard_Natural_Vectors.Vector;
               x : in QuadDobl_Complex_Vectors.Vector;
               idx : out Standard_Integer_Vectors.Vector;
               enz : out Standard_Natural_Vectors.Vector;
               xnz : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Extracts the nz nonzero elements of e into nze
  --   and the corresponding entries of x into nzx.

  -- REQUIRED : enz'range = xnz'range = idx'range = 1..nz,
  --   where nz equals the number of nonzero elements of e.

  -- ON ENTRY :
  --   e       exponent vector of a monomial may contain zeroes;
  --   x       values where to evaluate the monomial at.

  -- ON RETURN :
  --   idx     indices of the nonzero elements in e;
  --   enz     nonzero elements in e;
  --   xnz     values of x corresponding to nonzero exponents of e.

  function Indexed_Speel
             ( nx,nz : integer32;
               indnz : Standard_Integer_Vectors.Vector;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..nx = x'last with
  --   at position 0 the product of the values in x; and
  --   at position i the i-th derivative of the product.

  -- ON ENTRY :
  --   nx      the number of elements in x, x is of range 1..nx;
  --   nz      the number of participating indices (nonzero exponents);
  --   indnz   has range 1..nz indicates the participating variables
  --           in the product, those that occur with nonzero exponent;
  --   x       values for all elements, of range 1..nx.

  function Reverse_Speel
             ( e : Standard_Natural_Vectors.Vector;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..x'last, if y = Reverse_Speel(x),
  --   then y(0) contains the product of all x(k) if e(k) > 0
  --   and y(k) contains the k-th derivative of y(0) at x.
  --   This is the more efficient evaluation in reverse mode.

  -- REQUIRED : n = e'last = x'last > 1, and nz > 1
  --   where nz equals the number of nonzero elements in e.

  -- COST : at most nz-1 + nz-2 + nz-2  = 3*nz - 5 multiplications,
  --   with 2*nz-1 additional storage for intermediate results.

  function Indexed_Reverse_Speel
             ( idx : Standard_Integer_Vectors.Vector;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DECRIPTION :
  --   Returns a vector y of range 0..x'last where y(0) contains
  --   the product of the values in x as indexed by idx.
  --   The index vector idx defines which variables in x
  --   participate to the product of the variables, computed as
  --   idx = Standard_Speelpenning_Products.Nonzero_Indices(e),
  --   for an exponent vector that has the same range as x.
  --   The derivatives are in y(idx(i)) for i in idx'range.

end QuadDobl_Speelpenning_Products;
