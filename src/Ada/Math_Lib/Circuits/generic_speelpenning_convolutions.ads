with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Abstract_Ring;
with Generic_Vectors;
with Generic_VecVecs;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package VecVecs is new Generic_VecVecs(Ring,Vectors);

package Generic_Speelpenning_Convolutions is

-- DESCRIPTION :
--   The reverse mode of algorithmic differentiation computes
--   the value of a product and all its partial derivatives.
--   This package offers a vectorized version on the coefficients
--   of power series, all truncated to the same fixed degree.

  type VecVecVec is array ( integer32 range <> ) of VecVecs.Link_to_VecVec;
  -- A three dimensional structure to store the coefficient vectors
  -- of powers of series.
 
  type Link_to_VecVecVec is access VecVecVec; -- stores the power table

  function Create ( x : VecVecs.VecVec;
                    d : Standard_Integer_Vectors.Vector )
                  return Link_to_VecVecVec;

  -- DESCRIPTION :
  --   Stores all powers x(i)^k for k ranging from 2 to d(i),
  --   for i in d'range = x'range, in the power table.
  --   The i-th entry in the power table contains the powers of x(i),
  --   if d(i) > 1, starting with x(i)^2 at the first position.

  procedure Clear ( pwt : in out Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the power table pwt.

  function Allocate_Coefficients
             ( deg : integer32 ) return Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of the series,
  --   truncaged to degree deg and initialized to zero.

  function Allocate_Coefficients
             ( dim,deg : integer32 ) return VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of the series
  --   truncated to degree deg and initialized to zero.
  --   The vector on return has range 1..dim.

  procedure Update ( values : in Vectors.Link_to_Vector;
                     inc : in Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds the elements in inc to the values.

  -- REQUIRED : values'range = inc'range.

  procedure Multiply ( first,second,product : in Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Convolutes the vectors first and second into the product,
  --   corresponding to the multiplication of two series of the same degree.

  -- REQUIRED : first'last = second'last = product'last.

  procedure Speel ( x : in VecVecs.VecVec;
                    forward,backward,cross : in out VecVecs.VecVec );
  procedure Speel ( x : in VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    forward,backward,cross : in out VecVecs.VecVec );

  -- DESCRIPTION :
  --   Inline computation of the coefficients of the product
  --   of the series in x, which are all of the same degree.

  -- REQUIRED :
  --   The range of x starts at 1 and ends at 2 or a larger value.
  --   All vectors in x, forward, backward, and cross are allocated
  --   and of the range 0..degree, for some same fixed degree.

  -- ON ENTRY :
  --   x            a vector with range starting at 1, ending at 2 or higher;
  --   idx          if provided, then only those indices of x with values
  --                in idx participate and dim = idx'last,
  --                otherwise, all values of x participate and dim = x'last;
  --   forward      has space allocated for dim-1 coefficient vectors,
  --                for series of the same fixed degree;
  --   backward     has space reserved for dim-2 coefficient vectors,
  --                for series of the same fixed degree;
  --   cross        has space reserved for dim-2 coefficient vectors,
  --                for series of the same fixed degree.

  -- ON RETURN :
  --   forward      accumulates the forward products,
  --                forward(dim-1) holds the coefficients for the product,
  --                forward(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   backward     accumulates the backward products,
  --                backward(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   cross        stores the cross products, cross(k) contains the
  --                coefficients of the partial derivative w.r.t. k+1.

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    x : in VecVecs.VecVec;
                    forward,backward,cross,yd : in out VecVecs.VecVec );
  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    cff : in VecVecs.VecVec; x : in VecVecs.VecVec;
                    forward,backward,cross,yd : in out VecVecs.VecVec;
                    wrk : in Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluation and differentiation of the sum of products,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   idx          indexed representation of a sum of products of variables;
  --   cff          coefficients of the products, if not provided,
  --                then all coefficients are considered as one;
  --   x            coefficient vectors of power series of same degree;
  --   forward      work space allocated for x'last-1 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   backward     work space allocated for x'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   cross        work space allocated for x'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   yd           vector of range 0..x'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   wrk          work space for the coefficients of the same fixed degree.

  -- ON RETURN :
  --   yd           yd(x'last+1) contains the coefficient vector of the value
  --                of the sum of products, evaluated at x,
  --                yd(k) is the k-th partial derivative at x.

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    cff : in VecVecs.VecVec; x : in VecVecs.VecVec;
                    forward,backward,cross,yd : in out VecVecs.VecVec;
                    wrk,acc : in Vectors.Link_to_Vector;
                    pwt : in Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Evaluation and differentiation of a polynomial,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   xps          exponent vector of the monomials;
  --   idx          indexed representation of the variables in exponents;
  --   fac          factor index of the exponents;
  --   cff          coefficients of the products, if not provided,
  --                then all coefficients are considered as one;
  --   x            coefficient vectors of power series of same degree;
  --   forward      work space allocated for x'last-1 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   backward     work space allocated for x'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   cross        work space allocated for x'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   yd           vector of range 0..x'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   wrk          work space for the coefficients of the same fixed degree;
  --   acc          work space for the coefficients of the same fixed degree;
  --   pwt          power table of the values in x.

  -- ON RETURN :
  --   yd           yd(x'last+1) contains the coefficient vector of the value
  --                of the sum of products, evaluated at x,
  --                yd(k) is the k-th partial derivative at x.

end Generic_Speelpenning_Convolutions;
