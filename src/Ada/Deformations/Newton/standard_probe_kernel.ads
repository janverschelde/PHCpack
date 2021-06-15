with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors; 
with Standard_Complex_Matrices;         use Standard_Complex_Matrices; 
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;    
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;

package Standard_Probe_Kernel is

-- DESCRIPTION :
--   This package provides primitive operations used to determine
--   the minimal deflation order needed to recondition a singular root.

  function Maximal_Degree ( p : Poly_Sys ) return integer32;

  -- DESCRIPTION :
  --   Returns the maximal degree of the polynomials in p,
  --   as the natural upper bound for the order of the deflation.

  function Random_Vector_in_Kernel
              ( V : Matrix; corank : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns a random vector in the kernel of a matrix,
  --   given its numerical corank and the V matrix from its SVD. 

  -- ON ENTRY :
  --   V        the V matrix from the SVD output;
  --   corank   numerical corank of a matrix.

  -- ON RETURN :
  --   a random combination of the last corank columns of V.

  function Sample_Sum_on_Line
              ( f : Eval_Poly_Sys; z,w,t : Vector ) return Vector;

  -- DESCRIPTION :
  --   Evaluates the line z+t*w at the sum of the polynomials in f.

  -- ON ENTRY :
  --   f        a system of polynomial functions;
  --   z        the offset vector of the line;
  --   w        a vector in the kernel as the direction of the line;
  --   t        values of the parameter for points on the line.

  -- ON RETURN :
  --   the sum of the polynomials in f at z+t(i)*w, for i in t'range;
  --   i.e.: the vector on return has range t'range.

  function Interpolation_Coefficients ( x,y : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the polynomial interpolating y at x.
  
  -- ON ENTRY :
  --   x        x-coordinates of the interpolation points, range 0..d;
  --   y        y-coordinates of the interpolation points, range 0..d.

  -- ON RETURN :
  --   coefficients of the interpolating polynomial of degree d,
  --   in a vector of range 0..d.

  function Numerical_Order ( c : Vector; tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the first element in c whose magnitude
  --   is larger than the given tolerance.
  --   If all elements in c are in magnitude smaller than tol,
  --   then c'last+1 is returned.

end Standard_Probe_Kernel;
