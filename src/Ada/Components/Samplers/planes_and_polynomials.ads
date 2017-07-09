with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;

package Planes_and_Polynomials is

-- DESCRIPTION :
--   Provides operations on polynomials in connection with planes.

-- CONVERTORS :

  function Hyperplane ( cff : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Polynomials.Poly;
  function Hyperplane ( cff : DoblDobl_Complex_Vectors.Vector )
                      return DoblDobl_Complex_Polynomials.Poly;
  function Hyperplane ( cff : QuadDobl_Complex_Vectors.Vector )
                      return QuadDobl_Complex_Polynomials.Poly;
  function Hyperplane ( cff : Multprec_Complex_Vectors.Vector )
                      return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial representation of the hyperplane
  --     cff(0) + cff(1)*x1 + cff(2)*x2 + .. + cff(n)*xn.

  function Hyperplane ( cff : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Laurentials.Poly;
  function Hyperplane ( cff : DoblDobl_Complex_Vectors.Vector )
                      return DoblDobl_Complex_Laurentials.Poly;
  function Hyperplane ( cff : QuadDobl_Complex_Vectors.Vector )
                      return QuadDobl_Complex_Laurentials.Poly;
  function Hyperplane ( cff : Multprec_Complex_Vectors.Vector )
                      return Multprec_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns the representation of the hyperplane
  --     cff(0) + cff(1)*x1 + cff(2)*x2 + .. + cff(n)*xn,
  --   where n = cff'last, as a Laurent polynomial in n variables.

  function Hyperplane ( cff : Standard_Complex_Vectors.Vector;
                        tol : double_float )
                      return Standard_Complex_Polynomials.Poly;
  function Hyperplane ( cff : Multprec_Complex_Vectors.Vector;
                        tol : double_float )
                      return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial representation of the hyperplane
  --     cff(0) + cff(1)*x1 + cff(2)*x2 + .. + cff(n)*xn,
  --   ignoring all coefficients in magnitude less than tol.

  function Polynomial ( p : Standard_Complex_Polynomials.Poly )
                      return Standard_Complex_Vectors.Vector;
  function Polynomial ( p : Standard_Complex_Laurentials.Poly )
                      return Standard_Complex_Vectors.Vector;
  function Polynomial ( p : DoblDobl_Complex_Polynomials.Poly )
                      return DoblDobl_Complex_Vectors.Vector;
  function Polynomial ( p : DoblDobl_Complex_Laurentials.Poly )
                      return DoblDobl_Complex_Vectors.Vector;
  function Polynomial ( p : QuadDobl_Complex_Polynomials.Poly )
                      return QuadDobl_Complex_Vectors.Vector;
  function Polynomial ( p : QuadDobl_Complex_Laurentials.Poly )
                      return QuadDobl_Complex_Vectors.Vector;
  function Polynomial ( p : Multprec_Complex_Polynomials.Poly )
                      return Multprec_Complex_Vectors.Vector;
  function Polynomial ( p : Multprec_Complex_Laurentials.Poly )
                      return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector representation of the hyperplane in p.

-- ELIMINATORS FOR POLYNOMIALS :

  function Remove_Variable
             ( p : Standard_Complex_Polynomials.Poly; k : integer32 )
             return Standard_Complex_Polynomials.Poly;

  function Remove_Variable
             ( p : Multprec_Complex_Polynomials.Poly; k : integer32 )
             return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Removes the k-th entry out of all exponents in the polynomial.

  function Substituting_Polynomial
             ( nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Standard_Complex_Vectors.Vector; tol : double_float )
             return Standard_Complex_Polynomials.Poly;

  function Substituting_Polynomial
             ( nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Multprec_Complex_Vectors.Vector; tol : double_float )
             return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a polynomial in nvar variables that expresses the
  --   nonpivot variable as a polynomial in the pivot variables.
  --   The coefficients of the hyperplane are in hypcff in the format:
  --     hypcff(0) + hypcff(1)*x(1) + .. + hyp(n)*x(n) = 0
  --   and it is assumed the nonpivot variables do no interfere.
  --   Coefficients in absolute value smaller than tol are ignored.

  function Substituting_Polynomials
             ( nequ,nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Standard_Complex_VecVecs.VecVec; tol : double_float )
             return Standard_Complex_Poly_Systems.Poly_Sys;

  function Substituting_Polynomials
             ( nequ,nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Multprec_Complex_VecVecs.VecVec; tol : double_float )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns nequ polynomials in nvar variables that express each
  --   nonpivot variable as a polynomial in the pivot variables.

  function Restrict_to_Linear_Space
             ( p : Standard_Complex_Polynomials.Poly; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               subpols : Standard_Complex_Poly_Systems.Poly_Sys;
               tol : double_float ) return Standard_Complex_Polynomials.Poly;

  function Restrict_to_Linear_Space
             ( p : Multprec_Complex_Polynomials.Poly; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               subpols : Multprec_Complex_Poly_Systems.Poly_Sys;
               tol : double_float ) return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Restricts the polynomial to the linear space using the
  --   polynomials in subpols that express the nonpivot variables
  --   in terms of those indexed by the pivots.
  --   The last variables whose number equals level remain present.

  function Restrict_to_Linear_Space
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               hypcff : Standard_Complex_VecVecs.VecVec; tol : double_float )
             return Standard_Complex_Poly_Systems.Poly_Sys;

  function Restrict_to_Linear_Space
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               hypcff : Multprec_Complex_VecVecs.VecVec; tol : double_float )
             return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Restricts the polynomials to the linear space defined by hypcff
  --     hypcff(0) + hypcff(1)*x(1) + .. + hypcff(n)*x(n) = 0
  --   using the planes to remove of all nonpivot variables.
  --   The last variables whose number equals level remain present.

  -- REQUIRED : pivots'length + hypcff'length = Number_of_Unknowns(p).
  --   The coefficient matrix of the hyperplanes is diagonalized with
  --   respect to the nonpivot variables, which means that the nonpivot
  --   variables do not depend on each other.  The coefficients with
  --   the nonpivot variables are either zero or one.

-- ELIMINATORS FOR SOLUTIONS :

  function Remove_Variables
             ( v : Standard_Complex_Vectors.Vector; level,newdim : integer32;
               pivots : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;

  function Remove_Variables
             ( v : Multprec_Complex_Vectors.Vector; level,newdim : integer32;
               pivots : Standard_Integer_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..newdim that has only those entries from
  --   v as indicated in pivots and the last level number of values.

  function Remove_Variables
             ( s : Standard_Complex_Solutions.Solution;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Solutions.Solution;

  function Remove_Variables
             ( s : Multprec_Complex_Solutions.Solution;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Multprec_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Removes all nonpivot variables from the solution, but keeps the
  --   values for the added variables, whose number is indicated by level.

  function Remove_Variables
             ( sols : Standard_Complex_Solutions.Solution_List;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Solutions.Solution_List;

  function Remove_Variables
             ( sols : Multprec_Complex_Solutions.Solution_List;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Removes all nonpivot variable from all solutions in the list.

  function Restrict_Solution
             ( sols : Standard_Complex_Solutions.Solution_List;
               ind,level : integer32;
               pivots : Standard_Integer_Vectors.Vector ) 
             return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the same list of solutions as sols, except for the
  --   solution at position ind, which is restricted to the subspace.
  --   There is sharing between the lists.

end Planes_and_Polynomials;
