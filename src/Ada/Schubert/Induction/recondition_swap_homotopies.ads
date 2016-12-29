with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Matrices;

package Recondition_Swap_Homotopies is

-- DESCRIPTION :
--   To prevent columns in the moving coordinates of a k-plane from
--   drifting to zero in the swap cases of the Littlewood-Richardson
--   homotopies, an addition coordinate is added at the place of the
--   one in column s+1.  The solutions that fit in the k-plane are
--   then extended with an extra coordinate and the solution k-plane
--   is then reconditioned so that x(r+1,s+1) equals one.
--   One linear equations with random coefficients is added to the
--   system to select one affine patch in the homogenized problem.

  procedure Insert_One_Variable
              ( k : in integer32;
                t : in out Standard_Complex_Polynomials.Term );

  -- DESCRIPTION :
  --   Inserts one variable at position k into the term t,
  --   shifting all variables at position k and higher one
  --   position further.

  -- REQUIRED : k <= t.dg'last.

  procedure Insert_One_Variable
              ( k : in integer32;
                p : in out Standard_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Applies Insert_One_Variable(k,t) to all terms t of p.

  -- REQUIRED : k <= Number_of_Unknowns(p).

  procedure Insert_One_Variable
              ( k : in integer32;
                x : in out Standard_Complex_Poly_Matrices.Matrix );

  -- DESCRIPTION :
  --   Inserts one extra variable to each non null entry in X,
  --   at position k, shifting variables at position k and higher
  --   one position further in each exponent vector.
  --   The exponent of the extra variable equals zero.

  -- REQUIRED : k <= number of variables in each non null entry of X.

  procedure Insert_Variable_Pivot
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                i,j,k : in integer32 );

  -- DESCRIPTION :
  --   Inserts the variable at position k in x(i,j),
  --   with exponent equal to one.

  -- REQUIRED : Insert_One_Variable(k,x) has been executed.

  procedure Recondition
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                dim,s : in integer32 );

  -- DESCRIPTION :
  --   Prepares the moving plane x for reconditioning.

  -- ON ENTRY :
  --   x        symbolic form as defined by Swap_Moving_Plane
  --            in the package Checker_Homotopies;
  --   locmap   the localization map for the coordinates;
  --   dim      number of free variables computed from the localization map;
  --   s        index of the swap column.

  -- ON RETURN :
  --   x        to all non null polynomials in x, one variable has been
  --            added at position dim+1 and in the row of the pivot in
  --            column s+1 an extra variable has been added.

end Recondition_Swap_Homotopies;
