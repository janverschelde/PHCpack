with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

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
--   Every procedure/function appears in three precisions,
--   for double, double double, and quad double.

  procedure Insert_One_Variable
              ( t : in out Standard_Complex_Polynomials.Term;
                k : in integer32 );
  procedure Insert_One_Variable
              ( t : in out DoblDobl_Complex_Polynomials.Term;
                k : in integer32 );
  procedure Insert_One_Variable
              ( t : in out QuadDobl_Complex_Polynomials.Term;
                k : in integer32 );

  -- DESCRIPTION :
  --   Inserts one variable at position k into the term t,
  --   shifting all variables at position k and higher one
  --   position further.

  -- REQUIRED : k <= t.dg'last.

  procedure Remove_One_Variable
              ( t : in out Standard_Complex_Polynomials.Term;
                k : in integer32 );
  procedure Remove_One_Variable
              ( t : in out DoblDobl_Complex_Polynomials.Term;
                k : in integer32 );
  procedure Remove_One_Variable
              ( t : in out QuadDobl_Complex_Polynomials.Term;
                k : in integer32 );

  -- DESCRIPTION :
  --   Removes the k-th variable from the term t.
  --   This is the inverse of Insert_One_Variable(t,k).

  procedure Insert_One_Variable
              ( p : in out Standard_Complex_Polynomials.Poly;
                k : in integer32 );
  procedure Insert_One_Variable
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                k : in integer32 );
  procedure Insert_One_Variable
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                k : in integer32 );

  -- DESCRIPTION :
  --   Applies Insert_One_Variable(k,t) to all terms t of p.

  -- REQUIRED : k <= Number_of_Unknowns(p).

  procedure Remove_One_Variable
              ( p : in out Standard_Complex_Polynomials.Poly;
                k : in integer32 );
  procedure Remove_One_Variable
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                k : in integer32 );
  procedure Remove_One_Variable
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                k : in integer32 );

  -- DESCRIPTION :
  --   Removes the k-th variable from the polynomal p.
  --   This is the inverse of Insert_One_Variable(p,k).

  procedure Insert_One_Variable
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                k : in integer32 );
  procedure Insert_One_Variable
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 );
  procedure Insert_One_Variable
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 );

  -- DESCRIPTION :
  --   Inserts one extra variable to each non null entry in X,
  --   at position k, shifting variables at position k and higher
  --   one position further in each exponent vector.
  --   The exponent of the extra variable equals zero.

  -- REQUIRED : k <= number of variables in each non null entry of X.

  procedure Remove_One_Variable
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                k : in integer32 );
  procedure Remove_One_Variable
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 );
  procedure Remove_One_Variable
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 );

  -- DESCRIPTION :
  --   Removes the variable with index k from all polynomials in x.
  --   This is the inverse of Insert_One_Variable(x,k).

  procedure Insert_Variable_Pivot
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                i,j,k : in integer32 );
  procedure Insert_Variable_Pivot
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                i,j,k : in integer32 );
  procedure Insert_Variable_Pivot
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                i,j,k : in integer32 );

  -- DESCRIPTION :
  --   Inserts the variable at position k in x(i,j),
  --   with exponent equal to one.

  -- REQUIRED : Insert_One_Variable(k,x) has been executed.

  procedure Recondition
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                dim,s : in integer32 );
  procedure Recondition
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                dim,s : in integer32 );
  procedure Recondition
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
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

  function Random_Linear_Equation
              ( x : Standard_Complex_Poly_Matrices.Matrix;
                s : integer32 )
              return Standard_Complex_Polynomials.Poly;
  function Random_Linear_Equation
              ( x : DoblDobl_Complex_Poly_Matrices.Matrix;
                s : integer32 )
              return DoblDobl_Complex_Polynomials.Poly;
  function Random_Linear_Equation
              ( x : QuadDobl_Complex_Poly_Matrices.Matrix;
                s : integer32 )
              return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a random linear polynomial with the polynomials in x,
  --   in columns s and s+1.  The polynomial on return is the sum of
  --   the polynomials of x in columns s and s+1, multiplied each
  --   with a random complex constant.

  procedure Set_Exponent_to_Zero
              ( p : in out Standard_Complex_Polynomials.Poly;
                k : in integer32 );
  procedure Set_Exponent_to_Zero
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                k : in integer32 );
  procedure Set_Exponent_to_Zero
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                k : in integer32 );

  -- DESCRIPTION :
  --   For all terms in p, sets the exponent of the k-th variable
  --   to zero.  This procedure is applied to remove the variable
  --   (typically the last one) which is the continuation parameter.

  -- REQUIRED : k is in the range of 1..Number_of_Unknowns(p).

  procedure Add_Random_Constant
              ( p : in out Standard_Complex_Polynomials.Poly );
  procedure Add_Random_Constant
              ( p : in out DoblDobl_Complex_Polynomials.Poly );
  procedure Add_Random_Constant
              ( p : in out QuadDobl_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Adds a random complex constant to the polynomial p.

  -- REQUIRED : p /= Null_Poly.

  function Recondition_Target_Equation
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               s,t : integer32 )
             return Standard_Complex_Polynomials.Poly;
  function Recondition_Target_Equation
             ( x : DoblDobl_Complex_Poly_Matrices.Matrix;
               s,t : integer32 )
             return DoblDobl_Complex_Polynomials.Poly;
  function Recondition_Target_Equation
             ( x : QuadDobl_Complex_Poly_Matrices.Matrix;
               s,t : integer32 )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a random linear equation to recondition the columns
  --   s and s+1 in the matrix.  The variable t is the index of
  --   the variable used as the continuation parameter.

  function Recondition_Start_Equation
             ( n,k : integer32 )
             return Standard_Complex_Polynomials.Poly;
  function Recondition_Start_Equation
             ( n,k : integer32 )
             return DoblDobl_Complex_Polynomials.Poly;
  function Recondition_Start_Equation
             ( n,k : integer32 )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the start equation in the reconditioned swap homotopy,
  --   of the form x(k) - 1, where the polynomial on return is
  --   a polynomial in n variables.

  function Recondition_Equation
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               s,t,k : integer32 )
             return Standard_Complex_Polynomials.Poly;
  function Recondition_Equation
             ( x : DoblDobl_Complex_Poly_Matrices.Matrix;
               s,t,k : integer32 )
             return DoblDobl_Complex_Polynomials.Poly;
  function Recondition_Equation
             ( x : QuadDobl_Complex_Poly_Matrices.Matrix;
               s,t,k : integer32 )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the added linear equation in the reconditioned homotopy.

  -- ON ENTRY :
  --   x       matrix of indeterminates with the scaling variable;
  --   s       the column s in the swap homotopy;
  --   t       the index of the continuation parameter,
  --           which is the last variable in each polynomial;
  --   k       index of the variable x(r+1,s+1).

  function Recondition_Solution_Vector
             ( x : Standard_Complex_Vectors.Vector; k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector;
  function Recondition_Solution_Vector
             ( x : DoblDobl_Complex_Vectors.Vector; k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Vectors.Vector;
  function Recondition_Solution_Vector
             ( x : QuadDobl_Complex_Vectors.Vector; k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the scaled solution vector where every entry in x,
  --   in columns s and s+1 of the localization pattern locmap,
  --   has been divided by x(k) and 1/x(k) is the last added entry
  --   in the vector on return.
  --   Variables in column s which are multiplied by the pivot k in xp
  --   are not divided.

  -- REQUIRED : k in x'range and x(k) /= 0.

  function Recondition_Solution
             ( sol : Standard_Complex_Solutions.Solution;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Solutions.Solution;
  function Recondition_Solution
             ( sol : DoblDobl_Complex_Solutions.Solution;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Solutions.Solution;
  function Recondition_Solution
             ( sol : QuadDobl_Complex_Solutions.Solution;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns a solution with the same data as sol, but with its
  --   solution vector scaled with respect to its k-th coordinate.

  function Recondition_Solutions
             ( sols : Standard_Complex_Solutions.Solution_List;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Solutions.Solution_List;
  function Recondition_Solutions
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Solutions.Solution_List;
  function Recondition_Solutions
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   The list on return has all solution vectors scaled with
  --   respect to the k-th coordinate of the vector.

  function Rescale_Solution_Vector
             ( x : Standard_Complex_Vectors.Vector; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return Standard_Complex_Vectors.Vector;
  function Rescale_Solution_Vector
             ( x : DoblDobl_Complex_Vectors.Vector; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return DoblDobl_Complex_Vectors.Vector;
  function Rescale_Solution_Vector
             ( x : QuadDobl_Complex_Vectors.Vector; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   The vector on return has range x'first..x'last-1 and the 
  --   entries of x divided by the last coordinate in x,
  --   but only those entries of x which correspond to columns s and s+1
  --   with respect to the localization pattern in locmap.

  -- REQUIRED : x(x'last) /= 0.

  function Rescale_Solution
             ( sol : Standard_Complex_Solutions.Solution; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return Standard_Complex_Solutions.Solution;
  function Rescale_Solution
             ( sol : DoblDobl_Complex_Solutions.Solution; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return DoblDobl_Complex_Solutions.Solution;
  function Rescale_Solution
             ( sol : QuadDobl_Complex_Solutions.Solution; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns the solution sol with a rescaled solution vector,
  --   with respect to the swap column s in the localization pattern.

  -- REQUIRED : sol.v(sol.v'last) /= 0.

  function Rescale_Solutions
             ( sols : Standard_Complex_Solutions.Solution_List;
               s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return Standard_Complex_Solutions.Solution_List;
  function Rescale_Solutions
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return DoblDobl_Complex_Solutions.Solution_List;
  function Rescale_Solutions
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the list of rescaled solutions, rescaled sols,
  --   with respect to the swap column s in the localization pattern.

end Recondition_Swap_Homotopies;
