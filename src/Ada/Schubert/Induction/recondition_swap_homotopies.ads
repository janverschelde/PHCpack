with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Solutions;

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

  function Random_Linear_Equation
              ( x : Standard_Complex_Poly_Matrices.Matrix;
                s : integer32 )
              return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a random linear polynomial with the polynomials in x,
  --   in columns s and s+1.  The polynomial on return is the sum of
  --   the polynomials of x in columns s and s+1, multiplied each
  --   with a random complex constant.

  procedure Set_Exponent_to_Zero
              ( p : in out Standard_Complex_Polynomials.Poly;
                k : in integer32 );

  -- DESCRIPTION :
  --   For all terms in p, sets the exponent of the k-th variable
  --   to zero.  This procedure is applied to remove the variable
  --   (typically the last one) which is the continuation parameter.

  -- REQUIRED : k is in the range of 1..Number_of_Unknowns(p).

  procedure Add_Random_Constant
              ( p : in out Standard_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Adds a random complex constant to the polynomial p.

  -- REQUIRED : p /= Null_Poly.

  function Recondition_Target_Equation
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               s,t : integer32 )
             return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a random linear equation to recondition the columns
  --   s and s+1 in the matrix.  The variable t is the index of
  --   the variable used as the continuation parameter.

  function Recondition_Start_Equation
             ( n,k : integer32 )
             return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the start equation in the reconditioned swap homotopy,
  --   of the form x(k) - 1, where the polynomial on return is
  --   a polynomial in n variables.

  function Recondition_Equation
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               s,t,k : integer32 )
             return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the added linear equation in the reconditioned homotopy.

  -- ON ENTRY :
  --   x       matrix of indeterminates with the scaling variable;
  --   s       the column s in the swap homotopy;
  --   t       the index of the continuation parameter,
  --           which is the last variable in each polynomial;
  --   k       index of the variable x(r+1,s+1).

  function Recondition_Solution_Vector
             ( x : Standard_Complex_Vectors.Vector; k : integer32 )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the scaled solution vector where every entry in x
  --   has been divided by x(k) and 1/x(k) is the last added entry
  --   in the vector on return.

  -- REQUIRED : k in x'range and x(k) /= 0.

  function Recondition_Solution
             ( s : Standard_Complex_Solutions.Solution; k : integer32 )
             return Standard_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns a solution with the same data as s, but with its
  --   solution vector scaled with respect to its k-th coordinate.

  function Recondition_Solutions
             ( sols : Standard_Complex_Solutions.Solution_List;
               k : integer32 )
             return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   The list on return has all solution vectors scaled with
  --   respect to the k-th coordinate of the vector.

end Recondition_Swap_Homotopies;
