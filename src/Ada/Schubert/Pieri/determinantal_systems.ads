with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems; 
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Brackets;                           use Brackets;

package Determinantal_Systems is

-- DESCRIPTION :
--   This package constitutes the bridge between linear subspace intersections
--   expressed by determinantal equations and polynomial systems.

-- LOCALIZATION MAPS :

  function Standard_Coordinate_Frame
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               plane : Standard_Complex_Matrices.Matrix )
             return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a coordinate frame that corresponds to the polynomial matrix x.
  --   To determine the position of the one, the first nonzero element in the
  --   columns of the plane is taken.

  function Maximal_Coordinate_Frame
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               plane : Standard_Complex_Matrices.Matrix )
             return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a coordinate frame that corresponds to the polynomial matrix x.
  --   To determine the position of the one, the maximal element in every
  --   column is taken.

  function Localize ( locmap : Standard_Natural_Matrices.Matrix;
                      p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Applies the localization map to the polynomial system p,
  --   assuming a rowwise lexicographical order on the variables.

-- CONSTRUCT POLYNOMIAL SYSTEMS :

  procedure Concat ( L : in out Link_to_Poly_Sys; p : Poly_Sys );

  -- DESCRIPTION :
  --   Concatenates the nonzero polynomials in p to those in L.
  --   There is sharing.

  function Polynomial_Equations
              ( L : Standard_Complex_Matrices.Matrix;
                x : Standard_Complex_Poly_Matrices.Matrix ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial equations that express the intersection of
  --   the subspace spanned by the columns of l with the indeterminates
  --   in the variable matrix x.

  function Polynomial_Equations
              ( L : Standard_Complex_VecMats.VecMat;
                x : Standard_Complex_Poly_Matrices.Matrix ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system generated from expanding the 
  --   determinants that represent the intersections of x with the
  --   planes in L, with equations of the form det(L|x) = 0.

-- EVALUATORS AND DIFFERENTIATORS :

  function Eval ( L,x : Standard_Complex_Matrices.Matrix )
                return Complex_Number;

  -- DESCRIPTION :
  --   Returns the result of the evaluation of the determinant which
  --   has in its columns the columns of L and x.

  -- REQUIRED : x'length(1) = x'length(2) + L'length(2) = L'length(1).

  function Eval ( L,x : Standard_Complex_Matrices.Matrix ) return Vector;

  -- DESCRIPTION :
  --   Returns the vector of all maximal minors of the matrix that has
  --   in its columns the columns of l and x.

  -- REQUIRED : x'length(1) = L'length(1) >= x'length(2) + L'length(2).

  function Diff ( L,x : Standard_Complex_Matrices.Matrix; i : integer32 )
                return Complex_Number;
  function Diff ( L,x : Standard_Complex_Matrices.Matrix;
                  locmap : Standard_Natural_Matrices.Matrix; i : integer32 )
                return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the derivative of the i-th variable in the
  --   vector representation of the matrix x of unknowns, with coefficients
  --   obtained from the determinantal equation det(l|x) = 0.
  --   The same requirements as in Eval on the dimensions hold.
  --   The i-th variable is w.r.t. the localization map.

  function Eval ( L : Standard_Complex_VecMats.VecMat;
                  x : Standard_Complex_Matrices.Matrix ) return Vector;

  -- DESCRIPTION :
  --   Returns the vector that has in its components the determinant
  --   of l(i) with x, for i in l'range.

  function Diff ( L : Standard_Complex_VecMats.VecMat;
                  x : Standard_Complex_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix;
  function Diff ( L : Standard_Complex_VecMats.VecMat;
                  x : Standard_Complex_Matrices.Matrix; nvars : integer32;
                  locmap : Standard_Natural_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Jacobian matrix of the system det(l(i)|x) = 0.
  --   If the localization map is specificied, then nvars equals the
  --   number of free variables.

-- SOLUTIONS AND SYSTEMS FOR QUANTUM PIERI :

  function Solution_Plane
              ( top,bottom : Bracket;
                locmap : Standard_Natural_Matrices.Matrix;
                mat : Standard_Complex_Matrices.Matrix )
              return Solution;

  -- DESCRIPTION :
  --   Returns the representation of the solution plane as a vector.

  function Solution_Planes
              ( top,bottom : Bracket;
                locmap : Standard_Natural_Matrices.Matrix;
                vm : Standard_Complex_VecMats.VecMat ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the representation of the vector of planes as a solution list.

  function Create_Polynomial_System
              ( top,bottom : Bracket;
                locmap : Standard_Natural_Matrices.Matrix;
                xpm : Standard_Complex_Poly_Matrices.Matrix;
                svals : Standard_Complex_Vectors.Vector;
                planes : Standard_Complex_VecMats.VecMat ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system that collects the intersection
  --   conditions for meeting the given m-planes at the specified s-values.
  --   The system is localized according to the given localization map.

end Determinantal_Systems;
