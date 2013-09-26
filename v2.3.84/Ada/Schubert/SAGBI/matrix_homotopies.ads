with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;

package Matrix_Homotopies is

-- DESCRIPTION :
--   This package provides an abstraction for dealing with homotopies
--   of matrices.
--   In particular, we can store and evaluate n different couples of start and
--   target matrices, defining the matrix(t) := (1-t)*start + t*target.

-- CREATORS :

  procedure Init ( n : in natural32 );

  -- DESCRIPTION :
  --   Reserves space for n matrix homotopies.

  procedure Add ( start,target : in Matrix );

  -- DESCRIPTION :
  --   Adds a new couple to the matrix homotopies.

  -- REQUIRED : start'range(1) = target'range(1) = 1..n
  --        and start'range(2) = target'range(2) = 1..m,
  --        and not exceed initial capacity n.

  procedure Add_Start ( mapno : in natural32; start : in Matrix );
  procedure Add_Target ( mapno : in natural32; target : in Matrix );

  -- DESCRIPTION :
  --   Adds the start or target system for the indicated map.
  --   This can also be used to modify the homotopies.

  -- REQUIRED : mapno <= n, start and target should have same dimensions.

-- SELECTOR :

  function Empty ( mapno : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the matrix homotopy indicated by mapno is empty.

  function Cardinality return natural32;

  -- DESCRIPTION :
  --   Returns the number of matrix homotopies that have been added.

-- EVALUATOR :

  function Eval ( mapno : natural32; t : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   Evaluates the kth matrix homotopy at t, with k = mapno.

  -- REQUIRED : mapno <= n.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Destroys all matrices.

end Matrix_Homotopies;
