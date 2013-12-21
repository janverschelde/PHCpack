with Abstract_Ring;
with Abstract_Ring.Field;
with Generic_Vectors;
with Generic_Matrices;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Field is new Ring.Field(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);

package Generic_Norms_Equals is

-- DESCRIPTION :
--   Provides norms of vectors and decision routines for equalities.

  use Ring,Field,Vectors,Matrices;

  function Max_Norm ( v : Vector ) return number;

  -- DESCRIPTION :
  --   Returns the absolute value of the greatest element in v.

  function Sum_Norm ( v : Vector ) return number;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the elements in v.

  function Max_Norm ( m : Matrix ) return number;

  -- DESCRIPTION :
  --   Returns the maximal element in the matrix in absolute value.

  function Equal ( x,y,tol : number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if abs(x-y) < tol, otherwise false is returned.

  function Equal ( x,y : Vector; tol : number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Equal(x(i),y(i),tol), for i in x'range=y'range,
  --   otherwise false is returned.

  -- REQUIRED : x'range = y'range.

end Generic_Norms_Equals;
