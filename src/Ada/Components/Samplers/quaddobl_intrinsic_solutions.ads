with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;           use QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_Intrinsic_Solutions is

-- DESCRIPTION :
--   A solution in a k-space has a short coordinate vector of length k.
--   Moving from the long representation in the ambient space to the
--   intrinsic shorter k-vector is done by projection.  Expanding the
--   k-vector with the parametric representation of the linear space,
--   we obtain again the longer representation in the ambient space.
--   Calculations happen in quad double complex arithmetic.

-- NOTICE :
--   The projection operator forgets the values for the slack variables
--   (at the end of each coordinate vector), so that expanding after a
--   projecting a solution returns the original solution without slacks.

  function Project ( s : Solution; b : Vector; v : VecVec ) return Solution;
  function Project ( s : Solution; p : Matrix ) return Solution;

  -- DESCRIPTION :
  --   Given a linear space with offset vector b and directions in v,
  --   or by the columns in a matrix p, the solution on return has 
  --   coordinates projected in this linear space. 
  -- REQUIRED : The direction vectors form an orthonormal basis.

  function Project ( sols : Solution_List; b : Vector; v : VecVec )
                   return Solution_List;
  function Project ( sols : Solution_List; p : Matrix ) return Solution_List;

  -- DESCRIPTION :
  --   Projects every solution in the list onto the linear space.
  -- REQUIRED : The direction vectors are othonormal to each other.

  function Expand ( s : Solution; b : Vector; v : VecVec ) return Solution;
  function Expand ( s : Solution; p : Matrix ) return Solution;

  -- DESCRIPTION :
  --   Given coordinates of a solution in a linear space defined by b and v,
  --   or defined by p, Expand returns coordinates in the ambient space.

  function Expand ( sols : Solution_List; b : Vector; v : VecVec )
                  return Solution_List;
  function Expand ( sols : Solution_List; p : Matrix ) return Solution_List;

  -- DESCRIPTION :
  --   Expands every solution to the representation in the ambient space.

  function Transform ( s : Solution; p1,p2 : Matrix ) return Solution;

  -- DESCRIPTION :
  --   Transforms the coordinate system for the solution s,
  --   given on entry in intrinsic coordinates in p1, 
  --   the solution on return has intrinsic coordinates in p2.

  -- REQUIRED : the directions in p2 are orthonormal.

  function Transform ( sols : Solution_List; p1,p2 : Matrix )
                     return Solution_List;
  procedure Transform ( sols : in out Solution_List; p1,p2 : in Matrix );

  -- DESCRIPTION :
  --   Applies the Transform to every solution in the list sols.

end QuadDobl_Intrinsic_Solutions;
