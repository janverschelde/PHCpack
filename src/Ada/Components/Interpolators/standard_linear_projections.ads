with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Linear_Projections is

-- DESCRIPTION :
--   Provides a projection operator by evaluation in a general hyperplane
--   followed by truncation down to an n-dimensional vector.

  function Evaluate ( hyp,point : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the evaluation of the point into the equation for
  --   the hyperplane:
  --     hyp(0) + hyp(1)*point(1) + .. + hyp(n)*point(n),
  --   where n = point'last.

  function Evaluate ( hyps : VecVec; point : Vector; n : integer32 )
                    return Vector;

  -- DESCRIPTION :
  --   Returns the evaluation of the point in each of the hyperplanes
  --   given in hyps.  The vector on return has the range hyps'first..n.

  function Evaluate ( hyps,pts : VecVec; n : integer32 ) return VecVec;
  function Evaluate ( hyps : VecVec; arvv : Array_of_VecVecs; n : integer32 )
                   return Array_of_VecVecs;

  -- DESCRIPTION :
  --   Projects the vectors on an n-dimensional space by evaluation
  --   in the hyperplanes and truncation to n-dimensional vectors.

  function Evaluate ( hyps : VecVec; sols : Solution_List; n : integer32 )
                    return VecVec;
  
  -- DESCRIPTION :
  --   Projects the solution vectors on an n-dimensional space.

end Standard_Linear_Projections;
