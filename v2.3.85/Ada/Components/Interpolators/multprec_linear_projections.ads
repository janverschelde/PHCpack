with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;

package Multprec_Linear_Projections is

-- DESCRIPTION :
--   Provides a projection operator by evaluation in a general hyperplane
--   followed by truncation down to an n-dimensional vector.

  function Evaluate ( hyp : Multprec_Complex_Vectors.Vector;
                      point : Standard_Complex_Vectors.Vector )
                    return Complex_Number;

  function Evaluate ( hyp,point : Multprec_Complex_Vectors.Vector )
                    return Complex_Number;

  -- DESCRIPTION :
  --   Returns the evaluation of the point into the equation for
  --   the hyperplane:
  --     hyp(0) + hyp(1)*point(1) + .. + hyp(n)*point(n),
  --   where n = point'last.

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      point : Standard_Complex_Vectors.Vector;
                      n : integer32 ) return Multprec_Complex_Vectors.Vector;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      point : Multprec_Complex_Vectors.Vector;
                      n : integer32 ) return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the evaluation of the point in each of the hyperplanes
  --   given in hyps.  The vector on return has the range hyps'first..n.

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      points : Standard_Complex_VecVecs.VecVec;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec;
  function Evaluate ( hyps,points : Multprec_Complex_VecVecs.VecVec;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Applies the evaluation to every point in points.

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      arvv : Standard_Complex_VecVecs.Array_of_VecVecs;
                      n : integer32 )
                   return Multprec_Complex_VecVecs.Array_of_VecVecs;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      arvv : Multprec_Complex_VecVecs.Array_of_VecVecs;
                      n : integer32 )
                   return Multprec_Complex_VecVecs.Array_of_VecVecs;

  -- DESCRIPTION :
  --   Projects the vectors on an n-dimensional space by evaluation
  --   in the hyperplanes and truncation to n-dimensional vectors.

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      sols : Standard_Complex_Solutions.Solution_List;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      sols : Multprec_Complex_Solutions.Solution_List;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec;
  
  -- DESCRIPTION :
  --   Projects the solution vectors on an n-dimensional space.

end Multprec_Linear_Projections;
