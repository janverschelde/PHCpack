with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;

package Standard_Central_Projections is

-- DESCRIPTION :
--   Provides projection operators by intersecting a line with
--   a general hyperplane for standard complex numbers.

-- PROJECTION OF ONE POINT :

  function Intersect ( hyp,base,point : Vector; evabas : Complex_Number; 
                       dim : integer32 ) return Vector;
  function Intersect ( hyp,base,point : Vector; dim : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the intersection of the points "base" and "point" with
  --   the hyperplane:  hyp(0) + hyp(1)*x(1) + .. + hyp(n)*x(n) = 0,
  --   where n = point'last.
  --   The hyperplane hyp should be general enough so that the line
  --   between the two points has an intersection point.

  -- ON ENTRY :
  --   hyp      coefficients of hyperplane:
  --              hyp(0) + hyp(1)*x(1) + .. + hyp(n)*x(n) = 0;
  --   base     base point in the projection;
  --   point    point that is subject to projection;
  --   dim      length of the returning vector;
  --   evabas   evaluation of base point, is hyp(base'range)*base.

  procedure Intersect_Base_Points ( hyp : in VecVec; base : in out VecVec );

  -- DESCRIPTION :
  --   Projects the base points successively onto the hyperplanes:
  --   base(1) is used to project base(2..base'last) onto hyp(1),
  --   base(2) is used to project base(3..base'last) onto hyp(2), etc...
  --   The points in base(2..base'last) satisfy hyp(1),
  --   the points in base(3..base'last) satisfy hyp(2), etc...

  function Intersect ( hyp,base : VecVec; point : Vector; dim : integer32 )
                     return Vector;

  -- DESCRIPTION :
  --   Returns the result of successive projections of the point with
  --   several base points, using the intersection construction above.

  -- REQUIRED :
  --   The base points have been filtered through "Intersect_Base_Points".

-- PROJECTION OF A SEQUENCE OF POINTS :

  function Intersect ( hyp,base : Vector; evabas : Complex_Number;
                       points : VecVec; dim : integer32 ) return VecVec;
  function Intersect ( hyp,base : Vector; points : VecVec; dim : integer32 )
                     return VecVec;

  -- DESCRIPTION :
  --   Applies the Intersect projection operator to all points.

end Standard_Central_Projections;
