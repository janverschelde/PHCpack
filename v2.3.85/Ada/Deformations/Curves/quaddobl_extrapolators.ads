with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;

package QuadDobl_Extrapolators is

-- DESCRIPTION :
--   Application of divided differences for extrapolation
--   with double double precision arithmetic.

-- SCALAR VERSIONS :

  function Extrapolate ( t,t0,t1,x0,x1 : Complex_Number )
                       return Complex_Number;

  -- DESCRIPTION :
  --   Returns x(t), the value of the line through (t0,x0) and (t1,x1).

  -- REQUIRED : t0 /= t1.

  function Extrapolate ( t,t0,t1,t2,x0,x1,x2 : Complex_Number )
                       return Complex_Number;

  -- DESCRIPTION :
  --   Returns x(t), the value of the quadric through the points
  --   (t0,x0), (t1,x1), and (t2,x2).

  -- REQUIRED : the coordinates t0, t1, and t2 are mutually distinct.

  function Extrapolate ( t,t0,t1,t2,t3,x0,x1,x2,x3 : Complex_Number )
                       return Complex_Number;

  -- DESCRIPTION :
  --   Returns x(t), the value of the cubic through the points
  --   (t0,x0), (t1,x1), (t2,x2), and (t3,x3).

  -- REQUIRED : the coordinates t0, t1, t2, and t3 are mutually distinct.

  function Extrapolate ( t,t0,t1,t2,t3,t4,x0,x1,x2,x3,x4 : Complex_Number )
                       return Complex_Number;

  -- DESCRIPTION :
  --   Returns x(t), the value of the quartic through the points
  --   (t0,x0), (t1,x1), (t2,x2), (t3,x3), and (t4,x4).

  -- REQUIRED : the coordinates t0, t1, t2, t3, and t4 are mutually distinct.

  function Extrapolate
             ( t,t0,t1,t2,t3,t4,t5,x0,x1,x2,x3,x4,x5 : Complex_Number )
             return Complex_Number;

  -- DESCRIPTION :
  --   Returns x(t), the value of the quintic through the points
  --   (t0,x0), (t1,x1), (t2,x2), (t3,x3), (t4,x4), and (t5,x5).

  -- REQUIRED : the coordinates t0, t1, t2, t3, t4, and t5 
  --   are mutually distinct.

-- VECTOR VERSIONS :

  function Extrapolate ( t,t0,t1 : Complex_Number; x0,x1 : Vector )
                       return Vector;

  -- DESCRIPTION :
  --   Returns x(t), the value of the line through (t0,x0) and (t1,x1).

  -- REQUIRED : t0 /= t1.

  function Extrapolate ( t,t0,t1,t2 : Complex_Number; x0,x1,x2 : Vector )
                       return Vector;

  -- DESCRIPTION :
  --   Returns x(t), the value of the quadric through the points
  --   (t0,x0), (t1,x1), and (t2,x2).

  -- REQUIRED : the coordinates t0, t1, and t2 are mutually distinct.

  function Extrapolate ( t,t0,t1,t2,t3 : Complex_Number;
                         x0,x1,x2,x3 : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns x(t), the value of the cubic through the points
  --   (t0,x0), (t1,x1), (t2,x2), and (t3,x3).

  -- REQUIRED : the coordinates t0, t1, t2, and t3 are mutually distinct.

  function Extrapolate ( t,t0,t1,t2,t3,t4 : Complex_Number;
                         x0,x1,x2,x3,x4 : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns x(t), the value of the quartic through the points
  --   (t0,x0), (t1,x1), (t2,x2), (t3,x3), and (t4,x4).

  -- REQUIRED : the coordinates t0, t1, t2, t3, and t4 are mutually distinct.

  function Extrapolate ( t,t0,t1,t2,t3,t4,t5 : Complex_Number;
                         x0,x1,x2,x3,x4,x5 : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns x(t), the value of the quintic through the points
  --   (t0,x0), (t1,x1), (t2,x2), (t3,x3), (t4,x4), and (t5,x5).

  -- REQUIRED : the coordinates t0, t1, t2, t3, t4, and t5 
  --   are mutually distinct.

end QuadDobl_Extrapolators;
