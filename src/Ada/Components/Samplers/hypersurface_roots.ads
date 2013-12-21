with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Complex_Vectors;         use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;     use Standard_Complex_Polynomials;

package Hypersurface_Roots is

-- DESCRIPTION :
--   The problem of sampling points on a hypersurface is reduced to
--   that of finding roots of a univariate polynomial.

  function Substitute ( p : Poly; v : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the coefficient of the polynomial in one variable,
  --   obtained as p(t*v), substituting x(i) by v(i)*t.

  procedure Scale ( v : in out Vector );

  -- DESCRIPTION :
  --   Divides the vector by its max norm. 

  procedure Affine_Path_Tracker
              ( d : in natural32; p : in Poly; v0,v1 : in Vector; 
                s : in out Complex_Number; nbsteps : out natural32 );
  procedure Projective_Path_Tracker
              ( d : in natural32; p : in Poly; v0,v1 : in Vector; 
                s0,s1 : in out Complex_Number; nbsteps : out natural32 );

  -- DESCRIPTION :
  --   Tracks one path in affine or projective coordinates from v0 to v1,
  --   using the homotopy v0*t + v1*(1-t), for t going from 1 to 0.
  --   The intermediate polynomials are defined by p(v*s), where p
  --   is a polynomial of degree d.

  procedure Affine_Track_Moving_Line
              ( p : in Poly; v0,v1 : in Vector; s : in out Vector );
  procedure Projective_Track_Moving_Line
              ( p : in Poly; v0,v1 : in Vector; s0,s1 : in out Vector );

  -- DESCRIPTION :
  --   Tracks the solutions in s for p(v0*s) = 0 to p(v1*s) = 0,
  --   as the line moves from v0 to v1.

end Hypersurface_Roots;
