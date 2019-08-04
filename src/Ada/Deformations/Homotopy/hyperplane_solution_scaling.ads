with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;

package Hyperplane_Solution_Scaling is

-- DESCRIPTION :
--   Scales a solution vector and the corresponding last linear equation.

  procedure Sub ( p : in out Standard_Complex_Polynomials.Poly;
                  c : in Standard_Complex_Numbers.Complex_Number );
  procedure Sub ( p : in out DoblDobl_Complex_Polynomials.Poly;
                  c : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Sub ( p : in out QuadDobl_Complex_Polynomials.Poly;
                  c : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Subtracts the number c from the polynomial p.

  procedure Scale ( v : in out Standard_Complex_Vectors.Vector );
  procedure Scale ( v : in out DoblDobl_Complex_Vectors.Vector );
  procedure Scale ( v : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides every component in v by the largest element in v.

  procedure Adjust ( c : in Standard_Complex_Vectors.Link_to_Vector;
                     v : in Standard_Complex_Vectors.Vector );
  procedure Adjust ( c : in DoblDobl_Complex_Vectors.Link_to_Vector;
                     v : in DoblDobl_Complex_Vectors.Vector );
  procedure Adjust ( c : in QuadDobl_Complex_Vectors.Link_to_Vector;
                     v : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Adjusts the last coefficients of c so that the c*v equals zero,
  --   where * is the regular inner product (not Hermitian).

end Hyperplane_Solution_Scaling;
