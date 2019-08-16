with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
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

  function Max_Norm ( v : Standard_Complex_Vectors.Vector;
                      k : natural32;
                      z : Standard_Natural_Vectors.Vector ) 
                    return double_float;
  function Max_Norm ( v : DoblDobl_Complex_Vectors.Vector;
                      k : natural32;
                      z : Standard_Natural_Vectors.Vector ) 
                    return double_double;
  function Max_Norm ( v : QuadDobl_Complex_Vectors.Vector;
                      k : natural32;
                      z : Standard_Natural_Vectors.Vector ) 
                    return quad_double;

  -- DESCRIPTION :
  --   Returns the maximum of all coordinates of v
  --   in the k-th set of the partition z.

  procedure Scale ( v : in out Standard_Complex_Vectors.Vector;
                    m : in natural32;
                    z : in Standard_Natural_Vectors.Vector );
  procedure Scale ( v : in out DoblDobl_Complex_Vectors.Vector;
                    m : in natural32;
                    z : in Standard_Natural_Vectors.Vector );
  procedure Scale ( v : in out QuadDobl_Complex_Vectors.Vector;
                    m : in natural32;
                    z : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides every component in v by the largest element in its set,
  --   as defined by the m-homogeneous coordinates in z.

  -- ON ENTRY :
  --   v        a vector where the last m coordinates contain the
  --            values of the added variables in the m-homogenization;
  --   m        the number of sets in the partition of the variables;
  --   z        the index representation of the partition,
  --            z'last = n = the number of original affine variables,
  --            z(k) defines the set index i for which the (n+i)-th
  --            variable is the added one for the i-th set in z.

  -- ON RETURN :
  --   v        every coordinate is divided by the largest element
  --            in its set.

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
