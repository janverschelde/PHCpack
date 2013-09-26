with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Laur_Polys;        use Standard_Complex_Laur_Polys;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Weighted_Projective_Transfos is

-- DESCRIPTION :
--   Weighted projective embeddings are determined by an integral weight 
--   vector q : replace x^a by x^a*z^<a,q>.  The homogeneous coordinates
---  are then (l^(-1) x0, l^q1 x1 , l^q2 x2,.., l^qn x_n), for l /= 0.
--   Set l := x0 to move from projective to affine space.
--   Usual projective space corresponds to q = (-1,-1,..,-1).

  function Projective_Transformation ( p : Poly; q : Vector ) return Poly;
  function Projective_Transformation
             ( p : Laur_Sys; q : Vector ) return Laur_Sys;

  procedure Projective_Transformation ( p : in out Poly; q : in Vector );
  procedure Projective_Transformation ( p : in out Laur_Sys; q : in Vector );

  -- DESCRIPTION :
  --   An unknown using q as weight is added according to the recipe above.

  function Affine_Transformation ( s : Solution; q : Vector ) return Solution;

  function Affine_Transformation
             ( sols : Solution_List; q : Vector ) return Solution_List;

  procedure Affine_Transformation
              ( sols : in out Solution_List; q : in Vector );

  -- DESCRIPTION :
  --   All components of the solution vector will be divided by the last
  --   component, which is afterwards cut off.

end Weighted_Projective_Transfos;
