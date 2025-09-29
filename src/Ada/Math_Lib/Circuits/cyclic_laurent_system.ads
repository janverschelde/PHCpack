with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Cyclic_Laurent_System is

-- DESCRIPTION :
--   Reformulates the cyclic n-root problem as a Laurent system,
--   following the paper of Uffe Haagerup on "Cyclic p-roots of
--   prime length p and related complex Hadamard matrices.

  function Monomial_Support
             ( i,j,n : natural32 ) return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the support vector of the j-th monomial
  --   in the i-th equation, of cyclic n-roots, where i starts at 0,
  --   as the first equation is x0 = 1.
  --   The x0 = 1 is omitted and the reformulated cyclic n-roots system
  --   has n-1 equations in n-1 variables.

  function Cyclic_Monomial ( i,j,n : natural32 ) return Term;

  -- DESCRIPTION :
  --   Returns the j-th term in the i-th polynomial of the
  --   reformulated cyclic n-roots problem.

  function Polynomial_Support ( i,n : natural32 ) return List;

  -- DESCRIPTION :
  --   Returns the support of the i-th polynomial 
  --   of the cyclic n-roots system.

  function Cyclic_Polynomial ( i,n : natural32 ) return Poly;

  -- DESCRIPTION :
  --   Returns the i-th polynomial of the reformulated cyclic n-roots system.

  function System_Support ( n : natural32 ) return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the supports of the cyclic n-roots system.

  function Cyclic_System ( n : natural32 ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the reformulated cyclic n-roots system
  --   as a system of n-1 polynomials in n-1 variables.

end Cyclic_Laurent_System;
