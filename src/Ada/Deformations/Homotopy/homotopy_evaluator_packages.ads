with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Homotopy_Evaluator_Packages is

-- DESCRIPTION :
--   Provides routines to create package to evaluate the homotopy
--     h(x,t) = a*(1-t)^k*q(x) + t^k*p(x) = 0, for t in [0,1]. 
--   If the package name is not provided, then it will be read.

  procedure Create ( packname : in String; p,q : in Poly_Sys );
  procedure Create ( p,q : in Poly_Sys );

end Homotopy_Evaluator_Packages;
