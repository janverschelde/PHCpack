with Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_Polynomials;
with Standard_Floating_Poly_Functions;
with Standard_Floating_Poly_Systems;
with Generic_Poly_System_Functions;

package Standard_Floating_Poly_SysFun is
  new Generic_Poly_System_Functions(Standard_Floating_Ring,
                                    Standard_Floating_Vectors,
                                    Standard_Floating_VecVecs,
                                    Standard_Floating_Polynomials,
                                    Standard_Floating_Poly_Functions,
                                    Standard_Floating_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials 
--   with standard real floating-point coefficients.
