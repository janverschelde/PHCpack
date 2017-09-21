with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems; 
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;

package Reduction_of_Polynomial_Systems is

-- DESCRIPTION :
--   Linear and nonlinear reduction to reduce the total degree.

  procedure Reduce ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean );
  procedure Reduce ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean );
  procedure Reduce ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean );

  -- DESCRIPTION :
  --   This procedure tries to lower the total degree of p by means
  --   of linear reduction, in standard double, double double, or
  --   quad double precision.

  -- ON ENTRY : 
  --   p             a polynomial system.

  -- ON RETURN :
  --   p             a polynomial system with a possible lower total degree;
  --   diagonal      true if all leading terms in p are different;
  --   inconsistent  is true if the reduced system has equations `4=0';
  --   infinite      is true if some equations of the original system
  --                 disappeared during the reduction process.

  procedure Sparse_Reduce ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                            inconsistent,infinite : in out boolean );
  procedure Sparse_Reduce ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                            inconsistent,infinite : in out boolean );
  procedure Sparse_Reduce ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                            inconsistent,infinite : in out boolean );

  -- DESCRIPTION :
  --   This procedure makes the coefficient matrix of p as sparse as
  --   possible, in double, double double, or quad double precision.

  procedure Reduce ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     res : in out Standard_Complex_Poly_Systems.Poly_Sys;
                     cnt_eq : in out natural32; max_eq : in natural32;
                     cnt_sp : in out natural32; max_sp : in natural32;
                     cnt_rp : in out natural32; max_rp : in natural32 );

  -- DESCRIPTION : 
  --   This procedure tries to lower the total degree of the system p
  --   by means of nonlinear reduction. 

  -- REQUIRED : the counters must equal 0, on entry.

  -- ON ENTRY :
  --   p             a polynomial system;
  --   cnt_eq        counts the number of equal degree substitutions;
  --   max_eq        limit on the number of equal degree substitutions;
  --   cnt_sp        counts the number of S-polynomial computations;
  --   max_sp        limit on the number of S-polynomial computations.
  --   cnt_rp        counts the number of R-polynomial computations;
  --   max_rp        limit on the number of R-polynomial computations.

  -- ON RETURN :
  --   res           the reduced system;
  --   cnt_eq        the number of equal degree substitutions;
  --   cnt_sp        the number of computed S-polynomials;
  --   cnt_rp        the number of computed R-polynomials.

  procedure Sparse_Reduce
                   ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     res : in out Standard_Complex_Poly_Systems.Poly_Sys;
                     cnt_eq : in out natural32;
                     max_eq : in natural32 );

  -- DESCRIPTION :
  --   the polynomial system is reduced by computing S-polynomials.
  --   After each replacement, the coefficient matrix is made as sparse
  --   as possible.

end Reduction_of_Polynomial_Systems;
