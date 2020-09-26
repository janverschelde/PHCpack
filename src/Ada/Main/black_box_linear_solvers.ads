with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;

package Black_Box_Linear_Solvers is

-- DESCRIPTION :
--   Solvers for systems that are linear and
--   that have as many equations as unknowns.

  procedure Solve ( infilename,outfilename : in string;
                    p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                    n : in natural32; append_sols : in boolean;
                    fail : out boolean; verbose : in integer32 := 0 );
  procedure Solve ( infilename,outfilename : in string;
                    p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    n : in natural32; append_sols : in boolean;
                    fail : out boolean; verbose : in integer32 := 0 );
  procedure Solve ( infilename,outfilename : in string;
                    p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    n : in natural32; append_sols : in boolean;
                    fail : out boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Parses the system to see if it is square and linear.
  --   If square and linear, then it is solved.

  -- ON ENTRY :
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   p              a polynomial system;
  --   n              number of variables in the polynomials of p;
  --   append_sols    true if solutions need to be appended to input file;
  --   verbose        the verbose level.

  -- ON RETURN :
  --   fail           true if system p is nonlinear.

end Black_Box_Linear_Solvers;
