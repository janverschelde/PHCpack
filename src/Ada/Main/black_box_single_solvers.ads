with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;

package Black_Box_Single_Solvers is

-- DESCRIPTION :
--   A solver for one single polynomial applies a univariate root finder
--   if there is only one variable, or factors the polynomial if there
--   are many variables.

  procedure Solve ( infilename,outfilename : in string;
                    p : in Standard_Complex_Polynomials.Poly;
                    append_sols : in boolean; verbose : in integer32 := 0 );
  procedure Solve ( infilename,outfilename : in string;
                    p : in DoblDobl_Complex_Polynomials.Poly;
                    append_sols : in boolean; verbose : in integer32 := 0 );
  procedure Solve ( infilename,outfilename : in string;
                    p : in QuadDobl_Complex_Polynomials.Poly;
                    append_sols : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Solves one single polynomial p,
  --   distinguishing between 1 and several variables,
  --   with double, double double, or quad double arithmetic.

  -- ON ENTRY :
  --   infilename   name of the input file;
  --   outfilename  name of the output file;
  --   p            a polynomial in one variable;
  --   append_sols  whether the solutions should be added to the input file;
  --   verbose      the verbose level.

end Black_Box_Single_Solvers;
