with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;
with Symmetry_Group;                     use Symmetry_Group;

package Main_Verification is

-- DESCRIPTION :
--   Verification of solutions consists in running root refiners,
--   computing of multiplicity structure and winding numbers,
--   and the verification of polyhedral end games.

  procedure Display_Verification_Info;

  -- DESCRIPTION :
  --   Displays information about available verification methods on screen.

  procedure Refine_Roots
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                solsfile,invar,allperms,signsym : in boolean;
                v : in List_of_Permutations;
                epsxa,epsfa,tolsing : in double_float;
                maxit : in natural32; deflate : in out boolean;
                wout : in boolean;
                sols : in out Standard_Complex_Solutions.Solution_List;
                refsols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Refines the roots and computes generating solutions when required.

  -- ON ENTRY :
  --   file        for writing results on;
  --   p           the polynomial system under consideration;
  --   solsfile    whether refined solution have to go to separate file;
  --   invar       whether generating solutions have to be computed;
  --   allperms    whether invariant under all permutations;
  --   signsym     whether there is sign-symmetry;
  --   v           group representation, only needed when invar;
  --   sols        solutions that need to be refined.

  -- ON RETURN :
  --   sols        solutions after applying some Newton iteration;
  --   refsols     refined solutions, with the exception of failures and
  --               the non-generating solutions.

  procedure Refine_Roots
             ( file : in file_type;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               solsfile : in boolean;
               epsxa,epsfa,tolsing : in double_float; maxit : in natural32;
               deflate : in out boolean; wout : in boolean;
               sols : in out Standard_Complex_Solutions.Solution_List;
               refsols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION : 
  --   Root refinement without computing of generating solutions.

  procedure Refine_Roots
             ( file : in file_type;
               p : in Standard_Complex_Laur_Systems.Laur_Sys;
               solsfile : in boolean;
               epsxa,epsfa,tolsing : in double_float; maxit : in natural32;
               wout : in boolean;
               sols : in out Standard_Complex_Solutions.Solution_List;
               refsols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION : 
  --   Root refinement without computing of generating solutions.

-- VERIFICATION PROCEDURES :

  procedure Winding_Verification
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Verification by computing winding numbers by homotopy continuation.

  procedure Standard_Weeding_Verification
              ( infile : in out file_type; outfilename : in string;
                lp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sysonfile : in boolean );

  -- DESCRIPTION :
  --   Verification by refining the roots and weeding out the solution set.

  procedure Standard_Weeding_Verification
              ( infile : in out file_type; outfilename : in string;
                lp : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                sysonfile : in boolean );

  -- DESCRIPTION :
  --   Verification by refining the roots and weeding out the solution set.

  procedure Standard_Weeding_Verification
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Verification by refining the roots and weeding out the solution set.

  procedure Multprec_Residual_Evaluator
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluation of residuals using multi-precision arithmetic.

  procedure Call_Multprec_Root_Refiner
               ( file : in file_type; n,m : in natural32;
                 ls : in Link_to_Array_of_Strings;
                 sols : in out Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Calls the multiprecision root refiner.

  procedure Call_Varbprec_Root_Refiner
               ( file : in file_type;
                 ls : in Link_to_Array_of_Strings;
                 sols : in out Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Offers the user a menu to set the parameters and then calls
  --   the variable precision Newton's method on the list of solutions.

  procedure Multprec_Weeding_Verification
              ( infilename,outfilename : in string; varbprec : in boolean );

  -- DESCRIPTION :
  --   Newton's method using multi-precision arithmetic,
  --   or with variable precision if varbprec is true.

  procedure Polyhedral_End_Game_Verification;

  -- DESCRIPTION :
  --   Verification of the polyhedral end game.

  procedure Standard_Newton_with_Deflation
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads the polynomial system in standard double precision
  --   confirms the output to the file and calls the deflation method.
  --   The verbose level is in vrblvl.

  procedure DoblDobl_Newton_with_Deflation 
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads the polynomial system in double double precision
  --   confirms the output to the file and calls the deflation method.
  --   The verbose level is in vrblvl.

  procedure QuadDobl_Newton_with_Deflation
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads the polynomial system in quad double precision
  --   confirms the output to the file and calls the deflation method.
  --   The verbose level is in vrblvl.

  procedure Newton_with_Deflation
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and then calls the corresponding driver
  --   to apply Newton's method with deflation.

  procedure Standard_Multiplicity
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the multiplicity structure in standard double precision,
  --   after prompting for a polynomial system and solutions.
  --   The verbose level is in vrblvl.

  procedure DoblDobl_Multiplicity 
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the multiplicity structure in double double precision,
  --   after prompting for a polynomial system and solutions.
  --   The verbose level is in vrblvl.

  procedure QuadDobl_Multiplicity
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the multiplicity structure in quad double precision,
  --   after prompting for a polynomial system and solutions,
  --   The verbose level is in vrblvl.

  procedure Multiplicity_Structure 
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and then launches the corresponding driver
  --   to compute the multiplicity structure.
  --   The verbose level is in vrblvl.

  procedure Solution_Scanner
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Handles option #0, to deal with huge solution lists.

  procedure Main ( infilename,outfilename : in string;
                    verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   The main procedure to verify solutions takes on input two strings,
  --   which are the arguments passed at the command line: the names of
  --   the input and output file.  The last argument is the verbose level.
  --   A menu is shown and an answer is prompted.

end Main_Verification;
