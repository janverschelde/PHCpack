with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Coefficient_Circuits;      use Standard_Coefficient_Circuits;

package Standard_Refiner_Circuits is

-- DESCRIPTION :
--   Offers procedures to apply Newton's method on coefficient circuits
--   in double precision, with condition tables and cluster report.

  procedure Show_Parameters
              ( maxit : in natural32;
                tolres,tolerr,tolsing : in double_float );

  -- DESCRIPTION :
  --   Displays the values of the parameters to run several steps
  --   with Newton's method.

  -- ON ENTRY :
  --   maxit        maximum number of iterations;
  --   tolres       tolerance on the residual;
  --   tolerr       tolerance on the forward error;
  --   tolsing      tolerance on a singularity.

  procedure Set_Parameters
              ( maxit : out natural32;
                tolres,tolerr,tolsing : out double_float );

  -- DESCRIPTION :
  --   Sets default values for the parameters, and then runs
  --   and interactive loop to allow the user to modify the values.

  -- ON RETURN :
  --   maxit        maximum number of iterations;
  --   tolres       tolerance on the residual;
  --   tolerr       tolerance on the forward error;
  --   tolsing      tolerance on a singularity.

  procedure Monitor_Report
              ( idx : in integer32; fail,isreal : in boolean;
                err,rco,res,wgt,tolsing : in double_float );

  -- DESCRIPTION :
  --   Writes one line to screen with a report on a solution,
  --   to monitor the progress on the verification.

  -- ON ENTRY :
  --   idx      index number of the current solution;
  --   fail     true if Newton's method failed, false otherwise;
  --   isreal   true if real, false otherwise (only if not fail);
  --   err      forward error;
  --   rco      estimate for the inverse condition number;
  --   res      residual;
  --   wgt      weight of the coordinates;
  --   tolsing  tolerance to decide whether a solution is singular.

  procedure Run ( s : in Link_to_System; sols : in Solution_List;
                  maxit : in natural32;
                  tolres,tolerr,tolsing : in double_float;
                  cntfail,cntreal,cntcmplx : out natural32;
                  cntregu,cntsing,cntclus : out natural32;
                  t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                  verbose : in boolean; vrb : in integer32 := 0 );
  procedure Run ( file : in file_type;
                  s : in Link_to_System; sols : in Solution_List;
                  maxit : in natural32;
                  tolres,tolerr,tolsing : in double_float;
                  cntfail,cntreal,cntcmplx : out natural32;
                  cntregu,cntsing,cntclus : out natural32;
                  t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                  verbose : in boolean; vrb : in integer32 := 0 );
  procedure Inlined_Run
                ( s : in Link_to_System; sols : in Solution_List;
                  maxit : in natural32;
                  tolres,tolerr,tolsing : in double_float;
                  cntfail,cntreal,cntcmplx : out natural32;
                  cntregu,cntsing,cntclus : out natural32;
                  t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                  verbose : in boolean; vrb : in integer32 := 0 );
  procedure Inlined_Run
                ( file : in file_type;
                  s : in Link_to_System; sols : in Solution_List;
                  maxit : in natural32;
                  tolres,tolerr,tolsing : in double_float;
                  cntfail,cntreal,cntcmplx : out natural32;
                  cntregu,cntsing,cntclus : out natural32;
                  t_err,t_rco,t_res : out Standard_Natural_Vectors.Vector;
                  verbose : in boolean; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the Newton verification on s and the solutions in sols.
  --   The Inlined_Run uses inlined linear algebra solving.

  -- ON ENTRY :
  --   file       (optional) to write all output to file;
  --   s          a system of circuits;
  --   sols       approximate solutions;
  --   maxit      maximum number of Newton steps;
  --   tolres     tolerance on the residual;
  --   tolerr     tolerance on the forward error;
  --   tolsing    tolerance on the inverse condition number.
  --   t_err      vector of range 0..15;
  --   t_rco      vector of range 0..15;
  --   t_res      vector of range 0..15;
  --   verbose    if true, then one linear is written for every solution;
  --   vrb        the verbose level.

  -- ON RETURN :
  --   cntfail    the number of failures;
  --   cntreal    the number of real solutions;
  --   cntcmplx   the number of nonreal solutions;
  --   cntregu    the number of regular solutions;
  --   cntsing    the number of singular solutions;
  --   cntclus    the number of clustered solutions;
  --   t_err      frequency table with forward errors;
  --   t_rco      frequency table with condition numbers;
  --   t_res      frequency table with residuals.

  procedure Run ( s : in Link_to_System; sols : in Solution_List;
                  vrb : in integer32 := 0 );
  procedure Inlined_Run
                ( s : in Link_to_System; sols : in Solution_List;
                  vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs several steps of Newton's method on the system s,
  --   starting at the solutions in sols.
  --   For each solution writes one line to screen.
  --   This version does not write to file.
  --   The verbose level is given in vrb.
  --   The Inlined_Run uses inlined linear system solving.

  procedure Run ( file : in file_type;
                  s : in Link_to_System; sols : in Solution_List;
                  vrb : in integer32 := 0 );
  procedure Inlined_Run
                ( file : in file_type;
                  s : in Link_to_System; sols : in Solution_List;
                  vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs several steps of Newton's method on the system s,
  --   starting at the solutions in sols.
  --   Writes the output to file.
  --   The verbose level is given in vrb.
  --   The Inlined_Run uses inlined linear system solving.

  procedure Main ( vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts a user for a polynomial system with solutions
  --   and then runs the refiners.
  --   The verbose level is in vrb.

  procedure Main ( infilename,outfilename : in string;
                   vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   The strings infilename and outfilename contain the names
  --   of the input and output files.  If empty, then the user
  --   will be prompted for those file names,
  --   otherwise, the names will be used for input and output.
  --   The verbose level is in vrb.

end Standard_Refiner_Circuits;
