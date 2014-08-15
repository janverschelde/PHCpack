with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Process_io is

-- DESCRIPTION :
--   This package determines the output operations during the continuation.
--   The purpose is to concentrate all trivial output operations which could
--   possibly overload the coding of the continuation process.
--   Moreover, an uniform output format is achieved by this package.

  type output_code is ( nil,s,p,c,sp,sc,pc,spc );

  -- Explanation of the output_code :
  --   nil  :  nothing will be written during the continuation process
  --   s    :  all intermediate solutions are written
  --   p    :  predictor information is written
  --   c    :  corrector information is written
  --   sp, sc, pc and spc are combinations of s, p and c

  procedure Set_Output_Code ( u : in output_code );
 
  -- DESCRIPTION :
  --   Sets the status code for output during continuation.

  function Get_Output_Code return output_code;

  -- DESCRIPTION :
  --   Returns the value of the status code for output.

  function Contains_Output_Code ( u : output_code ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the output code set contains u.
  --   For example "Contains_Output_Code(p)" will be true if the
  --   current output code is p, sp, pc, or spc.

  procedure Write_path ( file : in file_type; n : in natural32 ); 

  -- DESCRIPTION :
  --   The number of the paths is written on file or on standard output.

  procedure Write_block ( file : in file_type; n : in natural32 );

  -- DESCRIPTION :
  --   The block number is written on the output device

  procedure sWrite ( file : in file_type;
                     sol : in Standard_Complex_Solutions.Solution );
  procedure sWrite ( file : in file_type;
                     sol : in DoblDobl_Complex_Solutions.Solution );
  procedure sWrite ( file : in file_type;
                     sol : in QuadDobl_Complex_Solutions.Solution );
  procedure sWrite ( file : in file_type;
                     sol : in Multprec_Complex_Solutions.Solution );

  -- DESCRIPTION :
  --   The solution is written on file or on standard output.

  procedure pWrite ( file : in file_type; step : in double_float; 
                     t : in Standard_Complex_Numbers.Complex_Number );
  procedure pWrite ( file : in file_type; step : in double_float; 
                     t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure pWrite ( file : in file_type; step : in double_float; 
                     t : in QuadDobl_Complex_Numbers.Complex_Number );
  procedure pWrite ( file : in file_type; step : in Floating_Number;
                     t : in Multprec_Complex_Numbers.Complex_Number );

  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in Standard_Complex_Numbers.Complex_Number;
                     sol : in Standard_Complex_Solutions.Solution );
  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in DoblDobl_Complex_Numbers.Complex_Number;
                     sol : in DoblDobl_Complex_Solutions.Solution );
  procedure pWrite ( file : in file_type; step : in double_float;
                     t : in QuadDobl_Complex_Numbers.Complex_Number;
                     sol : in QuadDobl_Complex_Solutions.Solution );
  procedure pWrite ( file : in file_type; step : in Floating_Number;
                     t : in Multprec_Complex_Numbers.Complex_Number;
                     sol : in Multprec_Complex_Solutions.Solution );

  procedure pWrite ( file : in file_type; step_number : in natural32;
                     step : in double_float;
                     t : in Standard_Complex_Numbers.Complex_Number );
  procedure pWrite ( file : in file_type; step_number : in natural32;
                     step : in double_double;
                     t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure pWrite ( file : in file_type; step_number : in natural32;
                     step : in quad_double;
                     t : in QuadDobl_Complex_Numbers.Complex_Number );


  -- DESCRIPTION :
  --   The predictor information is written on file or on standard output.

  -- ON ENTRY :
  --   file          must be opened for output;
  --   step_number   number of the current step;
  --   step          current step size;
  --   t             value of the continuation parameter;
  --   sol           predicted solution.

  procedure cWrite ( file : in file_type;
                     normax,normrx,normaf,normrf : in double_float );
  procedure cWrite ( file : in file_type;
                     normax,normrx,normaf,normrf : in Floating_Number );

  -- DESCRIPTION :
  --   The norm of the correction on x and residual is written.

  -- ON ENTRY :
  --   file       file type, must be created or opened for output,
  --              if not specified, then standard output will be taken;
  --   normax     absolute norm of the correction dx on the solution x: ||dx||;
  --   normrx     relative norm of the correction dx: ||dx||/||x||;
  --   normaf     absolute norm of the residual: ||f(x)||;
  --   normrf     relative norm of the residual: ||f(x)||/||x||.

  procedure cWrite ( file : in file_type;
                     rcond : in double_float; m : in natural32 );
  procedure cWrite ( file : in file_type;
                     rcond : in Floating_Number; m : in natural32 );

  -- DESCRIPTION :
  --   The estimate for the inverse condition number of the Jacobi matrix
  --   is written, jointly with the (estimated) multiplicity of the solution.

  procedure Write_Statistics ( file : in file_type;
                               nstep,nfail,niter,nsyst : in natural32 );

  -- DESCRIPTION :
  --   This procedure writes statistical information after the
  --   computation of parts of the results.

  -- ON ENTRY :
  --   nstep     the number of predictor steps;
  --   nfail     the number of failures;
  --   niter     the number of corrector iterations;
  --   nsyst     the number of linear systems solved.

  procedure Write_Total_Statistics
               ( file : in file_type;
                 tnstep,tnfail,tniter,tnsyst : in natural32 );

  -- DESCRIPTION
  --   This procedure writes statistical information after the 
  --   solution of the problem.

  -- ON ENTRY :
  --   tnstep    the total number of predictor steps;
  --   tnfail    the total number of failures;
  --   tniter    the total number of corrector iterations;
  --   tnsyst    the total number of linear systems solved.

  procedure sWrite_Solutions
               ( file : in file_type;
                 sols : in Standard_Complex_Solutions.Solution_List );
  procedure sWrite_Solutions
               ( file : in file_type;
                 sols : in Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes down the computed solutions on standard output or on file.

end Process_io;
