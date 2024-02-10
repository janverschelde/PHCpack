with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package PHCpack_Operations is

-- DESCRIPTION :
--   Provides the most common operations of PHCpack,
--   with data management for polynomial systems and solutions.

  file_okay : boolean := false;    -- if output file defined
  output_file : file_type;         -- output file defined by user

  procedure Define_Output_File;
  procedure Define_Output_File ( name : in string );
  function Is_File_Defined return boolean;
 -- function Retrieve_Output_File return file_type;
  procedure Close_Output_File;

  -- DESCRIPTION :
  --   The "Determine_Output_File" lets the user define an output file,
  --   which can be retrieved by "Retrieve_Output_File".
  --   If "Is_File_Defined" returns true, then it is okay to write to
  --   "Retrieve_Output_File", otherwise "standard_output" is returned.

  procedure Tuned_Continuation_Parameters;
  function Are_Continuation_Parameters_Tuned return boolean;

  -- DESCRIPTION :
  --   The "Tuned_Continuation_Parameters" indicates the values
  --   of the continuation parameters have already been set,
  --   and no automatic tuning of these parameters is necessary.
  --   After calling "Tuned_Continuation_Parameters", the result
  --   of "Are_Continuation_Parameters_Tuned" is true.

-- MANAGING PERSISTENT DATA :
--   the "store" operations make copies of the systems and solutions,
--   while the "retrieve" procedures give pointers to the stored objects.

  procedure Store_Start_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Start_System
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Store_Start_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Start_System
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Store_Start_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Start_System
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Store_Start_System
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Start_Solutions
              ( sols : in Standard_Complex_Solutions.Solution_List );
  procedure Store_Start_Solutions
              ( sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Store_Start_Solutions
              ( sols : in QuadDobl_Complex_Solutions.Solution_List );
  procedure Store_Start_Solutions
              ( sols : in Multprec_Complex_Solutions.Solution_List );
  procedure Store_Target_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Target_System
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Store_Target_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Target_System
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Store_Target_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Target_System
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Store_Target_System
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys );
  procedure Store_Target_Solutions
              ( sols : in Standard_Complex_Solutions.Solution_List );
  procedure Store_Target_Solutions
              ( sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Store_Target_Solutions
              ( sols : in QuadDobl_Complex_Solutions.Solution_List );
  procedure Store_Target_Solutions
              ( sols : in Multprec_Complex_Solutions.Solution_List );

  procedure Retrieve_Start_System
              ( p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Start_System
              ( p : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys );
  procedure Retrieve_Start_System
              ( p : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Start_System
              ( p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys );
  procedure Retrieve_Start_System
              ( p : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Start_System
              ( p : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys );
  procedure Retrieve_Start_System
              ( p : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Start_Solutions
              ( sols : out Standard_Complex_Solutions.Solution_List );
  procedure Retrieve_Start_Solutions
              ( sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Retrieve_Start_Solutions
              ( sols : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Retrieve_Start_Solutions
              ( sols : out Multprec_Complex_Solutions.Solution_List );
  procedure Retrieve_Target_System
              ( p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Target_System
              ( p : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys );
  procedure Retrieve_Target_System
              ( p : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Target_System
              ( p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys );
  procedure Retrieve_Target_System
              ( p : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Target_System
              ( p : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys );
  procedure Retrieve_Target_System
              ( p : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Retrieve_Target_Solutions
              ( sols : out Standard_Complex_Solutions.Solution_List );
  procedure Retrieve_Target_Solutions
              ( sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Retrieve_Target_Solutions
              ( sols : out QuadDobl_Complex_Solutions.Solution_List );
  procedure Retrieve_Target_Solutions
              ( sols : out Multprec_Complex_Solutions.Solution_List );

  procedure Store_Gamma_Constant
              ( gamma : in Standard_Complex_Numbers.Complex_Number );
  procedure Store_Gamma_Constant
              ( gamma : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Store_Gamma_Constant
              ( gamma : in QuadDobl_Complex_Numbers.Complex_Number );
  procedure Store_Gamma_Constant
              ( gamma : in Multprec_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Stores the gamma constant.

  procedure Retrieve_Gamma_Constant
              ( gamma : out Standard_Complex_Numbers.Complex_Number );
  procedure Retrieve_Gamma_Constant
              ( gamma : out DoblDobl_Complex_Numbers.Complex_Number );
  procedure Retrieve_Gamma_Constant
              ( gamma : out QuadDobl_Complex_Numbers.Complex_Number );
  procedure Retrieve_Gamma_Constant
              ( gamma : out Multprec_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Retrieves the gamma constant.
  --   If no gamma constant was stored, the gamma on return is zero.

  procedure Create_Standard_Homotopy;
  procedure Create_Standard_Homotopy
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                pwrt : in natural32 := 2 );
  procedure Create_DoblDobl_Homotopy;
  procedure Create_DoblDobl_Homotopy
              ( gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                pwrt : in natural32 := 2 );
  procedure Create_QuadDobl_Homotopy;
  procedure Create_QuadDobl_Homotopy
              ( gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                pwrt : in natural32 := 2 );
  procedure Create_Multprec_Homotopy;
  procedure Create_Multprec_Homotopy
              ( gamma : in Multprec_Complex_Numbers.Complex_Number;
                pwrt : in natural32 := 2 );

  -- DESCRIPTION :
  --   Creates a linear homotopy between target and start system,
  --   which must be already stored.  If no gamma is provided, then
  --   a random complex constant will be generated.
  --   The default value for the power of t is 2.

  -- REQUIRED : start and target system have been stored.

  procedure Clear_Standard_Homotopy;
  procedure Clear_DoblDobl_Homotopy;
  procedure Clear_QuadDobl_Homotopy;

  -- DESCRIPTION :
  --   Allocated memory for the linear homotopy is cleared.

  procedure Create_Standard_Laurent_Homotopy;
  procedure Create_Standard_Laurent_Homotopy
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                pwrt : in natural32 := 2 );
  procedure Create_DoblDobl_Laurent_Homotopy;
  procedure Create_DoblDobl_Laurent_Homotopy
              ( gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                pwrt : in natural32 := 2 );
  procedure Create_QuadDobl_Laurent_Homotopy;
  procedure Create_QuadDobl_Laurent_Homotopy
              ( gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                pwrt : in natural32 := 2 );

  -- DESCRIPTION :
  --   Creates a linear homotopy between target and start system,
  --   which must be already stored.  If no gamma is provided, then
  --   a random complex constant will be generated.

  -- REQUIRED : start and target system have been stored.

  procedure Clear_Standard_Laurent_Homotopy;
  procedure Clear_DoblDobl_Laurent_Homotopy;
  procedure Clear_QuadDobl_Laurent_Homotopy;

  -- DESCRIPTION :
  --   Allocated memory for the linear homotopy 
  --   for Laurent systems is cleared.

  procedure Standard_Cascade_Homotopy;
  procedure DoblDobl_Cascade_Homotopy;
  procedure QuadDobl_Cascade_Homotopy;

  -- DESCRIPTION :
  --   Creates the homotopy to go down one level in the cascade
  --   to remove one slice from the start system,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The start system has been stored and contains at least one
  --   slack variable and at least one slicing hyperplane at the end.
  --   Otherwise, if there is no start system yet, then the stored 
  --   target system will be used as start system.

  procedure Standard_Cascade_Laurent_Homotopy;
  procedure DoblDobl_Cascade_Laurent_Homotopy;
  procedure QuadDobl_Cascade_Laurent_Homotopy;

  -- DESCRIPTION :
  --   Creates the homotopy to go down one level in the cascade
  --   to remove one slice from the start system,
  --   in standard double, double double, or quad double precision,
  --   for Laurent polynomial systems.

  -- REQUIRED :
  --   The start system has been stored and contains at least one
  --   slack variable and at least one slicing hyperplane at the end.
  --   Otherwise, if there is no start system yet, then the stored 
  --   target system will be used as start system.

  procedure Standard_Diagonal_Homotopy ( a,b : in natural32 );
  procedure DoblDobl_Diagonal_Homotopy ( a,b : in natural32 );
  procedure QuadDobl_Diagonal_Homotopy ( a,b : in natural32 );

  -- DESCRIPTION :
  --   Creates the homotopy to start the cascade in a diagonal homotopy
  --   to intersect two positive dimensional solution sets,
  --   defined by ordinary polynomial systems,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The target and start system stored internally have their symbols
  --   matched and the dimensions are sorted: a >= b.

  -- ON ENTRY :
  --   a        dimension of the first set stored as the target system;
  --   b        dimension of the second set stored as the start system.

  -- ON RETURN :
  --   The start system is the system to start the cascade and
  --   the target system is the system read to perform a cascade on.

  procedure Standard_Diagonal_Laurent_Homotopy ( a,b : in natural32 );
  procedure DoblDobl_Diagonal_Laurent_Homotopy ( a,b : in natural32 );
  procedure QuadDobl_Diagonal_Laurent_Homotopy ( a,b : in natural32 );

  -- DESCRIPTION :
  --   Creates the homotopy to start the cascade in a diagonal homotopy
  --   to intersect two positive dimensional solution sets,
  --   defined by Laurent polynomial systems,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The target and start system stored internally have their symbols
  --   matched and the dimensions are sorted: a >= b.

  -- ON ENTRY :
  --   a        dimension of the first set stored as the target system;
  --   b        dimension of the second set stored as the start system.

  -- ON RETURN :
  --   The start system is the system to start the cascade and
  --   the target system is the system read to perform a cascade on.

  procedure Standard_Diagonal_Cascade_Solutions ( a,b : in natural32 );
  procedure DoblDobl_Diagonal_Cascade_Solutions ( a,b : in natural32 );
  procedure QuadDobl_Diagonal_Cascade_Solutions ( a,b : in natural32 );

  -- DESCRIPTION :
  --   Makes the start solutions to start the cascade homotopy to
  --   intersect the witness sets stored in target and start system
  --   with corresponding witness points in target and start solutions,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The target and start system stored represent witness sets
  --   of dimensions a and b respectively.

  -- ON ENTRY :
  --   a        dimension of the first set stored as the target system;
  --   b        dimension of the second set stored as the start system.

  -- ON RETURN :
  --   The start solutions are replaced by the solutions to start the
  --   diagonal cascade.  The target solutions are cleared.

  procedure Silent_Path_Tracker
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Uses the created homotopy to track one path starting at the
  --   given solution.  The silent version does not write any output.

  -- REQUIRED : the homotopy has been initialized.

  -- ON ENTRY :
  --   ls        solution of the start system.

  -- ON RETURN :
  --   ls        solution at the end of the solution path;
  --   length    measure for the length of the path;
  --   nbstep    number of steps along the path;
  --   nbfail    number of failures along the path;
  --   nbiter    number of corrector iterations along the path;
  --   nbsyst    number of linear systems solved along the path;
  --   crash     true if some exception occurred.

  procedure Silent_Path_Tracker 
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 wnd : out integer32;
                 dir : out Standard_Floating_Vectors.Link_to_Vector;
                 err,length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker 
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 wnd : out integer32;
                 dir : out Double_Double_Vectors.Link_to_Vector;
                 err,length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Silent_Path_Tracker 
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 wnd : out integer32;
                 dir : out Quad_Double_Vectors.Link_to_Vector;
                 err,length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Uses the created homotopy to track one path starting at the
  --   given solution.  The silent version does not write any output.
  --   For the extra parameters wnd, dir, and err to have meaning,
  --   the polyhedral end games must be invoked by setting the
  --   order of the extrapolation number of a positive number.

  -- REQUIRED : the homotopy has been initialized.

  -- ON ENTRY :
  --   ls        solution of the start system.

  -- ON RETURN :
  --   ls        solution at the end of the solution path;
  --   wnd       estimated winding number;
  --   dir       computed direction of the path;
  --   err       error on the estimated direction;
  --   length    measure for the length of the path;
  --   nbstep    number of steps along the path;
  --   nbfail    number of failures along the path;
  --   nbiter    number of corrector iterations along the path;
  --   nbsyst    number of linear systems solved along the path;
  --   crash     true if some exception occurred.

  procedure Reporting_Path_Tracker
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Reporting_Path_Tracker
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Reporting_Path_Tracker
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Uses the created homotopy to track one path starting at the
  --   given solution.  The reporting version writes extra output
  --   to the defined output file.

  -- REQUIRED : the homotopy has been initialized.

  -- ON ENTRY :
  --   ls        solution of the start system.

  -- ON RETURN :
  --   ls        solution at the end of the solution path;
  --   length    measure for the length of the path;
  --   nbstep    number of steps along the path;
  --   nbfail    number of failures along the path;
  --   nbiter    number of corrector iterations along the path;
  --   nbsyst    number of linear systems solved along the path;
  --   crash     true if some exception occurred.

  procedure Reporting_Path_Tracker
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 wnd : out integer32;
                 dir : out Standard_Floating_Vectors.Link_to_Vector;
                 err,length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Reporting_Path_Tracker
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 wnd : out integer32;
                 dir : out Double_Double_Vectors.Link_to_Vector;
                 err,length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Reporting_Path_Tracker
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 wnd : out integer32;
                 dir : out Quad_Double_Vectors.Link_to_Vector;
                 err,length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Uses the created homotopy to track one path starting at the
  --   given solution.  The reporting version writes extra output
  --   to the defined output file.
  --   For the extra parameters wnd, dir, and err to have meaning,
  --   the polyhedral end games must be invoked by setting the
  --   order of the extrapolation number of a positive number.

  -- REQUIRED : the homotopy has been initialized.

  -- ON ENTRY :
  --   ls        solution of the start system.

  -- ON RETURN :
  --   ls        solution at the end of the solution path;
  --   wnd       estimated winding number;
  --   dir       computed direction of the path;
  --   err       error on the estimated direction;
  --   length    measure for the length of the path;
  --   nbstep    number of steps along the path;
  --   nbfail    number of failures along the path;
  --   nbiter    number of corrector iterations along the path;
  --   nbsyst    number of linear systems solved along the path;
  --   crash     true if some exception occurred.

  procedure Silent_Laurent_Path_Tracker
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Silent_Laurent_Path_Tracker
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Silent_Laurent_Path_Tracker
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );

  procedure Reporting_Laurent_Path_Tracker
               ( ls : in Standard_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Reporting_Laurent_Path_Tracker
               ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );
  procedure Reporting_Laurent_Path_Tracker
               ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 length : out double_float;
                 nbstep,nbfail,nbiter,nbsyst : out natural32;
                 crash : out boolean; nbq : in integer32 := 0 );

  -- DESCRIPTION :
  --   Uses the created homotopy to track one path starting at the
  --   given solution.  The reporting version writes extra output
  --   to the defined output file.

  -- REQUIRED : {Standard,DoblDobl,QuadDobl}_Laurent_Homotopy has
  --   been initialized.

  -- ON ENTRY :
  --   ls        solution of the start system.

  -- ON RETURN :
  --   ls        solution at the end of the solution path;
  --   length    measure for the length of the path;
  --   nbstep    number of steps along the path;
  --   nbfail    number of failures along the path;
  --   nbiter    number of corrector iterations along the path;
  --   nbsyst    number of linear systems solved along the path;
  --   crash     true if some exception occurred.

  function Solve_by_Standard_Homotopy_Continuation 
             ( number_of_tasks : natural32;
               vrblvl : integer32 := 0 ) return integer32;
  function Solve_by_DoblDobl_Homotopy_Continuation
             ( number_of_tasks : natural32;
               vrblvl : integer32 := 0 ) return integer32;
  function Solve_by_QuadDobl_Homotopy_Continuation
             ( number_of_tasks : natural32;
               vrblvl : integer32 := 0 ) return integer32;
  function Solve_by_Multprec_Homotopy_Continuation
             ( decimals : natural32;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks the paths starting at the current start solutions.
  --   A normal termination returns 0.  Otherwise, the return value
  --   of this function signals an exception.
  --   If the number_of_tasks equals 0, then no multicore will be used,
  --   otherwise as many tasks as the value of number_of_tasks will be
  --   launched to do the path tracking.
  --   The number of decimal places in the working precision for the
  --   multiprecision version equals the value of decimals.

  function Solve_by_Standard_Laurent_Homotopy_Continuation 
             ( number_of_tasks : natural32;
               vrblvl : integer32 := 0 ) return integer32;
  function Solve_by_DoblDobl_Laurent_Homotopy_Continuation
             ( number_of_tasks : natural32;
               vrblvl : integer32 := 0 ) return integer32;
  function Solve_by_QuadDobl_Laurent_Homotopy_Continuation
             ( number_of_tasks : natural32;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Tracks the paths starting a the stored start solutions,
  --   in standard double, double double, or quad double precision,
  --   on the defined homotopy for Laurent polynomial systems.
  --   Multitasking will be applied if number_of_tasks > 0,
  --   using as many tasks as the value of number_of_tasks.

  procedure Standard_Clear;

  -- DESCRIPTION :
  --   Deallocation of the persistent data for standard homotopies,
  --   target, start system, and the solutions.

  procedure Standard_Laurent_Clear;

  -- DESCRIPTION :
  --   Deallocation of the persistent data for standard Laurent homotopies,
  --   target, start system, and the solutions.

  procedure DoblDobl_Clear;

  -- DESCRIPTION :
  --   Deallocation of the persistent data for double double homotopies,
  --   target, start system, and the solutions.

  procedure DoblDobl_Laurent_Clear;

  -- DESCRIPTION :
  --   Deallocation of the persistent data for double double Laurent
  --   homotopies, target, start system, and the solutions.

  procedure QuadDobl_Clear;

  -- DESCRIPTION :
  --   Deallocation of the persistent data for quad double homotopies,
  --   target, start system, and the solutions.

  procedure QuadDobl_Laurent_Clear;

  -- DESCRIPTION :
  --   Deallocation of the persistent data for quad double Laurent
  --   homotopies, target, start system, and the solutions.

  procedure Multprec_Clear;

  -- DESCRIPTION :
  --   Deallocation of the persistent data for multiprecision homotopies.

end PHCpack_Operations;
