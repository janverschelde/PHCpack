with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Standard_Continuation_Data;         use Standard_Continuation_Data;
with Exponent_Vectors;                   use Exponent_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Polyhedral_Coefficient_Trackers is

-- DESCRIPTION :
--   In a polyhedral coefficient homotopy, the coefficients are random
--   complex numbers multiplied with some (possibly) high powers of the
--   continuation parameter.  To deal with the numerical instabilities
--   caused by these high powers, an exponentional transformation of
--   on the continuation parameter is applied in the path trackers
--   provided by this package.

  procedure Silent_Track_One_Path
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solu_Info; fail : out boolean  );

  procedure Reporting_Track_One_Path
              ( file : in file_type;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solu_Info; fail : out boolean  );

  -- DESCRIPTION :
  --   Tracks one path, starting at one solution "sol" of start system
  --   corresponding to a mixed cell with inner normal given in "nor".

  -- ON ENTRY :
  --   file     to write all kinds of diagnostics of the path;
  --   hq       evaluable form of coefficient Laurent system;
  --   ctm      work space to hold coefficient for fixed value
  --            of the continuation parameter;
  --   pow      powers of the continuation parameter in the homotopy
  --            defined by one mixed cell;
  --   coeffv   coefficients in the polyhedral coefficient homotopy;
  --   jacmat   coefficient Jacobian matrix;
  --   mulfac   multiplication factors for Jacobian matrix;
  --   sol      information for one start solution.

  -- ON RETURN :
  --   sol      solution info at the end of the path, if not fail;
  --   fail     true if path failed to converged in the allowed #steps.

  procedure Track_Paths_for_Cell
              ( file : in file_type; report,monitor : in boolean;
                lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                ind : in integer32; mic : in Mixed_Cell;
                cnt : in out integer32 );

  -- DESCRIPTION :
  --   Tracks the paths using the polyhedral homotopies defined by one
  --   mixed cell, appending to the list of solutions.

  -- ON ENTRY :
  --   file     for writing diagnostics and solutions;
  --   report   reporting or silent trackers depending on true or false; 
  --   monitor  allow user to monitor progress on scree if true;
  --   lq       random coefficient start system;
  --   ctm      work space for the coefficients of the system;
  --   coeffv   coefficients in the polyhedral coefficient homotopy;
  --   jacmat   coefficient Jacobian matrix;
  --   mulfac   multiplication factors for Jacobian matrix;
  --   mix      type of mixture;
  --   mic      a mixed cell;
  --   ind      number to the current mixed cell;
  --   cnt      number of paths tracked so far;

  -- ON RETURN :
  --   ctm      modified work space;
  --   cnt      updated counter for the number of paths.

  procedure Track_Paths_for_Subdivision
              ( file : in file_type; report,monitor : in boolean;
                lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Tracks the paths using the polyhedral homotopies defined by the
  --   mixed cells in the subdivision.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   report   reporting or silent trackers depending on true or false;
  --   monitor  allows user to monitor progress on screen if true;
  --   lq       a random coefficient Laurent polynomial system;
  --   ls       lifted supports of the polynomial system lq;

  -- ON RETURN :
  --   the solutions of the start system are immediately written to file.

  procedure Track_Paths_for_Subdivision
              ( infile,outfile : in file_type;
                report,monitor : in boolean;
                m : in integer32; lq : in Laur_Sys;
                vs : in Standard_Floating_VecVecs.Array_of_VecVecs;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Tracks the paths using the polyhedral homotopies defined by the
  --   mixed cells on the input file.

  -- REQUIRED :
  --   The current position at the input file must be at the beginning
  --   of a mixed cell.
  
  -- ON ENTRY :
  --   infile   input file with a regular mixed-cell configuration,
  --            properly positioned at the beginning of a mixed cell;
  --   outfile  output file for intermediate results;
  --   report   reporting or silent trackers depending on true or false;
  --   monitor  allows user to track progress on screen if true;
  --   m        total number of mixed cells on the input file;
  --   lq       a random coefficient Laurent polynomial system;
  --   vs       lifted supports of the polynomial system lq;
  --   ls       lifted supports of the polynomial system lq;
  --   hq       polyhedral coefficient homotopy;

  -- ON RETURN :
  --   The solutions to the system lq are immediately written to file.

  procedure Polyhedral_Continuation
              ( file : in file_type; report,monitor : in boolean;
                n : in integer32; mv : in natural32; p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision; q : out Poly_Sys );

  -- DESCRIPTION :
  --   Creates a random coefficient system with the vertices used in the
  --   mixed-cell configuration, constructs a polyhedral homotopy, and
  --   then calls the path trackers to solve the random coefficient system.
 
  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics;
  --   report   reporting or silent trackers depending on true or false;
  --   monitor  allows user to monitor progress on screen if true;
  --   n        ambient dimension before the lifting;
  --   mv       the mixed volume (or total number of expected solutions);
  --   p        polynomial system for which a random coefficient
  --            start system will be constructed;
  --   mix      type of mixture read from infile counts the number
  --            of occurrences of each support in the subdivision.

  -- ON RETURN :
  --   q        a random coefficient start system to solve p;
  --            its solutions are written to file.

  procedure Polyhedral_Continuation
              ( infile,outfile : in file_type;
                report,monitor : in boolean;
                n : in integer32; mv : in natural32;
                mix : in Standard_Integer_Vectors.Vector; q : out Poly_Sys );

  -- DESCRIPTION :
  --   Creates a random coefficient system with the vertices used in the
  --   mixed-cell configuration, constructs a polyhedral homotopy, and
  --   then calls the path trackers to solve the random coefficient system.
  --   This jumpstarting version does not require for all mixed cells
  --   to be in the internal memory at the same time.

  -- REQUIRED :
  --   The regular mixed-cell configuration on the input file is in
  --   its labeled presentation, with the lifted sets in front.

  -- ON ENTRY :
  --   infile   input file with a regular mixed-cell configuration,
  --            positioned at the start of the lifted support sets;
  --   outfile  output file for intermediate results and diagnostics;
  --   report   reporting or silent trackers depending on true or false;
  --   monitor  allows user to monitor progress on screen if true;
  --   n        ambient dimension before the lifting;
  --   mv       the mixed volume (or total number of expected solutions);
  --   p        polynomial system for which a random coefficient
  --            start system will be constructed;
  --   mix      type of mixture read from infile counts the number
  --            of occurrences of each support in the subdivision.

  -- ON RETURN :
  --   q        a random coefficient start system to solve p;
  --            its solutions are written to file.
  
end Polyhedral_Coefficient_Trackers;
