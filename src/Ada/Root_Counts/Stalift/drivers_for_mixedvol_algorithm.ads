with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Drivers_for_MixedVol_Algorithm is

-- DESCRIPTION :
--   This package offers an interface to the MixedVol Algorithm.
--   The polyhedral homotopies come in several flavors:
--   (0) silent or reporting intermediate results to file;
--   (1) precision: standard double, double double, or quad double;
--   (2) multithreaded or not, with nt tasks or serial when nt = 0.

  procedure MixedVol_Algorithm_Info;

  -- DESCRIPTION :
  --   Displays information on the MixedVol Algorithm to screen.

  procedure Mixed_Volume_Computation
              ( n : in integer32; 
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float; r : out integer32;
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                sub : out Mixed_Subdivision; mixvol : out natural32;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system and computes its
  --   mixed volume.  A regular mixed-cell configuration is returned.

  -- ON ENTRY :
  --   n        number of variables and number of supports;
  --   s        supports of a polynomial system;
  --   stlb     lifting bound to use for stable mixed volumes,
  --            equals 0.0 if no stable mixed volumes are requested;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            must be used for the Hermite normal form;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   r        number of different supports;
  --   mix      array of range 1..r with occurences of supports;
  --   perm     permutation of the supports, starting to count at zero,
  --            if r < n then the equations must be permuted;
  --   sub      a regular mixed-cell configuration;
  --   mixvol   mixed volume of the polytopes spanned by the supports.
 
  procedure Random_Coefficient_System
              ( nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies polyhedral continuation methods to solve a random coefficient
  --   start system, using a regular mixed-cell configuration.
  --   This procedure produces no intermediate output.

  -- ON ENTRY :
  --   nt       number of tasks to run the polyhedral homotopies,
  --            for sequential runs, set nt to zero;
  --   n        ambient dimension of the system;
  --   mix      type of mixture;
  --   ls       lifted supports of the polynomial system,
  --            contains only those points that actually occur in a cell;
  --   sub      a regular mixed-cell configuration;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            must be used for the Hermite normal form;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   q        a random coefficient system with same supports as in s;
  --   qsols    solutions to q, as many as its mixed volume.
 
  procedure Random_Coefficient_System
              ( file : in file_type;
                nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( file : in file_type;
                nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( file : in file_type;
                nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( file : in file_type;
                nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( file : in file_type;
                nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Random_Coefficient_System
              ( file : in file_type;
                nt,n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                sub : in Mixed_Subdivision;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies polyhedral continuation methods to solve a random coefficient
  --   start system, using a regular mixed-cell configuration.
  --   Intermediate output is written to file, according to the output level
  --   requested by the user.

  -- ON ENTRY :
  --   file     file for intermediate output and diagnostics;
  --   nt       number of tasks to run the polyhedral homotopies,
  --            for sequential runs, set nt to zero;
  --   n        ambient dimension of the system;
  --   mix      type of mixture;
  --   ls       lifted supports of the polynomial system,
  --            contains only those points that actually occur in a cell;
  --   sub      a regular mixed-cell configuration;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            must be used for the Hermite normal form;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   q        a random coefficient system with same supports as in s;
  --   qsols    solutions to q, as many as its mixed volume.

  procedure Polyhedral_Homotopies
              ( cfile,qfile : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( cfile,qfile : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( cfile,qfile : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( cfile,qfile : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( cfile,qfile : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( cfile,qfile : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Interactive driver to the polyhedral homotopies to create
  --   a random coefficient start system to solve the system p.

  -- ON ENTRY :
  --   cfile    the file to write a regular mixed-cell configuration on;
  --   qfile    the file to write a random coefficient start system on;
  --   p        a polynomial system.

  procedure Polyhedral_Continuation
              ( file : in file_type; nt : in integer32;
                stable,contrep : in boolean;
                n,r : in integer32; stlb : in double_float;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                p : in Standard_Complex_Poly_Systems.Poly_Sys; 
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                sub,mcc,stbmcc : in Mixed_Subdivision;
               -- mcc,stbmcc : in Mixed_Subdivision;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Continuation
              ( file : in file_type; nt : in integer32;
                stable,contrep : in boolean;
                n,r : in integer32; stlb : in double_float;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys; 
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                sub,mcc,stbmcc : in Mixed_Subdivision;
               -- mcc,stbmcc : in Mixed_Subdivision;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Continuation
              ( file : in file_type; nt : in integer32;
                stable,contrep : in boolean;
                n,r : in integer32; stlb : in double_float;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys; 
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                sub,mcc,stbmcc : in Mixed_Subdivision;
               -- mcc,stbmcc : in Mixed_Subdivision;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Continuation
              ( file : in file_type; nt : in integer32;
                stable,contrep : in boolean;
                n,r : in integer32; stlb : in double_float;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                p : in Standard_Complex_Laur_Systems.Laur_Sys; 
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                sub,mcc,stbmcc : in Mixed_Subdivision;
               -- mcc,stbmcc : in Mixed_Subdivision;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Continuation
              ( file : in file_type; nt : in integer32;
                stable,contrep : in boolean;
                n,r : in integer32; stlb : in double_float;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys; 
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                sub,mcc,stbmcc : in Mixed_Subdivision;
               -- mcc,stbmcc : in Mixed_Subdivision;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Continuation
              ( file : in file_type; nt : in integer32;
                stable,contrep : in boolean;
                n,r : in integer32; stlb : in double_float;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys; 
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                sub,mcc,stbmcc : in Mixed_Subdivision;
               -- mcc,stbmcc : in Mixed_Subdivision;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs polyhedral continuation to solve a random coefficient
  --   system, using the given regular mixed-cell configuration.

  -- ON ENTRY :
  --   file     to write intermediate output and diagnostics;
  --   nt       number of tasks to execute the polyhedral homotopies,
  --            for sequential code, set nt to zero;
  --   stable   if zero component solutions have to be computed as well;
  --   contrep  indicates whether extra output is needed;
  --   n        ambient dimension;
  --   r        number of different supports;
  --   stlb     lifting bound used for stable mixed volumes;
  --   mix      type of mixture;
  --   perm     permutation of the equations in case r < n;
  --   p        the original polynomial system;
  --   s        supports of the system p;
  --   sub      all mixed cells;
  --   mcc      a regular mixed-cell configuration;
  --   stbmcc   extra stable mixed cells if stable is on;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            must be used for the Hermite normal form;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   q        a random coefficient system with same supports as in s;
  --   qsols    solutions of q.
  --   qsols0   solutions with zero components of q.

  procedure Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                stable,misufile,ranstart,contrep : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mv,smv,tmv : out natural32;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                stable,misufile,ranstart,contrep : in boolean;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                mv,smv,tmv : out natural32;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                stable,misufile,ranstart,contrep : in boolean;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                mv,smv,tmv : out natural32;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                stable,misufile,ranstart,contrep : in boolean;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mv,smv,tmv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                stable,misufile,ranstart,contrep : in boolean;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mv,smv,tmv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                stable,misufile,ranstart,contrep : in boolean;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mv,smv,tmv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Interactive driver to the polyhedral homotopies to create
  --   a random coefficient start system to solve the system p.

  -- ON ENTRY :
  --   file     to write timings and other results;
  --   cfile    the file to write a regular mixed-cell configuration on,
  --            will be used only if misufile equals true;
  --   qfile    the file to write a random coefficient start system on,
  --            will be used only if ranstart equals true;
  --   nt       number of tasks to use in polyhedral homotopies,
  --            if 0, then plain sequential code is applied;
  --   stable   indicates whether the user wants stable mixed volumes;
  --   misufile indicates whether to write the mixed-cell configuration
  --            to a separate file;
  --   ranstart indicates whether user wants a random coefficient system;
  --   contrep  indicates whether output is wanted during continuation;
  --   p        a polynomial system;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            must be used for the Hermite normal form;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   mv       mixed volume of the tuple of Newton polytopes
  --            spanned by the supports of p;
  --   smv      stable mixed volume counts also solutions with zeroes;
  --   tmv      total mixed volume includes volumes of spurious cells;
  --   q        if ranstart, then q is a random coefficient system,
  --            with as many solutions as the mixed volume;
  --   qsols    solution to a random coefficient system, if ranstart;
  --   qsols0   if stable, then solutions with zero components.

  procedure Ask_for_Stable_and_Cells_File
              ( stable,onfile : out boolean; file : in out file_type );

  -- DESCRIPTION :
  --   Asks the user if the stable mixed volume is wanted and whether
  --   the mixed-cell configuration needs to be written to a separate file.

  -- ON RETURN :
  --   stable   true if the stable mixed volume is wanted,
  --            false if the user does not wants the stable mixed volume;
  --   onfile   true if the mixed-cell configuration needs to written
  --            to a separate file, false otherwise;
  --   file     opened for output if onfile is true.

  procedure Ask_only_if_Stable_and_Cells_File
              ( stable : in out boolean; onfile : out boolean;
                file : in out file_type );

  -- DESCRIPTION :
  --   Asks the user if the stable mixed volume is wanted, 
  --   only if stable is true on entry, and if the mixed-cell
  --   configuration needs to be written to a separate file.

  -- ON ENTRY :
  --   stable   if false, then no question about stable mixed volume is asked,
  --            as is useful when the system is a genuine Laurent system and
  --            therefore, because of negative exponents, a stable mixed
  --            volumes does not make sense, otherwise,
  --            if true, then the user is asked about stable mixed volumes.

  -- ON RETURN :
  --   stable   true if the stable mixed volume is wanted,
  --            false if the user does not wants the stable mixed volume;
  --   onfile   true if the mixed-cell configuration needs to written
  --            to a separate file, false otherwise;
  --   file     opened for output if onfile is true.

  procedure Driver_for_MixedVol_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                byebye,nostart : in boolean;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Driver_for_MixedVol_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                byebye,nostart : in boolean;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Driver_for_MixedVol_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                byebye,nostart : in boolean;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Driver_for_MixedVol_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                byebye,nostart : in boolean;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Driver_for_MixedVol_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                byebye,nostart : in boolean;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );
  procedure Driver_for_MixedVol_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                byebye,nostart : in boolean;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Interactive driver to the MixedVol Algorithm,
  --   as called by phc -m.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   nt       number of tasks to use in polyhedral homotopies,
  --            if nt = 0, then the plain sequential code is used;
  --   p        a polynomial system;
  --   byebye   true if a bye-bye message needs to be given;
  --   nostart  true if the user does not want a start system;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            must be used for the Hermite normal form;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   q        start system (if requested by user);
  --   qsols    solutions to q (if mv > 0);
  --   qsols0   solutions with zeroes;
  --   mv       mixed volume;
  --   smv      stable mixed volume;
  --   tmv      total mixed volume;

end Drivers_for_MixedVol_Algorithm;
