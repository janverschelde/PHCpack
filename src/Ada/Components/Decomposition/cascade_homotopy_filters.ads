with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Cascade_Homotopy_Filters is

-- DESCRIPTION :
--   A cascade homotopy defines a sequence of homotopies to compute
--   to compute generic points on all components of the solution set.
--   Versions of the cascade homotopies are provided
--   1) for ordinary and for Laurent polynomial systems, and
--   2) in double, double double, and quad double precision.
--   The six versions are combined with various levels of output:
--   1) the output is written to one single file, or
--   2) each superwitness set is written to a separate file.
--   The junk points are removed from the superwitness sets
--   by the application of homotopy membership filters

  procedure Witness_Filter
               ( outfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );

  -- DESCRIPTION :
  --   Calculates candidate witness points on every component,
  --   starting at the component of dimension k,
  --   and removes the junk points from the superwitness sets,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   outfile   file for intermediate results and diagnostics;
  --   nt        number of tasks, set to zero for no tasking;
  --   ep        an embedded polynomial system for the top dimension,
  --             with as many slack variables and random hyperplanes added
  --             as the value of topdim;
  --   sols      solutions to the system ep, for all slack variables,
  --             with zero an nonzero values;
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   rcotol    tolerance on inverse of the condition number estimate,
  --             to bypass the homotopy membership test for regular solutions;
  --   restol    tolerance on the residual;
  --   homtol    tolerance for the homotopy membership test.

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );

  -- DESCRIPTION :
  --   Calculates candidate witness points on every component,
  --   starting at the component of dimension k,
  --   and removes the junk points from the superwitness sets,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   outfile   file for intermediate results and diagnostics;
  --   resfile   file for final results:
  --               1) system at each level, for topdim down to 0; and
  --               2) solutions with zero slack variables;
  --   nt        number of tasks, set to zero for no tasking;
  --   ep        an embedded polynomial system for the top dimension,
  --             with as many slack variables and random hyperplanes added
  --             as the value of topdim;
  --   sols      solutions to the system ep, for all slack variables,
  --             with zero an nonzero values;
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   rcotol    tolerance on inverse of the condition number estimate,
  --             to bypass the homotopy membership test for regular solutions;
  --   restol    tolerance on the residual;
  --   homtol    tolerance for the homotopy membership test.

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );
  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float );

  -- DESCRIPTION :
  --   This witness filter writes the witness sets to files,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   name      file name for the top embedded system;
  --   outfile   file for all intermediate and final results;  
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        an embedded polynomial system for the top dimension,
  --             with as many slack variables and random hyperplane added
  --             as the value of topdim;
  --   sols      solutions to the system ep for the top dimension,
  --             with all zero and nonzero values for the slack variables;
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   rcotol    tolerance on inverse of the condition number estimate,
  --             to bypass the homotopy membership test for regular solutions;
  --   restol    tolerance on the residual;
  --   homtol    tolerance for the homotopy membership test.

-- SILENT VERSIONS OF FILTERS :

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );

  -- DESCRIPTION :
  --   Silent version of the witness generate to compute witness supersets,
  --   and of the witness filter to remove the junk points from the supersets,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        an embedded polynomial system for the top dimension,
  --             with as many slack variables and random hyperplanes added
  --             as the top dimension topdim;
  --   sols      solutions to the system ep with all values for the slack
  --             variables, the zero and the nonzero values;
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   rcotol    tolerance on inverse of the condition number estimate,
  --             to bypass the homotopy membership test for regular solutions;
  --   restol    tolerance on the residual;
  --   homtol    tolerance for the homotopy membership test.

  -- ON RETURN :
  --   embsys    sequence of embedded polynomial systems;
  --   esols0    witness points at each dimension;
  --   pathcnts  table with path counts during the cascade homotopies;
  --   filtcnts  counts of the witness points after each junk removal;
  --   castms    CPU time at each stage in the cascade homotopy;
  --   filtms    CPU time at each stage of the homotopy membership filters;
  --   totcas    total CPU time for running the cascade homotopies;
  --   totfil    total CPU time for the homotopy memberhip filters;
  --   alltime   the total elapsed CPU time.

-- FILTERS WITH OUTPUT TO CALLBACK PROCEDURES :

  procedure Witness_Filter_Callback
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Filter_Callback
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Filter_Callback
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Filter_Callback
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Filter_Callback
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Filter_Callback
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );

  -- DESCRIPTION :
  --   Silent version of the witness generate to compute witness supersets,
  --   and of the witness filter to remove the junk points from the supersets,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        an embedded polynomial system for the top dimension,
  --             with as many slack variables and random hyperplanes added
  --             as the top dimension topdim;
  --   sols      solutions to the system ep with all values for the slack
  --             variables, the zero and the nonzero values;
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   rcotol    tolerance on inverse of the condition number estimate,
  --             to bypass the homotopy membership test for regular solutions;
  --   restol    tolerance on the residual;
  --   homtol    tolerance for the homotopy membership test.

  -- ON RETURN :
  --   embsys    sequence of embedded polynomial systems;
  --   esols0    witness points at each dimension;
  --   pathcnts  table with path counts during the cascade homotopies;
  --   filtcnts  counts of the witness points after each junk removal;
  --   castms    CPU time at each stage in the cascade homotopy;
  --   filtms    CPU time at each stage of the homotopy membership filters;
  --   totcas    total CPU time for running the cascade homotopies;
  --   totfil    total CPU time for the homotopy memberhip filters;
  --   alltime   the total elapsed CPU time.

-- FILTERS WITH OUTPUT TO FILE :

  procedure Witness_Filter
              ( file : in file_type; nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( file : in file_type; nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( file : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( file : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( file : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );
  procedure Witness_Filter
              ( file : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration );

  -- DESCRIPTION :
  --   Reporting version of the witness generate to compute witness supersets,
  --   and of the witness filter to remove the junk points from the supersets,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        an embedded polynomial system for the top dimension,
  --             with as many slack variables and random hyperplanes added
  --             as the top dimension topdim;
  --   sols      solutions to the system ep with all values for the slack
  --             variables, the zero and the nonzero values;
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   rcotol    tolerance on inverse of the condition number estimate,
  --             to bypass the homotopy membership test for regular solutions;
  --   restol    tolerance on the residual;
  --   homtol    tolerance for the homotopy membership test.

  -- ON RETURN :
  --   embsys    sequence of embedded polynomial systems;
  --   esols0    witness points at each dimension;
  --   pathcnts  table with path counts during the cascade homotopies;
  --   filtcnts  counts of the witness points after each junk removal;
  --   castms    CPU time at each stage in the cascade homotopy;
  --   filtms    CPU time at each stage of the homotopy membership filters;
  --   totcas    total CPU time for running the cascade homotopies;
  --   totfil    total CPU time for the homotopy memberhip filters;
  --   alltime   the total elapsed CPU time.

end Cascade_Homotopy_Filters;
