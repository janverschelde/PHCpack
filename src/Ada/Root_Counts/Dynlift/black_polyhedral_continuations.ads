with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Integer_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions;

package Black_Polyhedral_Continuations is

-- DESCRIPTION :
--   This package provides wrappers to the polyhedral continuation methods.

  procedure Black_Box_Polyhedral_Continuation
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Polyhedral_Continuation
               ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Creates and solves a random coefficient start system, based on 
  --   a regular mixed-cell configuration, induced by integer lifting.

  -- ON ENTRY :
  --   p         a (Laurent) polynomial system;
  --   mix       type of mixture;
  --   lifsup    lifted supports of the system;
  --   mixsub    regular mixed-cell configuration;
  --   mv        mixed volume;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   q         random coefficient start system;
  --   qsols     solutions of q.

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out DoblDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out QuadDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0
                   : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0
                   : in out DoblDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0
                   : in out QuadDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Creates and solves a random coefficient start system, based on 
  --   a regular mixed-cell configuration, induced by floating lifting.

  -- ON ENTRY :
  --   nt        number of tasks for the polyhedral continuation,
  --             if 0, then sequential execution;
  --   p         polynomial system;
  --   mix       type of mixture;
  --   stlb      lifting used for the artificial origin;
  --   lifsup    lifted supports of the system;
  --   orgmcc    a regular mixed-cell configuration with only
  --             the original mixed cells, without artificial origin;
  --   stbmcc    a regular mixed-cell configuration, with only
  --             the extra stable mixed cells, with aritifical origin;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   q         random coefficient start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q with zero components.

  procedure Read_System
              ( file : in out file_type; filename : in string;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Attempts to read a polynomial system from file.

  procedure Polyhedral_Solver
               ( file : in file_type; nt : in natural32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Calls the MixedVol algorithm and solves a random coefficient system.

  procedure Main ( nt : in natural32; infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Defines phc -b -m, calls MixedVol for a mixed cell configuration
  --   and solve a random coefficient system, written to the output file.
  
  -- ON ENTRY :
  --   nt            the number of tasks for the polyhedral path trackers,
  --                 set nt to zero for the sequential code;
  --   infilename    the name of the input file;
  --   outfilename   the name of the output file.

end Black_Polyhedral_Continuations;
