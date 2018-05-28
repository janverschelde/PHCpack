with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Drivers_for_DEMiCs_Algorithm is

-- DESCRIPTION :
--   This package offers an interface to the DEMiCs Algorithm,
--   for dynamic enumeration of all mixed cells.

  procedure DEMiCs_Algorithm_Info;

  -- DESCRIPTION :
  --   Displays information on the DEMiCs Algorithm to screen.

  procedure BlackBox_DEMiCs_Algorithm
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
              lif : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 );
  procedure BlackBox_DEMiCs_Algorithm
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
              lif : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 );

  -- DESCRIPTION :
  --   Calls the DEMiCs algorithm to compute the mixed volume of the
  --   Newton polytopes spanned by the supports of the system p.

  -- ON ENTRY :
  --   p        an ordinary polynomial or a Laurent polynomial system.

  -- ON RETURN :
  --   mix      type of mixture of the supports;
  --   lif      the lifted supports;
  --   mcc      a regular mixed-cell configuration;
  --   mv       the mixed volume.

  procedure Process_DEMiCs_Output
              ( file : in file_type;
                mcc2file,ranstart : in boolean;
                subfile,ranfile : in out file_type;
                dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                timer : in Timing_Widget );

  -- DESCRIPTION :
  --   Processes the output of the DEMiCs algorithm.

  -- ON ENTRY :
  --   file     for output;
  --   mcc2file if the mixed subdivision has to be written on file;
  --   ranstart if a random coefficient system is needed;
  --   subfile  file opened for output if mcc2file;
  --   ranfile  file to write the random coefficient system on;
  --   dim      dimension of the points;
  --   mix      type of mixture of the supports;
  --   sup      points in the supports;
  --   timer    contains timings of the DEMiCs algorithm.

  procedure Run_DEMiCs_Algorithm
              ( file : in file_type; nt : in integer32;
                mcc2file,ranstart : in boolean;
                subfile,ranfile : in out file_type;
                dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists );

  -- DESCRIPTION :
  --   Runs the DEMiCs Algorithm, followed by the application of polyhedral
  --   homotopies to solve a random coefficient system,
  --   either with no multitasking or with a 2-stage pipeline.

  -- ON ENTRY :
  --   file     for output;
  --   nt       number of tasks, if < 2, then no multitasking;
  --   mcc2file if the mixed subdivision has to be written on file;
  --   ranstart if a random coefficient system is needed;
  --   subfile  file opened for output if mcc2file;
  --   ranfile  file to write the random coefficient system on;
  --   dim      dimension of the points;
  --   mix      type of mixture of the supports;
  --   sup      points in the supports.

  procedure Driver_for_DEMiCs_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Interactive driver to the MixedVol Algorithm,
  --   as called by phc -m.
  --   All output is written to file.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   nt       number of tasks, 0 for no multitasking;
  --   p        a polynomial system.

end Drivers_for_DEMiCs_Algorithm;
