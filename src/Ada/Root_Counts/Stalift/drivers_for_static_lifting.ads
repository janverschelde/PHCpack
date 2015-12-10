with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Integer_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions;

package Drivers_for_Static_Lifting is

  procedure Static_Lifting_Info;

  -- DESCRIPTION :
  --   Displays information on static lifting on screen.

  procedure Integer_Volume_Computation
               ( file : in file_type;
                 n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 compmisu : in boolean;
                 lifpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32 );

  -- DESCRIPTION :
  --   Computes the volumes of the mixed cells in mixsub.
  --   If some cells are not fine yet, then they will be refined.

  -- ON ENTRY :
  --   file      for output of the mixed-cell configuration;
  --   n         ambient dimension, length of the point vectors;
  --   mix       type of mixture;
  --   compmisu  true if the mixed-cell configuration was computed,
  --             false if the user provided the mixed-cell configuration;
  --   lifpts    lifted points;
  --   mixsub    a regular mixed-cell configuration,
  --             induced by integer-valued lifting.

  -- ON RETURN :
  --   mixsub    if compmisu, then nonfine mixed cells are refined;
  --   mv        sum of the volume of all mixed cells.

  procedure Floating_Volume_Computation
               ( n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32;
                 multprec_hermite : in boolean := false );
  procedure Floating_Volume_Computation
              ( file : in file_type; n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32;
                multprec_hermite : in boolean := false );
  procedure Floating_Volume_Computation
               ( n : in integer32; stlb : in double_float;
                 mix : in Standard_Integer_Vectors.Vector;
                 mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv,smv,tmv : out natural32;
                 multprec_hermite : in boolean := false );
  procedure Floating_Volume_Computation
               ( file : in file_type;
                 n : in integer32; stlb : in double_float;
                 mix : in Standard_Integer_Vectors.Vector;
                 mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv,smv,tmv : out natural32;
                 multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   Writes the mixed-cell configuration to file (if provided)
  --   and computes the volumes of all mixed cells.
  --   Remains silent if no output file is provided.

  -- ON ENTRY :
  --   file      for output of the mixed-cell configuration;
  --   n         ambient dimension, length of the points;
  --   stlb      lifting bound used for the artificial origins,
  --             if zero, then no stable mixed volume was requested;
  --   mix       type of mixture;
  --   mixsub    a regular mixed-cell configuration,
  --             induced by floating-valued lifting.

  -- ON RETURN :
  --   mixsub    may be altered for refinement;
  --   mv        the mixed volume;
  --   smv       stable mixed volume (only if stlb is nonzero);
  --   tmv       total mixed volume (only if stlb is nonzero).

  procedure Integer_Polyhedral_Homotopy_Continuation
               ( file : in file_type; contrep : in boolean;
                 n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 q : in out Poly_Sys; qsols : in out Solution_List;
                 lifted : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision );
  procedure Integer_Polyhedral_Homotopy_Continuation
               ( file : in file_type; contrep : in boolean;
                 n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 q : in out Laur_Sys; qsols : in out Solution_List;
                 lifted : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision );

  -- DESCRIPTION :
  --   Solves a random coefficient (Laurent) polynomial systems, using a
  --   polyhedral homotopy induced by integer-valued lifting function.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   contrep   if true, then there will be output during continuation,
  --             otherwise, the continuation will run silently;
  --   n         ambient dimension = number of equations;
  --   mix       type of mixture;
  --   q         a random coefficient (Laurent) polynomial system;
  --   lifted    the lifted supports of q;
  --   mixsub    a regular mixed-cell configuration.

  -- ON RETURN :
  --   qsols     the solutions of the random coefficient system q.

  procedure Floating_Polyhedral_Homotopy_Continuation
               ( file : in file_type; nt : in integer32; contrep : in boolean;
                 n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 q : in Poly_Sys; qsols : in out Solution_List;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 fltsub : in Floating_Mixed_Subdivisions.Mixed_Subdivision );
  procedure Floating_Polyhedral_Homotopy_Continuation
               ( file : in file_type; nt : in integer32; contrep : in boolean;
                 n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 q : in Laur_Sys; qsols : in out Solution_List;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 fltsub : in Floating_Mixed_Subdivisions.Mixed_Subdivision );

  -- DESCRIPTION :
  --   Solves a random coefficient (Laurent) polynomial system, using a
  --   polyhedral homotopy induced by a floating-point lifting function.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks, if 0, then no multitasking, otherwise
  --             nt tasks will work on polyhedral homotopies;
  --   contrep   if true, then there will be output during continuation,
  --             otherwise, the continuation will be silent;
  --   n         ambient dimension = number of equations;
  --   mix       records the number of occurrences of each support;
  --   q         a random coefficient (Laurent) polynomial system;
  --   lifsup    the lifted supports;
  --   fltsub    a mixed-cell configuration induced by the lifting.

  -- ON RETURN :
  --   qsols     the solutions of the random coefficient system q.

  procedure Driver_for_Mixed_Volume_Computation 
               ( file : in file_type; nt : in integer32;
                 p : in Poly_Sys; byebye : in boolean;
                 q : out Poly_Sys; qsols,qsols0 : out Solution_List;
                 mv,smv,tmv : out natural32 );
  procedure Driver_for_Mixed_Volume_Computation 
               ( file : in file_type; nt : in integer32;
                 p : in Laur_Sys; byebye : in boolean;
                 q : out Laur_Sys; qsols,qsols0 : out Solution_List;
                 mv,smv,tmv : out natural32 );

  -- DESCRIPTION :
  --   Interactive driver for the computation of the mixed volume.

  -- ON ENTRY :
  --   file      output file, must be opened for output;
  --   nt        number of tasks, if 0, then no multitasking is used,
  --             otherwise nt tasks will do polyhedral continuation;
  --   p         a polynomial system;
  --   byebye    if true, then a bye-bye message will appear on screen,
  --             if false, then no bye-bye;

  -- ON RETURN :
  --   q         a start system with randomly choosen coefficients,
  --             which can be used in a coefficient homotopy;
  --   qsols     the solutions of q;
  --   qsols0    solutions of q with zero components;
  --   mv        mixed volume of p and the number of solutions of q;
  --   smv       stable mixed volume of p;
  --   tmv       total mixed volume of entire mixed-cell configuration.

end Drivers_for_Static_Lifting;
