with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Solutions;

package DoblDobl_Sampling_Laurent_Machine is

-- DESCRIPTION :
--   The sampling machine provides the lower level interface to the
--   path trackers and root refiners in double double precision,
--   for witness sets define by Laurent polynomial systems.
--   The routines are listed in the order in which they should be used:
--   initialize, tune, sample, and clear.

-- INITIALIZATION AND TUNING :

  procedure Initialize ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Initialization of the internal states of the machine with the
  --   embedded original polynomial system.

  -- ON ENTRY :
  --   ep         embedded system.

  function Embedded_System return DoblDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Return the embedded system.

  procedure Default_Tune_Sampler ( level : in natural32 );
  procedure Default_Tune_Sampler ( file : in file_type; level : in natural32 );

  -- DESCRIPTION :
  --   Adjusts the settings of the path trackers with the defaulted level,
  --   the higher level, the more conservative the path tracking.
  --   If a file is provided as argument, then the selected settings
  --   will be written to the file.

  procedure Interactive_Tune_Sampler;
  procedure Interactive_Tune_Sampler ( file : in file_type );

  -- DESCRIPTION :
  --   Allows the user to interactively tune the settings of the path
  --   trackers, starting to display the current default settings.
  --   If a file is provided as argument, then the selected settings
  --   will be written to the file.

  procedure Default_Tune_Refiner;
  procedure Default_Tune_Refiner ( file : in file_type );

  -- DESCRIPTION :
  --   Adjusts the settings of the root refiner.
  --   The settings are written to the file if it is given as argument.

  procedure Interactive_Tune_Refiner;
  procedure Interactive_Tune_Refiner ( file : in file_type );

  -- DESCRIPTION :
  --   Allows the user to tune the settings of the root refiner,
  --   starting with the display of the defaulted settings.
  --   The settings are written to the file if it is given as argument.

-- MODIFIER :

  procedure Change_Slices ( hyp : DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   After the machine is initialized with a system, only samples
  --   starting at solution that satisfy that system can be made.
  --   This routine changes the slices of the system that was initialized,
  --   thus allowing samples to be made from solutions on other slices.

  -- REQUIRED :
  --   The machine has been initialized with a system.
  --   The number of slices matches the number of slices of that system.

-- DIAGNOSTICS :

  function Satisfies ( sol : in DoblDobl_Complex_Solutions.Solution )
                     return boolean;
  function Satisfies ( file : in file_type;
                       sol : in DoblDobl_Complex_Solutions.Solution )
                     return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution satisfies the precision requirements.
  --   When the file is supplied, then its diagnostics are written.
  --   This function allows to monitor the quality of the samples.

-- SAMPLERS : 

  procedure Sample ( startsol : in DoblDobl_Complex_Solutions.Solution;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     newsol : out DoblDobl_Complex_Solutions.Solution );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     startsol : in DoblDobl_Complex_Solutions.Solution;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     newsol : out DoblDobl_Complex_Solutions.Solution );

  -- DESCRIPTION :
  --   Given a start solution and slices, a new sample is computed.

  -- REQUIRED :
  --   The machine is initialized and the sampler has been tuned.
  --   The startsol satisfies the system given on initialization.

  -- ON ENTRY :
  --   file          for diagnostics and intermediate output;
  --   full_output   if true, then all diagnostics from path tracking
  --                 and refining are given, otherwise output is summary;
  --   startsol      current sample to start from;
  --   newhyp        hyperplane sections for the new sample.

  -- ON RETURN : 
  --   newsol     a new sample.

  procedure Sample ( startsols : in DoblDobl_Complex_Solutions.Solution_List;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     newsols : out DoblDobl_Complex_Solutions.Solution_List);
  procedure Sample ( startsols : in DoblDobl_Complex_Solutions.Solution_List;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     gamma : in DoblDobl_Complex_Vectors.Vector;
                     newsols : out DoblDobl_Complex_Solutions.Solution_List);
  procedure Sample ( file : in file_type;
                     startsols : in DoblDobl_Complex_Solutions.Solution_List;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     newsols : out DoblDobl_Complex_Solutions.Solution_List);
  procedure Sample ( file : in file_type;
                     startsols : in DoblDobl_Complex_Solutions.Solution_List;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     gamma : in DoblDobl_Complex_Vectors.Vector;
                     newsols : out DoblDobl_Complex_Solutions.Solution_List);

  -- DESCRIPTION :
  --   Generates a new list of samples for the given hyperplane sections.

  -- REQUIRED :
  --   The machine is initialized and the sampler has been tuned.
  --   The startsols satisfy the system given on initialization.

  -- ON ENTRY :
  --   file          for diagnostics and intermediate output;
  --   startsols     samples to start from;
  --   newhyp        hyperplane sections for the new sample;
  --   gamma         vector of gamma constants for homotopy,
  --                 by default, it will be generated at random.

  -- ON RETURN : 
  --   newsols       new samples.

  procedure Sample_with_Stop
                   ( startsols : in DoblDobl_Complex_Solutions.Solution_List;
                     x : in DoblDobl_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     newsols : out DoblDobl_Complex_Solutions.Solution_List);
  procedure Sample_with_Stop
                   ( file : in file_type;
                     startsols : in DoblDobl_Complex_Solutions.Solution_List;
                     x : in DoblDobl_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in DoblDobl_Complex_VecVecs.VecVec;
                     newsols : out DoblDobl_Complex_Solutions.Solution_List);

  -- DESCRIPTION :
  --   Generates a new list of samples for the given hyperplane sections.
  --   Stops when x is found, i.e.: there is a solution which equals x
  --   within the given tolerance.

  -- REQUIRED :
  --   The machine is initialized and the sampler has been tuned.
  --   The startsols satisfy the system given on initialization.

  -- ON ENTRY :
  --   file         for diagnostics and intermediate output;
  --   startsols    samples to start from;
  --   x            vectors that needs to be found;
  --   tol          tolerance to decide equality of numbers;
  --   newhyp       hyperplane sections for the new sample.

  -- ON RETURN : 
  --   newsols      new samples.

-- DEALLOCATION :

  procedure Clear;

  -- DESCRIPTION :
  --   Destruction of the internal states of the machine,
  --   with deallocation of the occupied memory resources.

end DoblDobl_Sampling_Laurent_Machine;
