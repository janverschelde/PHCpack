with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;

package Sampling_Machine is

-- DESCRIPTION :
--   The sampling machine provides the lower level interface to the
--   path trackers and root refiners.  One could view this as a driver
--   to the facilities available in continuation library.
--   The routines are listed in the order in which they should be used:
--   initialize, tune, sample (evt tune), refine (evt tune), and clear.

-- INITIALIZATION AND TUNING :

  procedure Initialize ( ep : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Initialize ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                         mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                         k : in integer32; size : in natural32 );

  -- DESCRIPTION :
  --   Initialization of the internal states of the machine with the
  --   embedded original polynomial system(s).  If mp is not provided,
  --   the no multi-precision refinement will be possible.

  -- ON ENTRY :
  --   ep         embedded system without slices;
  --   mp         multi-precision version of the original system;
  --   k          number of hyperplane sections added to ep;
  --   size       size of the multi-precision numbers.

  function Embedded_System return Standard_Complex_Poly_Systems.Poly_Sys;
  function Original_System return Multprec_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   These two functions above return the embedded and orginal system
  --   with multi-precision coefficients.  These functions are useful
  --   for restriction to the linear span of components.

  procedure Initialize_Restricted
                ( ep : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Initialize_Restricted
                ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  k : in integer32; size : natural32 );

  -- DESCRIPTION :
  --   Initializes the machine with a polynomial system restricted
  --   to the span of a component.  The sampling machine does not
  --   really have to know this.  All that happens is that there two
  --   different system can be present at the same time.  Once an
  --   initialization with a restricted system has been done, then
  --   that system will be used for refinement, until it is cleared.

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
  procedure Default_Tune_Refiner ( size : in natural32 );
  procedure Default_Tune_Refiner ( file : in file_type; size : in natural32 );

  -- DESCRIPTION :
  --   Adjusts the settings of the root refiner for the size of the numbers.
  --   The settings are written to the file if it is given as argument.
  --   If no size is given, then standard arithmetic is assumed.

  procedure Interactive_Tune_Refiner;
  procedure Interactive_Tune_Refiner ( file : in file_type );
  procedure Interactive_Tune_Refiner ( size : in natural32 );
  procedure Interactive_Tune_Refiner
                ( file : in file_type; size : in natural32 );

  -- DESCRIPTION :
  --   Allows the user to tune the settings of the root refiner, starting
  --   with the display of the defaulted settings for the size of numbers.
  --   The settings are written to the file if it is given as argument.
  --   If no size is given, then standard arithmetic is assumed.

-- MODIFIER :

  procedure Change_Slices ( hyp : Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   After the machine is initialized with a system, only samples
  --   starting at solution that satisfy that system can be made.
  --   This routine changes the slices of the system that was initialized,
  --   thus allowing samples to be made from solutions on other slices.
  --   This procedure does not apply to the restricted system and should
  --   not be applied in conjunction with multi-precision refinement.

  -- REQUIRED :
  --   The machine has been initialized with a system.
  --   The number of slices matches the number of slices of that system.

-- DIAGNOSTICS :

  function Satisfies ( sol : in Standard_Complex_Solutions.Solution )
                     return boolean;
  function Satisfies ( file : in file_type;
                       sol : in Standard_Complex_Solutions.Solution )
                     return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution satisfies the precision requirements.
  --   When the file is supplied, then its diagnostics are written.
  --   This function allows to monitor the quality of the samples.

-- SAMPLERS : 

  procedure Sample ( startsol : in Standard_Complex_Solutions.Solution;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsol : out Standard_Complex_Solutions.Solution );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     startsol : in Standard_Complex_Solutions.Solution;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsol : out Standard_Complex_Solutions.Solution );

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

  procedure Sample ( startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List);
  procedure Sample ( startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     gamma : in Standard_Complex_Vectors.Vector;
                     newsols : out Standard_Complex_Solutions.Solution_List);
  procedure Sample ( file : in file_type;
                     startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List);
  procedure Sample ( file : in file_type;
                     startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     gamma : in Standard_Complex_Vectors.Vector;
                     newsols : out Standard_Complex_Solutions.Solution_List);

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
                   ( startsols : in Standard_Complex_Solutions.Solution_List;
                     x : in Standard_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List);
  procedure Sample_with_Stop
                   ( file : in file_type;
                     startsols : in Standard_Complex_Solutions.Solution_List;
                     x : in Standard_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List);

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

-- REFINERS :

  procedure Refine ( stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution;
                     mphyp : out Multprec_Complex_VecVecs.VecVec );
  procedure Refine ( file : in file_type; full_output : in boolean;
                     stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution;
                     mphyp : out Multprec_Complex_VecVecs.VecVec );
  procedure Refine ( mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in out Multprec_Complex_VecVecs.VecVec );
  procedure Refine ( file : in file_type; full_output : in boolean;
                     mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in out Multprec_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes a multi-precision refinement of a standard precision sample,
  --   or refines a given multi-precision sample.

  -- REQUIRED :
  --   The machine is initialized and the refiner has been tuned.

  -- ON ENTRY :
  --   file          for diagnostics and intermediate output;
  --   full_output   if true, then the errors and residuals of all
  --                 intermediate refining steps are written,
  --                 otherwise, only a summary is written to the file;
  --   stsol         current sample in standard machine precision;
  --   sthyp         hyperplane sections with the current sample.

  -- ON RETURN : 
  --   mpsol         refined solution with multi-precision arithmetic;
  --   mphyp         multi-precision representation of the hyperplane sections.

  procedure Refine_on_Slices
                   ( stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution );
  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution );
  procedure Refine_on_Slices
                   ( mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in Multprec_Complex_VecVecs.VecVec );
  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in Multprec_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Instead of converting the representation of the slices with
  --   standard precision into multi-precision, the "Refine_on_Slices"
  --   allow to refine on given slices with multi-precision numbers.

  -- REQUIRED :
  --   The machine is initialized and the refiner has been tuned.

  -- ON ENTRY :
  --   file          for diagnostics and intermediate output;
  --   full_output   if true, then the errors and residuals of all
  --                 intermediate refining steps are written,
  --                 otherwise, only a summary is written to the file;
  --   stsol         current sample in standard machine precision;
  --   sthyp         hyperplane sections with the current sample;
  --   mphyp         multi-precision representation of the hyperplane sections.

  -- ON RETURN : 
  --   mpsol         refined solution with multi-precision arithmetic.

-- DEALLOCATION :

  procedure Clear;
  procedure Clear_Restricted;

  -- DESCRIPTION :
  --   Destruction of the internal states of the machine, with deallocation
  --   of the occupied memory resources.
  --   Both operations corresponds to the states set by their corresponding
  --   initialize operations.

end Sampling_Machine;
