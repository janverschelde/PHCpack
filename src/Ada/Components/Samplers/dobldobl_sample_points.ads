with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers; 
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Solutions;

package DoblDobl_Sample_Points is

-- DESCRIPTION :
--   This package provides an abstraction of and operations on points
--   sampled from a component of solutions to a polynomial system,
--   in double double complex arithmetic.

-- DATA STRUCTURES :

  type DoblDobl_Sample is private;

  type Array_of_DoblDobl_Samples is
    array ( integer32 range <> ) of DoblDobl_Sample;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine.

-- CREATORS and COPY :

  function Create ( sol : in DoblDobl_Complex_Solutions.Solution;
                    hyp : in DoblDobl_Complex_VecVecs.VecVec )
                  return DoblDobl_Sample;

  -- DESCRIPTION :
  --   A sample is created from a solution to an embedded
  --   polynomial system and from hyperplane sections.

  -- WARNING : data sharing occurs, use Copy below to duplicate everything.

  procedure Copy ( s1 : in DoblDobl_Sample; s2 : out DoblDobl_Sample );

  -- DESCRIPTION :
  --   Makes a deep copy of the internal data involved.

-- SAMPLERS :

  procedure Sample ( s1 : in DoblDobl_Sample; s2 : out DoblDobl_Sample );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     s1 : in DoblDobl_Sample; s2 : out DoblDobl_Sample );
  procedure Parallel_Sample
                   ( s1 : in DoblDobl_Sample; s2 : out DoblDobl_Sample );
  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s1 : in DoblDobl_Sample; s2 : out DoblDobl_Sample );

  -- DESCRIPTION :
  --   A new sample is generated and stored in s2.
  --   As an option, diagnostics are written to a file, completely if
  --   full_output is true and only summarized if full_output is false.
  --   With Parallel_Sample, samples are generated on parallel slices.

  -- REQUIRED :
  --   The sampling machine is initialized and has been tuned.

  procedure Sample_on_Slices
                   ( s1 : in DoblDobl_Sample;
                     hyp : in DoblDobl_Complex_VecVecs.VecVec;
                     s2 : out DoblDobl_Sample );
  procedure Sample_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s1 : in DoblDobl_Sample;
                     hyp : in DoblDobl_Complex_VecVecs.VecVec;
                     s2 : out DoblDobl_Sample );

  -- DESCRIPTION :
  --   Instead of generating random slices, a new multi-precision
  --   sample will be computed on the given slices.

-- SELECTORS :

  function Number_of_Variables ( s : DoblDobl_Sample ) return integer32;

  -- DESCRIPTION :
  --   Returns the attribute n of the solution representing the point.
  --   For an empty sample, 0 is returned.

  function Number_of_Slices ( s : DoblDobl_Sample ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of hyperplane sections of the sample.
  --   For an empty sample, 0 is returned.

  function Sample_Point ( s : DoblDobl_Sample )
                        return DoblDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns the sample point represented as a solution.
  --   Warning, data sharing!

  -- REQUIRED : s /= null.

  function Hyperplane_Sections ( s : DoblDobl_Sample )
                               return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the hyperplane sections that cut out the sample.
  --   Warning, data sharing!

  -- REQUIRED : s /= null.

-- DESTRUCTORS :

  procedure Shallow_Clear ( s : in out DoblDobl_Sample );

  -- DESCRIPTION :
  --   Only deallocation of the encapsulating data structures.

  procedure Deep_Clear ( s : in out DoblDobl_Sample );

  -- DESCRIPTION :
  --   Deallocation of all the occupied memory resources.

  -- WARNING : mind data sharing with a deep clear.

private

  type DoblDobl_Sample_Rep;
  type DoblDobl_Sample is access DoblDobl_Sample_Rep;

end DoblDobl_Sample_Points;
