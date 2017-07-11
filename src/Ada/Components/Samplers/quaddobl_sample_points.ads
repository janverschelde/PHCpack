with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers; 
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Solutions;

package QuadDobl_Sample_Points is

-- DESCRIPTION :
--   This package provides an abstraction of and operations on points
--   sampled from a component of solutions to a polynomial system,
--   in quad double complex arithmetic.

-- DATA STRUCTURES :

  type QuadDobl_Sample is private;

  type Array_of_QuadDobl_Samples is
    array ( integer32 range <> ) of QuadDobl_Sample;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine.

-- CREATORS and COPY :

  function Create ( sol : in QuadDobl_Complex_Solutions.Solution;
                    hyp : in QuadDobl_Complex_VecVecs.VecVec )
                  return QuadDobl_Sample;

  -- DESCRIPTION :
  --   A sample is created from a solution to an embedded
  --   polynomial system and from hyperplane sections.

  -- WARNING : data sharing occurs, use Copy below to duplicate everything.

  procedure Copy ( s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample );

  -- DESCRIPTION :
  --   Makes a deep copy of the internal data involved.

-- SAMPLERS :

  procedure Sample ( s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample );
  procedure Parallel_Sample
                   ( s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample );
  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s1 : in QuadDobl_Sample; s2 : out QuadDobl_Sample );

  -- DESCRIPTION :
  --   A new sample is generated and stored in s2.
  --   As an option, diagnostics are written to a file, completely if
  --   full_output is true and only summarized if full_output is false.
  --   With Parallel_Sample, samples are generated on parallel slices.

  -- REQUIRED :
  --   The sampling machine is initialized and has been tuned.

  procedure Sample_on_Slices
                   ( s1 : in QuadDobl_Sample;
                     hyp : in QuadDobl_Complex_VecVecs.VecVec;
                     s2 : out QuadDobl_Sample );
  procedure Sample_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s1 : in QuadDobl_Sample;
                     hyp : in QuadDobl_Complex_VecVecs.VecVec;
                     s2 : out QuadDobl_Sample );

  -- DESCRIPTION :
  --   Instead of generating random slices, a new multi-precision
  --   sample will be computed on the given slices.

-- SELECTORS :

  function Number_of_Variables ( s : QuadDobl_Sample ) return integer32;

  -- DESCRIPTION :
  --   Returns the attribute n of the solution representing the point.
  --   For an empty sample, 0 is returned.

  function Number_of_Slices ( s : QuadDobl_Sample ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of hyperplane sections of the sample.
  --   For an empty sample, 0 is returned.

  function Sample_Point ( s : QuadDobl_Sample )
                        return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns the sample point represented as a solution.
  --   Warning, data sharing!

  -- REQUIRED : s /= null.

  function Hyperplane_Sections ( s : QuadDobl_Sample )
                               return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the hyperplane sections that cut out the sample.
  --   Warning, data sharing!

  -- REQUIRED : s /= null.

-- DESTRUCTORS :

  procedure Shallow_Clear ( s : in out QuadDobl_Sample );

  -- DESCRIPTION :
  --   Only deallocation of the encapsulating data structures.

  procedure Deep_Clear ( s : in out QuadDobl_Sample );

  -- DESCRIPTION :
  --   Deallocation of all the occupied memory resources.

  -- WARNING : mind data sharing with a deep clear.

private

  type QuadDobl_Sample_Rep;
  type QuadDobl_Sample is access QuadDobl_Sample_Rep;

end QuadDobl_Sample_Points;
