with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers; 
with Standard_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;

package Sample_Points is

-- DESCRIPTION :
--   This package provides an abstraction of and operations on points
--   sampled from a component of solutions to an ordinary polynomial system,
--   or to a Laurent polynomial system.

-- DATA STRUCTURES :

  type Standard_Sample is private;
  type Multprec_Sample is private;

  type Array_of_Standard_Samples is
    array ( integer32 range <> ) of Standard_Sample;
  type Array_of_Multprec_Samples is
    array ( integer32 range <> ) of Multprec_Sample;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine.

-- CREATORS and COPY :

  function Create ( sol : in Standard_Complex_Solutions.Solution;
                    hyp : in Standard_Complex_VecVecs.VecVec )
                  return Standard_Sample;
  function Create ( sol : in Multprec_Complex_Solutions.Solution;
                    hyp : in Multprec_Complex_VecVecs.VecVec )
                  return Multprec_Sample;

  -- DESCRIPTION :
  --   A sample is created from a solution to an embedded
  --   polynomial system and from hyperplane sections.

  -- WARNING : data sharing occurs, use Copy below to duplicate everything.

  procedure Copy ( s1 : in Standard_Sample; s2 : out Standard_Sample );
  procedure Copy ( s1 : in Multprec_Sample; s2 : out Multprec_Sample );

  -- DESCRIPTION :
  --   Makes a deep copy of the internal data involved.

-- SAMPLERS and REFINERS :

  procedure Sample ( s1 : in Standard_Sample; s2 : out Standard_Sample );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Standard_Sample );
  procedure Sample ( s1 : in Standard_Sample; s2 : out Multprec_Sample );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Multprec_Sample );
  procedure Parallel_Sample
                   ( s1 : in Standard_Sample; s2 : out Standard_Sample );
  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Standard_Sample );
  procedure Parallel_Sample
                   ( s1 : in Standard_Sample; s2 : out Multprec_Sample );
  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample; s2 : out Multprec_Sample );

  -- DESCRIPTION :
  --   A new sample is generated and, when s2 is a multi-precision sample,
  --   refined with multi-precision arithmetic.
  --   As an option, diagnostics are written to a file, completely if
  --   full_output is true and only summarized if full_output is false.
  --   With Parallel_Sample, samples are generated on parallel slices.

  -- REQUIRED :
  --   The sampling machine is initialized and has been tuned.

  procedure Sample_on_Slices
                   ( s1 : in Standard_Sample;
                     hyp : in Standard_Complex_VecVecs.VecVec;
                     s2 : out Standard_Sample );
  procedure Sample_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample;
                     hyp : in Standard_Complex_VecVecs.VecVec;
                     s2 : out Standard_Sample );
  procedure Sample_on_Slices
                   ( s1 : in Standard_Sample;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample );
  procedure Sample_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s1 : in Standard_Sample;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample );

  -- DESCRIPTION :
  --   Instead of generating random slices, a new multi-precision
  --   sample will be computed on the given slices.

  procedure Refine ( s1 : in out Standard_Sample; s2 : out Multprec_Sample );
  procedure Refine ( file : in file_type; full_output : in boolean;
                     s1 : in out Standard_Sample; s2 : out Multprec_Sample );
  procedure Refine ( s : in out Multprec_Sample );
  procedure Refine ( file : in file_type; full_output : in boolean;
                     s : in out Multprec_Sample );

  -- DESCRIPTION :
  --   Refines a sample with multi-precision arithmethic.
  --   As an option, diagnostics are written to a file.
  --   If full_output, then the errors and residuals of all intermediate
  --   refinement steps are listed, otherwise only a summary is written.

  -- REQUIRED :
  --   The sampling machine is initialized and has been tuned.

  procedure Refine_on_Slices
                   ( s1 : in out Standard_Sample;
                     hyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample );
  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s1 : in out Standard_Sample;
                     hyp : in Multprec_Complex_VecVecs.VecVec;
                     s2 : out Multprec_Sample );
  procedure Refine_on_Slices
                   ( s : in out Multprec_Sample );
  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s : in out Multprec_Sample );

  -- DESCRIPTION :
  --   Does multi-precision root refinement on the slices given in
  --   the parameter hyp or in the hyperplane sections of the sample s,
  --   working directly on the multi-precision coefficients.

-- SELECTORS :

  function Number_of_Variables ( s : Standard_Sample ) return integer32;
  function Number_of_Variables ( s : Multprec_Sample ) return integer32;

  -- DESCRIPTION :
  --   Returns the attribute n of the solution representing the point.
  --   For an empty sample, 0 is returned.

  function Number_of_Slices ( s : Standard_Sample ) return integer32;
  function Number_of_Slices ( s : Multprec_Sample ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of hyperplane sections of the sample.
  --   For an empty sample, 0 is returned.

  function Sample_Point ( s : Standard_Sample )
                        return Standard_Complex_Solutions.Solution;
  function Sample_Point ( s : Multprec_Sample )
                        return Multprec_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns the sample point represented as a solution.
  --   Warning, data sharing!

  -- REQUIRED : s /= null.

  function Hyperplane_Sections ( s : Standard_Sample )
                               return Standard_Complex_VecVecs.VecVec;
  function Hyperplane_Sections ( s : Multprec_Sample )
                               return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the hyperplane sections that cut out the sample.
  --   Warning, data sharing!

  -- REQUIRED : s /= null.

  function Refined ( s : Standard_Sample ) return Multprec_Sample;

  -- DESCRIPTION :
  --   If a refinement for the sample exists, then a link to the refined 
  --   sample is returned, otherwise a zero pointer is returned.

  function Original ( s : Multprec_Sample ) return Standard_Sample;

  -- DESCRIPTION :
  --   Returns the original sample from which the sample was refined.

-- DESTRUCTORS :

  procedure Shallow_Clear ( s : in out Standard_Sample );
  procedure Shallow_Clear ( s : in out Multprec_Sample );

  -- DESCRIPTION :
  --   Only deallocation of the encapsulating data structures.

  procedure Deep_Clear ( s : in out Standard_Sample );
  procedure Deep_Clear ( s : in out Multprec_Sample );

  -- DESCRIPTION :
  --   Deallocation of all the occupied memory resources.

  -- WARNING : mind data sharing with a deep clear.

private

  type Standard_Sample_Rep;
  type Standard_Sample is access Standard_Sample_Rep;
  type Multprec_Sample_Rep;
  type Multprec_Sample is access Multprec_Sample_Rep;

end Sample_Points;
