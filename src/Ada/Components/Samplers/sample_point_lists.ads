with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Generic_Lists;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;
with Sample_Points;                      use Sample_Points;

package Sample_Point_Lists is

-- DESCRIPTION :
--   This package provides lists of points sampled from components of
--   solutions of a polynomial system.  A grid of samples is an array
--   of lists of samples.

-- DATA STRUCTURES :

  package Lists_of_Standard_Samples is new Generic_Lists(Standard_Sample);
  type Standard_Sample_List is new Lists_of_Standard_Samples.List;
  package Lists_of_Multprec_Samples is new Generic_Lists(Multprec_Sample);
  type Multprec_Sample_List is new Lists_of_Multprec_Samples.List;

  type Array_of_Standard_Sample_Lists is
    array ( integer32 range <> ) of Standard_Sample_List;
  type Array_of_Multprec_Sample_Lists is
    array ( integer32 range <> ) of Multprec_Sample_List;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine.

-- CREATORS :

  function Create ( sols : Standard_Complex_Solutions.Solution_List;
                    hyps : Standard_Complex_VecVecs.VecVec )
                  return Standard_Sample_List;
  function Create ( sols : Multprec_Complex_Solutions.Solution_List;
                    hyps : Multprec_Complex_VecVecs.VecVec )
                  return Multprec_Sample_List;

  -- DESCRIPTION :
  --   Returns the list of samples, where all solutions in the list
  --   sols all have the same hyperplane sections hyps.

-- SAMPLERS and REFINERS :

  procedure Sample ( s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Standard_Sample_List );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Standard_Sample_List );
  procedure Sample ( s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Multprec_Sample_List );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Multprec_Sample_List );
  procedure Parallel_Sample 
                   ( s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Standard_Sample_List );
  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Standard_Sample_List );
  procedure Parallel_Sample
                   ( s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Multprec_Sample_List );
  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s : in Standard_Sample; m : in natural32;
                     samples,samples_last : in out Multprec_Sample_List );

  -- DESCRIPTION :
  --   Generates m additional samples starting from s and appends those
  --   to the samples list.  Intermediate output may be written on file,
  --   completely if full_output is true and only summarized otherwise.
  --   With Parallel_Sample, the slices are all parallel to each other.

  -- REQUIRED :
  --   The sampling machine is initialized and tuned.

  procedure Sample ( s1 : in Standard_Sample_List;
                     hyps : in Standard_Complex_VecVecs.VecVec;
                     s2,s2_last : in out Standard_Sample_List );
  procedure Sample ( s1 : in Standard_Sample_List;
                     hyps : in Standard_Complex_VecVecs.VecVec;
                     s2,s2_last : in out Multprec_Sample_List );
  procedure Sample_on_Slices
                   ( s1 : in Standard_Sample_List;
                     sthyps : in Standard_Complex_VecVecs.VecVec;
                     mphyps : in Multprec_Complex_VecVecs.VecVec;
                     s2,s2_last : in out Multprec_Sample_List );
  procedure Sample ( file : in file_type;
                     s1 : in Standard_Sample_List;
                     hyps : in Standard_Complex_VecVecs.VecVec;
                     s2,s2_last : in out Standard_Sample_List );
  procedure Sample_with_Stop
                   ( s1 : in Standard_Sample_List;
                     x : in Standard_Complex_Vectors.Vector;
                     tol : in double_float;
                     hyps : in Standard_Complex_VecVecs.VecVec;
                     s2,s2_last : in out Standard_Sample_List );

  -- DESCRIPTION :
  --   Starting at the samples in s1, new samples on the slices hyps
  --   are added to the list s2.  In the "with_stop", the sampling
  --   stops when the vector x has been found.

  procedure Refine ( s1 : in out Standard_Sample_List;
                     s2,s2_last : in out Multprec_Sample_List );
  procedure Refine ( file : in file_type; full_output : in boolean;
                     s1 : in out Standard_Sample_List;
                     s2,s2_last : in out Multprec_Sample_List );
  procedure Refine ( samples : in out Multprec_Sample_List );
  procedure Refine ( file : in file_type; full_output : in boolean;
                     samples : in out Multprec_Sample_List );

  -- DESCRIPTION :
  --   Refines a list of samples with multi-precision arithmethic.
  --   As an option, diagnostics are written to a file.
  --   If full_output, then the errors and residuals of all refinement
  --   steps are listed, otherwise, only a summary is written to file.

  -- REQUIRED :
  --   The sampling machine is initialized and has been tuned.

  procedure Refine_on_Slices
                   ( s1 : in out Standard_Sample_List;
                     hyps : in Multprec_Complex_VecVecs.VecVec;
                     s2,s2_last : in out Multprec_Sample_List );
  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     s1 : in out Standard_Sample_List;
                     hyps : in Multprec_Complex_VecVecs.VecVec;
                     s2,s2_last : in out Multprec_Sample_List );
  procedure Refine_on_Slices ( samples : in out Multprec_Sample_List );
  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     samples : in out Multprec_Sample_List );

  -- DESCRIPTION :
  --   The multi-precision root refinement works directly on the
  --   slices given with the parameter hyps or in the hyperplane
  --   sections of the samples.

-- SAMPLING MEMBERSHIP TEST :

  procedure Membership_Test
                ( s1 : in Standard_Sample_List; 
                  x : in Standard_Complex_Vectors.Vector;
                  tol : in double_float; isin : out natural32;
                  s2,s2_last : in out Standard_Sample_List );

  -- DESCRIPTION :
  --   Performs a sample test to decide whether a point belongs to
  --   a solution component.
  -- REQUIRED : not Is_Null(s1).

  -- ON ENTRY :
  --   s1         all generic points on positive dimensional component;
  --   x          point to test whether lies on that component;
  --   tol        tolerance to decide whether element is zero or not.

  -- ON RETURN :
  --   isin       0 if x does not belong to any of the components,
  --              otherwise this is the index of the matching sample in s;
  --   s2         list of new samples, computed for membership test,
  --              if isin = k /= 0, then x is k-th element in s2;
  --   s2_last    is pointer to last element in the list newsp.

-- SELECTORS :

  function Is_In ( s : Standard_Sample_List; tol : double_float;
                   x : Standard_Complex_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the element in s that matches with x.
  --   If x does not belong to s, then 0 is returned.

  function Is_In ( s : Standard_Sample_List; tol : double_float;
                   spt : Standard_Sample ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the element in s that matches with spt,
  --   within the given tolerance.  Returns 0 if spt does not belong to s.

  function Sample_Points ( s : Standard_Sample_List )
                         return Standard_Complex_Solutions.Solution_List;
  function Sample_Points ( s : Multprec_Sample_List )
                         return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the list of sample points, represented as solutions.

  function Select_Sub_List ( l : Standard_Sample_List;
                             f : Standard_Natural_Vectors.Vector )
                           return Standard_Sample_List;

  -- DESCRIPTION :
  --   Selects from the list only those samples whose position
  --   occurs in the vector f of labeled witness points.

  -- REQUIRED :
  --   The list must have sufficiently many samples and
  --   the points in f are sorted in increasing order.

  function Select_Sub_Grid ( grid : Array_of_Standard_Sample_Lists;
                             f : Standard_Natural_Vectors.Vector )
                           return Array_of_Standard_Sample_Lists;

  -- DESCRIPTION :
  --   Selects from the grid only those samples from the lists whose
  --   position occurs in the vector f of labeled witness points.

  -- REQUIRED :
  --   The lists in the grid have sufficiently many samples
  --   and the points in f are sorted in increasing order.

-- DESTRUCTORS :

  procedure Shallow_Clear ( samples : in out Standard_Sample_List );
  procedure Shallow_Clear ( samples : in out Multprec_Sample_List );
  procedure Shallow_Clear ( samples : in out Array_of_Standard_Sample_Lists );
  procedure Shallow_Clear ( samples : in out Array_of_Multprec_Sample_Lists );

  -- DESCRIPTION :
  --   Deallocates only the encapsulating list structure.

  procedure Deep_Clear ( samples : in out Standard_Sample_List );
  procedure Deep_Clear ( samples : in out Multprec_Sample_List );
  procedure Deep_Clear ( samples : in out Array_of_Standard_Sample_Lists );
  procedure Deep_Clear ( samples : in out Array_of_Multprec_Sample_Lists );

  -- DESCRIPTION :
  --   Deallocates all occupied memory by the structures.

end Sample_Point_Lists;
