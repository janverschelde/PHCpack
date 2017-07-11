with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Generic_Lists;
with Standard_Natural_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Solutions;
with QuadDobl_Sample_Points;             use QuadDobl_Sample_Points;

package QuadDobl_Sample_Lists is

-- DESCRIPTION :
--   This package provides lists of points sampled from components of
--   solutions of a polynomial system, in double double precision.
--   A grid of samples is an array of lists of samples.

-- DATA STRUCTURES :

  package Lists_of_QuadDobl_Samples is new Generic_Lists(QuadDobl_Sample);
  type QuadDobl_Sample_List is new Lists_of_QuadDobl_Samples.List;

  type Array_of_QuadDobl_Sample_Lists is
    array ( integer32 range <> ) of QuadDobl_Sample_List;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine.

-- CREATORS :

  function Create ( sols : QuadDobl_Complex_Solutions.Solution_List;
                    hyps : QuadDobl_Complex_VecVecs.VecVec )
                  return QuadDobl_Sample_List;

  -- DESCRIPTION :
  --   Returns the list of samples, where all solutions in the list
  --   sols all have the same hyperplane sections hyps.

-- SAMPLERS and REFINERS :

  procedure Sample ( s : in QuadDobl_Sample; m : in natural32;
                     samples,samples_last : in out QuadDobl_Sample_List );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     s : in QuadDobl_Sample; m : in natural32;
                     samples,samples_last : in out QuadDobl_Sample_List );
  procedure Parallel_Sample 
                   ( s : in QuadDobl_Sample; m : in natural32;
                     samples,samples_last : in out QuadDobl_Sample_List );
  procedure Parallel_Sample
                   ( file : in file_type; full_output : in boolean;
                     s : in QuadDobl_Sample; m : in natural32;
                     samples,samples_last : in out QuadDobl_Sample_List );

  -- DESCRIPTION :
  --   Generates m additional samples starting from s and appends those
  --   to the samples list.  Intermediate output may be written on file,
  --   completely if full_output is true and only summarized otherwise.
  --   With Parallel_Sample, the slices are all parallel to each other.

  -- REQUIRED :
  --   The sampling machine is initialized and tuned.

  procedure Sample ( s1 : in QuadDobl_Sample_List;
                     hyps : in QuadDobl_Complex_VecVecs.VecVec;
                     s2,s2_last : in out QuadDobl_Sample_List );
  procedure Sample ( file : in file_type;
                     s1 : in QuadDobl_Sample_List;
                     hyps : in QuadDobl_Complex_VecVecs.VecVec;
                     s2,s2_last : in out QuadDobl_Sample_List );
  procedure Sample_with_Stop
                   ( s1 : in QuadDobl_Sample_List;
                     x : in QuadDobl_Complex_Vectors.Vector;
                     tol : in double_float;
                     hyps : in QuadDobl_Complex_VecVecs.VecVec;
                     s2,s2_last : in out QuadDobl_Sample_List );

  -- DESCRIPTION :
  --   Starting at the samples in s1, new samples on the slices hyps
  --   are added to the list s2.  In the "with_stop", the sampling
  --   stops when the vector x has been found.

  -- REQUIRED :
  --   The sampling machine is initialized and has been tuned.

-- SAMPLING MEMBERSHIP TEST :

  procedure Membership_Test
                ( s1 : in QuadDobl_Sample_List; 
                  x : in QuadDobl_Complex_Vectors.Vector;
                  tol : in double_float; isin : out natural32;
                  s2,s2_last : in out QuadDobl_Sample_List );

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

  function Is_In ( s : QuadDobl_Sample_List; tol : double_float;
                   x : QuadDobl_Complex_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the element in s that matches with x.
  --   If x does not belong to s, then 0 is returned.

  function Is_In ( s : QuadDobl_Sample_List; tol : double_float;
                   spt : QuadDobl_Sample ) return natural32;

  -- DESCRIPTION :
  --   Returns the index of the element in s that matches with spt,
  --   within the given tolerance.  Returns 0 if spt does not belong to s.

  function Sample_Points ( s : QuadDobl_Sample_List )
                         return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the list of sample points, represented as solutions.

  function Select_Sub_List ( l : QuadDobl_Sample_List;
                             f : Standard_Natural_Vectors.Vector )
                           return QuadDobl_Sample_List;

  -- DESCRIPTION :
  --   Selects from the list only those samples whose position
  --   occurs in the vector f of labeled witness points.

  -- REQUIRED :
  --   The list must have sufficiently many samples and
  --   the points in f are sorted in increasing order.

  function Select_Sub_Grid ( grid : Array_of_QuadDobl_Sample_Lists;
                             f : Standard_Natural_Vectors.Vector )
                           return Array_of_QuadDobl_Sample_Lists;

  -- DESCRIPTION :
  --   Selects from the grid only those samples from the lists whose
  --   position occurs in the vector f of labeled witness points.

  -- REQUIRED :
  --   The lists in the grid have sufficiently many samples
  --   and the points in f are sorted in increasing order.

-- DESTRUCTORS :

  procedure Shallow_Clear ( samples : in out QuadDobl_Sample_List );
  procedure Shallow_Clear ( samples : in out Array_of_QuadDobl_Sample_Lists );

  -- DESCRIPTION :
  --   Deallocates only the encapsulating list structure.

  procedure Deep_Clear ( samples : in out QuadDobl_Sample_List );
  procedure Deep_Clear ( samples : in out Array_of_QuadDobl_Sample_Lists );

  -- DESCRIPTION :
  --   Deallocates all occupied memory by the structures.

end QuadDobl_Sample_Lists;
