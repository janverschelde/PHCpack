with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Matrices;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Irreducible_Decompositions;         use Irreducible_Decompositions;

package Witness_Generate_and_Classify is

-- DESCRIPTION :
--   This package offers three routines to implement the witness generate
--   and classify algorithm :
--    1) Witness_Generate runs down the homotopy embedding sequence and
--       finds lists of solutions that contain generic points and junk.
--    2) Witness_Classify classifies the solution lists found by 
--       Witness_Generate into generic points, grouped per irreducible
--       component, and junk points on higher dimensional components.
--    3) Witness_Generate_Classify interlaces Generate and Classify.

  procedure Witness_Generate
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32; degroco : in boolean;
                 zerotol : in double_float;
                 gentims : out Array_of_Duration;
                 flowtab : out Standard_Natural_Matrices.Matrix;
                 dc : out Standard_Irreducible_Decomposition );

  -- DESCRIPTION :
  --   Builds a sequence of homotopies between embedded systems to
  --   find generic points on each i-dimensional solution set of p,
  --   for i ranging from 0 to k.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   p         original system, must be square;
  --   k         top dimension;
  --   degroco   true if only root counting based on degrees is wanted,
  --             otherwise also polyhedral methods will be used as well;
  --   zerotol   tolerance to decide whether number is zero or not.

  -- ON RETURN :
  --   gentims   array of range 0..k, timings(i) lists time it took
  --             to solve the i-th system in the embedding;
  --   flowtab   matrix of range 0..k, 1..4 where the i-th row lists
  --               1) #paths traced at level i;
  --               2) #nonsolutions = solutions with zz /= 0, for i>0;
  --               3) #solutions with zz = 0, lie on components;
  --               4) #diverging paths;
  --   dc        contains embedded systems and solution lists.

  procedure Witness_Classify
               ( file : in file_type; full_output : in boolean;
                 dc : in out Standard_Irreducible_Decomposition;
                 method : in natural32; stoptol,membtol : in double_float;
                 clatims : out Array_of_Duration;
                 fp,fp_last : in out List );

  procedure Witness_Classify
               ( file : in file_type; full_output : in boolean;
                 dc : in out Multprec_Irreducible_Decomposition;
                 method,size : in natural32; stoptol,membtol : in double_float;
                 clatims : out Array_of_Duration;
                 fp,fp_last : in out List );

  -- DESCRIPTION :
  --   Classifies the solutions stored in dc into junk and generic
  --   points grouped according to the irreducible components.

  -- REQUIRED : there ought to be no junk at the top dimension.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   full_output indicates whether all sampling diagnostics (when true)
  --             or whether only summary (when false) will be written;
  --   dc        output of Witness_Generate;
  --   method    interpolating method used in breaking up
  --              = 0 : massive interpolate, comparison only at end,
  --              = 1 : incremental interpolate with linear projections,
  --              = 2 : interpolate with exploitation of span of component,
  --              = 3 : interpolate with central projections;
  --              = 4 : breakup with monodromy group actions;
  --   size      size of the multi-precision numbers;
  --   stoptol   tolerance to decide to stop interpolating,
  --             when method = 4, integer(stoptol) = stabilizing threshold;
  --   membtol   tolerance to determine membership.

  -- ON RETURN :
  --   dc        contains interpolation filters;
  --   clatims   stores timings at each level, for i in 0..Top_Dimension(dc);
  --   fp        list of characteristic data for each component, with
  --             its degree, dimension, and number of junk points on it;
  --   fp_last   pointer to last element of the list fp.

  procedure Witness_Generate_Classify
               ( file : in file_type; full_output : in boolean;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32; method : in natural32;
                 degroco : in boolean;
                 zerotol,stoptol,membtol : in double_float;
                 gentims,clatims : out Array_of_Duration;
                 flowtab : out Standard_Natural_Matrices.Matrix;
                 dc : out Standard_Irreducible_Decomposition;
                 fp,fp_last : in out List );

  procedure Witness_Generate_Classify
               ( file : in file_type; full_output : in boolean;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mp : in Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                 k : in integer32; method,size : in natural32;
                 degroco : in boolean;
                 zerotol,stoptol,membtol : in double_float;
                 gentims,clatims : out Array_of_Duration;
                 flowtab : out Standard_Natural_Matrices.Matrix;
                 dc : out Multprec_Irreducible_Decomposition;
                 fp,fp_last : in out List );

  -- DESCRIPTION :
  --   Interleaves Witness_Generate and Witness_Classify.
  --   The parameter mp is the original system with multi-precision
  --   coefficients.  For the other parameters, see the above routines.

end Witness_Generate_and_Classify;
