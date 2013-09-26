with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Sample_Point_Lists;                 use Sample_Point_Lists;

package Divided_Differences is

-- DESCRIPTION :
--   The Newton form of an interpolating polynomial in two variables
--   requires the interpolation points to lie on a grid in the plane.
--   Given a grid of samples on a planar curve, the operations of
--   this package create, evaluate and expand the Newton form of the
--   interpolating polynomial through those samples.
--   The main intent of this package is to shadow the multi-precision
--   calculations with standard arithmetic for debugging purposes.

  type Newton_Interpolator is private;

-- CREATORS :

  function Create ( file : file_type; size : natural32;
                    stgrid : Array_of_Standard_Sample_Lists;
                    mpgrid : Array_of_Multprec_Sample_Lists )
                  return Newton_Interpolator;
  function Create ( file : file_type; size : natural32;
                    stgrid : Array_of_Standard_Sample_Lists;
                    mpgrid : Array_of_Multprec_Sample_Lists;
                    c : Standard_Complex_Numbers.Complex_Number )
                  return Newton_Interpolator;

  -- DESCRIPTION :
  --   Creates the Newton form of the interpolating polynomial through
  --   the samples in the grid.  If c is not supplied as parameter,
  --   a random number for c is generated.  In the interpolation, c is
  --   used for normalization:  p(x(0),c) = 1 where x(0) is the first 
  --   coordinate of the rotated grid of samples (p is the interpolator).

  -- REQUIRED :
  --   1) All samples in grid(i) lie in hyperplanes parallel to the
  --      hyperplanes used for the samples in grid(j);
  --   2) grid'range = 0..d, where d is the degree of the interpolator;
  --   3) Length_Of(grid(i)) = d, for all i in grid'range.

-- EVALUATOR :

  procedure Eval ( file : in file_type; q : in Newton_Interpolator;
                   stx : in Standard_Complex_Vectors.Vector;
                   mpx : in Multprec_Complex_Vectors.Vector;
                   steva : out Standard_Complex_Numbers.Complex_Number;
                   mpeva : out Multprec_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Evaluates the Newton interpolator both with standard and
  --   multi-precision arithmetic and writes diagnostics on file.

-- DESTRUCTOR :

  procedure Clear ( q : in out Newton_Interpolator );

  -- DESCRIPTION :
  --   Releases all memory allocated to store the Newton form
  --   which ceases to exist after this Clear operation.

private

  type Newton_Interpolator_Rep;
  type Newton_Interpolator is access Newton_Interpolator_Rep;

end Divided_Differences;
