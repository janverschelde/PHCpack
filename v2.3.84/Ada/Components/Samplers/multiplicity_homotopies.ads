with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;

package Multiplicity_Homotopies is

-- DESCRIPTION :
--   This package provides homotopies to locally provide samples of
--   positive dimensional solution components of multiplicity > 1.
--   The basic data are an embedded polynomial system and a list of
--   clustered points, with the multiplicity as the length of the list.
--   The embedding contains two special types of moving equations:
--    1) the regular slices to get to generic points; and
--    2) the equation zz = z(eps), with slack variable zz,
--       to control the distance to the component.

  procedure Reconditioning_Homotopy
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 q : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 k : in integer32;
                 a : in Standard_Complex_Numbers.Complex_Number;
                 ind : in integer32;
                 ze : out Standard_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Does reconditioning of the initial approximations in the cluster.

  -- ON ENTRY :
  --   file      to write intermediate results and diagnostics;
  --   p         embedding for a space curve;
  --   q         start system used to solve p;
  --   sols      list of solutions leading up to the cluster;
  --   k         relaxation parameter used in homotopy to solve p;
  --   a         accessibility constant used in the homotopy to solve p;
  --   ind       index to the equation p(ind) : zz = 0.

  -- ON RETURN :
  --   p         ze is substracted from p(ind);
  --   sols      solutions for the value zz = ze;
  --   ze        random complex value with modulus 10^(-16/m),
  --             where m is the multiplicity = Length_Of(sols).

  procedure Incoming_Homotopy
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 ind : in integer32;
                 ze : in Standard_Complex_Numbers.Complex_Number;
                 mpsols : out Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Improves the estimates of the clusters towards the multiple point.

  -- ON ENTRY :
  --   file      to write intermediate results and diagnostics;
  --   p         embedding for space curve;
  --   sols      list of clustered solutions, solutions for zz = zze;
  --   ind       index to the equation p(ind) : zz - ze = 0;
  --   ze        value that will become (close to) zero.

  -- ON RETURN :
  --   mpsols    refined solutions for the slack variable zz closer to zero,
  --             the distance between the solutions in the cluster should
  --             have become much smaller.

  procedure Sampling_Homotopy
               ( file : in file_type; eps : in double_float;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 hyp : out Standard_Complex_VecVecs.VecVec;
                 sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Perturbs the slice, adding eps to real and imaginary part of the
  --   coefficients of the slice and refines the solutions.

  -- REQUIRED : the slice is the last equation of the polynomial system.

  -- ON ENTRY :
  --   file      to write intermediate results and diagnostics;
  --   p         embedding for space curve;
  --   sols      solutions to p.

  -- ON RETURN :
  --   p         system with modified last equation;
  --   hyp       new hyperplane sections to cut out generic points;
  --   sols      refined solutions to p.

end Multiplicity_Homotopies;
