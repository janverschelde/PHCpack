with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Partitions_of_Sets_Of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Black_Box_Root_Counters is

-- DESCRIPTION :
--   This package provides an interface to the various root counting
--   capabilities in PHCpack, called by phc -b.

  procedure Write_Root_Counts
               ( file : in file_type; no_mv : in boolean;
                 d : in natural64; mp_d : in Natural_Number;
                 m : in natural32; bz,bs : in natural64;
                 mv,smv : in natural32; z : in Partition );

  -- DESCRIPTION :
  --   Writes root counts and set structures to file.

  -- ON ENTRY :
  --   file      to write root counts on (could be standard_output);
  --   no_mv     if no mixed volume was computed;
  --   d         total degree;
  --   mp_d      multiprecision version of total degree (if overflow);
  --   m         number of sets in the m-homogeneous Bezout number;
  --   bz        m-homogeneous Bezout number;
  --   bs        set structure Bezout bound;
  --   mv        mixed volume;
  --   smv       stable mixed volume;
  --   z         partition of variables for m-homogeneous Bezout number.

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Solution_List;
                 rocotime,hocotime : out duration );

  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Solution_List;
                 rocotime,hocotime : out duration );

  -- DESCRIPTION :
  --   Calculates four different root counts: total degree, m-homogeneous
  --   Bezout number, generalized Bezout number based on set structure,
  --   and mixed volume.  Heuristics are used for the Bezout numbers.
  --   Returns the start system with lowest root count and least amount
  --   of work, which means that linear-product start systems are preferred,
  --   when Bezout numbers equal the mixed volume.
  --   In case the mixed volume is zero because of a equation with only
  --   one monomial, the smallest Bezout bound is returned.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   silent    if silent, then the computed root counts will not be shown
  --             on screen, otherwise, the root counter remains silent;
  --   p         a polynomial system.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   deg       if true, then only degree bounds are used,
  --             otherwise also the mixed volume is computed;
  --   rc        root count, Bezout number or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Solution_List;
                 rocotime,hocotime : out duration );

  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Solution_List;
                 rocotime,hocotime : out duration );

  -- DESCRIPTION :
  --   Wrapper to a mixed-volume computation for a Laurent polynomial system.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   silent    if silent, then the mixed volume will not be written
  --             to screen, otherwise, if silent, then the user will
  --             be shown the mixed volume on screen;
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   p         a Laurent polynomial system.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        the root count is the mixed volume, equals Lenght_Of(qsols);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

end Black_Box_Root_Counters;
