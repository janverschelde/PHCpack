with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Multitasking_Sampling is

-- DESCRIPTION :
--   Offers routines to takes samples of positive dimensional solution sets,
--   using multiple tasks, as a stepping stone to implement factorization
--   with monodromy loops.

  function Change_Slices
               ( p : in Poly_Sys;
                 h : in Standard_Complex_VecVecs.VecVec ) return Poly_Sys;

  -- DESCRIPTION :
  --   Given an embedded system in p for a d-dimensional solution set,
  --   on return is a copy of p, except for the last d polynomials.
  --   The coefficients of the last d polynomials are given in h.

  -- REQUIRED : h'range = 1..d.

  procedure Silent_Sampler
               ( s : in out Solution_List; n,d : in integer32;
                 h : in Standard_Complex_VecVecs.VecVec;
                 gamma : in Complex_Number; q : in Poly_Sys;
                 p : out Poly_Sys );
  procedure Reporting_Sampler
               ( s : in out Solution_List; n,d : in integer32;
                 h : in Standard_Complex_VecVecs.VecVec;
                 gamma : in Complex_Number;
                 q : in Poly_Sys; p : out Poly_Sys );

  -- DESCRIPTION :
  --   Moves the witness points in s to the new hyperplanes in h.
  --   The reporting version allows the user to monitor the progress.

  -- ON ENTRY :
  --   s         list of embedded solutions of q;
  --   n         number of tasks;
  --   d         dimension of the solution set;
  --   h         new set of hyperplanes for the new witness set;
  --   gamma     random constant for the homotopy;
  --   q         embedded polynomial system for the witness set.

  -- ON RETURN :
  --   s         solutions of p, on the hyperplanes h;
  --   p         embedded system for the new witness set,
  ---            in particular, we have: p = Change_Slices(q,d,h).

  procedure Driver_to_Sampler
               ( file : in file_type; n,d : in integer32;
                 ep : in Poly_Sys; esols : in out Solution_List );

  -- DESCRIPTION :
  --    Interactive driver to use n tasks to sample a d-dimensional solution
  --    set given by an initial witness set with embedded system in ep and
  --    solutions in esols.  Generates a new random set of d hyperplanes.
  --    Writes the new witness set to file.

  -- ON ENTRY :
  --    file     output file to write the new witness set on;
  --    n        number of tasks;
  --    d        dimension of the solution set;
  --    ep       embedding of a polynomial system, with d slack variables,
  --             and d generic hyperplanes;
  --    esols    embedded list of solution for ep.

  procedure Driver_to_Monodromy
               ( file : in file_type; n,d : in integer32;
                 ep : in Poly_Sys; esols : in out Solution_List );

  -- DESCRIPTION :
  --   The witness set defined by (ep,esols) is the starting point
  --   for a number of monodromy loops.

end Multitasking_Sampling;
