with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Multitasking_Polyhedral_Starters is

-- DESCRIPTION :
--   This package provides procedures to solve all start systems with
--   multithreading (shared memory parallel computing) as used in
--   polyhedral homotopies.  These solvers are written in preparation 
--   for the actual polyhedral path trackers to verify whether all
--   start systems are solved correctly.

  procedure Silent_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );
  procedure Reporting_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );
  procedure Silent_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector );
  procedure Reporting_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector );
  procedure Silent_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector );
  procedure Reporting_Multithreaded_Solve_Start_Systems
              ( nt : in integer32;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the start systems of q defined by the mixed cells in mcc,
  --   using nt threads of execution.  The silent version does not
  --   write intermediate output to screen, which the reporting one does.

  -- REQUIRED : the order of the polynomials in the start system must
  --   agree with the order of the supports in the cells, and
  --   sols is an array of range 1..nt = res'range.

  -- ON ENTRY :
  --   nt       the number of tasks;
  --   q        a random coefficient start system;
  --   r        number of distinct supports;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN :
  --   mixvol   sum of the mixed volumes of the cells in mcc;
  --   sols     the start solutions computed by the tasks,
  --            in particular sols(i) is computed by task i,
  --            following a static workload assignment.
  --   res      res(i) contains the sums of all residuals of
  --            the start solutions computed by task i.

end Multitasking_Polyhedral_Starters;
