with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Integer_Matrices;
with Standard_Complex_Matrices;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Exponent_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Multitasking_Polyhedral_Trackers is

-- DESCRIPTION :
--   This package offers routines to implement polyhedral continuation
--   methods for multithreading (shared memory parallel computing).

  procedure Silent_Multitasking_Path_Tracker
              ( q : in Laur_Sys; nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision; h : in Eval_Coeff_Laur_Sys;
                c : in Standard_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in Eval_Coeff_Jaco_Mat; mf : in Mult_Factors;
                sols : out Solution_List );
  procedure Reporting_Multitasking_Path_Tracker
              ( q : in Laur_Sys; nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision; h : in Eval_Coeff_Laur_Sys;
                c : in Standard_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in Eval_Coeff_Jaco_Mat; mf : in Mult_Factors;
                sols : out Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q.
  --   The reporting version allows the user to monitor the progress 
  --   of the computations on screen.

  -- ON ENTRY :
  --   q        a random coefficient system;
  --   nt       number of tasks;
  --   n        ambient dimension;
  --   r        number of different supports;
  --   mix      type of mixture;
  --   lif      supports of q, with a lifting as occurred in mcc;
  --   mcc      a regular mixed-cell configuration;
  --   h        coefficient homotopy;
  --   c        coefficients of the homotopy h;
  --   j        coefficient Jacobian matrix;
  --   mf       multiplication factors in coefficient Jacobian matrix.
 
  -- ON RETURN :
  --   sols     all solutions of the system q.

  procedure Silent_Multitasking_Path_Tracker
              ( q : in Laur_Sys; nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision; sols : out Solution_List );
  procedure Reporting_Multitasking_Path_Tracker
              ( file : in file_type; q : in Laur_Sys; nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision; sols : out Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q.
  --   The reporting version allows monitoring the computations.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   q        a random coefficient system;
  --   nt       number of tasks;
  --   n        ambient dimension, before the lifting;
  --   m        number of different supports;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.
 
  -- ON RETURN :
  --   sols     all solutions of the system q.

end Multitasking_Polyhedral_Trackers;
