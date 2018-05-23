with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_JacoMats;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with Exponent_Vectors;
with Semaphore;

package Pipelined_Cell_Trackers is

-- DESCRIPTION :
--   The granularity of the path trackers is such that one task is reponsible
--   for tracking all solution paths defined by one mixed cell in a regular
--   mixed cell configuration for a polyhedral homotopy.

  procedure Standard_Track_Cell
              ( sem : in out Semaphore.Lock; idtask,n,r : in integer32; 
                mix : in Standard_Integer_Vectors.Vector;
                mic : in Mixed_Cell;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : in Standard_Complex_VecVecs.VecVec;
                dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                cft : in Standard_Complex_VecVecs.Link_to_VecVec;
                epv : in Exponent_Vectors.Exponent_Vectors_Array;
                hom : in Standard_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                ejf : in Standard_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                jmf : in Standard_Complex_Laur_JacoMats.Mult_Factors;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                tmv : in out natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                last : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   One task tracks all paths defined by a polyhedral homotopy determined
  --   by one mixed cell in a regular mixed cell configuration.
  --   As many paths as the volume of the cell are tracked,
  --   in standard double precision.

  -- ON ENTRY :
  --   sem      semaphore for the critical sections;
  --   idtask   identification number of the task;
  --   n        ambient dimension of the points, before lifting;
  --   r        number of different supports;
  --   mix      type of mixture, counts the number of occurrences;
  --   mic      a mixed cell in a regular mixed cell configuration;
  --   lif      lifted supports, in an array of range 1..r;
  --   cff      coefficients of the system q;
  --   dpw      pointer to work space for the powers in the homotopy;
  --   cft      pointer to work space for coefficient in t;
  --   epv      exponents of the homotopy;
  --   hom      evaluable form of the coefficient polyhedral homotopy;
  --   ejf      evaluable form of the Jacobian matrix;
  --   jmf      multiplication factors of the Jacobian matrix;
  --   q        a random coefficient system;
  --   tmv      current number of paths tracked by the task;
  --   sols     pointer to the first solution in a list;
  --   last     pointer to the last solution in a list.

  -- ON RETURN :
  --   tmv      updated number of paths tracked by the task;
  --   sols     if sols was null on entry, then sols points to the first
  --            solution contributed by this first cell;
  --   last     pointer has moved as many positions as the mixed volume
  --            of the cell and the list sols contains as many additional
  --            solutions of q as that mixed volume.

  procedure DoblDobl_Track_Cell
              ( sem : in out Semaphore.Lock; idtask,n,r : in integer32; 
                mix : in Standard_Integer_Vectors.Vector;
                mic : in Mixed_Cell;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : in DoblDobl_Complex_VecVecs.VecVec;
                dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                cft : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                epv : in Exponent_Vectors.Exponent_Vectors_Array;
                hom : in DoblDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                ejf : in DoblDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                jmf : in DoblDobl_Complex_Laur_JacoMats.Mult_Factors;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                tmv : in out natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                last : in out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   One task tracks all paths defined by a polyhedral homotopy determined
  --   by one mixed cell in a regular mixed cell configuration.
  --   As many paths as the volume of the cell are tracked.

  -- ON ENTRY :
  --   sem      semaphore for the critical sections;
  --   idtask   identification number of the task;
  --   n        ambient dimension of the points, before lifting;
  --   r        number of different supports;
  --   mix      type of mixture, counts the number of occurrences;
  --   mic      a mixed cell in a regular mixed cell configuration;
  --   lif      lifted supports, in an array of range 1..r;
  --   cff      coefficients of the system q;
  --   dpw      pointer to work space for the powers in the homotopy;
  --   cft      pointer to work space for coefficient in t;
  --   epv      exponents of the homotopy;
  --   hom      evaluable form of the coefficient polyhedral homotopy;
  --   ejf      evaluable form of the Jacobian matrix;
  --   jmf      multiplication factors of the Jacobian matrix;
  --   q        a random coefficient system;
  --   tmv      current number of paths tracked by the task;
  --   sols     pointer to the first solution in a list;
  --   last     pointer to the last solution in a list.

  -- ON RETURN :
  --   tmv      updated number of paths tracked by the task;
  --   sols     if sols was null on entry, then sols points to the first
  --            solution contributed by this first cell;
  --   last     pointer has moved as many positions as the mixed volume
  --            of the cell and the list sols contains as many additional
  --            solutions of q as that mixed volume.

  procedure QuadDobl_Track_Cell
              ( sem : in out Semaphore.Lock; idtask,n,r : in integer32; 
                mix : in Standard_Integer_Vectors.Vector;
                mic : in Mixed_Cell;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : in QuadDobl_Complex_VecVecs.VecVec;
                dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                cft : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                epv : in Exponent_Vectors.Exponent_Vectors_Array;
                hom : in QuadDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                ejf : in QuadDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                jmf : in QuadDobl_Complex_Laur_JacoMats.Mult_Factors;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                tmv : in out natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                last : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   One task tracks all paths defined by a polyhedral homotopy determined
  --   by one mixed cell in a regular mixed cell configuration.
  --   As many paths as the volume of the cell are tracked.

  -- ON ENTRY :
  --   sem      semaphore for the critical sections;
  --   idtask   identification number of the task;
  --   n        ambient dimension of the points, before lifting;
  --   r        number of different supports;
  --   mix      type of mixture, counts the number of occurrences;
  --   mic      a mixed cell in a regular mixed cell configuration;
  --   lif      lifted supports, in an array of range 1..r;
  --   cff      coefficients of the system q;
  --   dpw      pointer to work space for the powers in the homotopy;
  --   cft      pointer to work space for coefficient in t;
  --   epv      exponents of the homotopy;
  --   hom      evaluable form of the coefficient polyhedral homotopy;
  --   ejf      evaluable form of the Jacobian matrix;
  --   jmf      multiplication factors of the Jacobian matrix;
  --   q        a random coefficient system;
  --   tmv      current number of paths tracked by the task;
  --   sols     pointer to the first solution in a list;
  --   last     pointer to the last solution in a list.

  -- ON RETURN :
  --   tmv      updated number of paths tracked by the task;
  --   sols     if sols was null on entry, then sols points to the first
  --            solution contributed by this first cell;
  --   last     pointer has moved as many positions as the mixed volume
  --            of the cell and the list sols contains as many additional
  --            solutions of q as that mixed volume.

end Pipelined_Cell_Trackers;
