with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Multitasked_DD_QD_Refiners is

-- DESCRIPTION :
--   Applies Newton's method in double double and quad double arithmetic
--   to refine roots with multitasking.

  procedure Refine_Solution
               ( id,nbt,solno : in integer32; output : in boolean;
                 ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
                 f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 epsxa,epsfa : in double_float; maxit : natural32 );
  procedure Refine_Solution
               ( id,nbt,solno : in integer32; output : in boolean;
                 ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
                 f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 epsxa,epsfa : in double_float; maxit : natural32 );

  -- DESCRIPTION :
  --   Task with identification number id reports the receipt of
  --   solution with number solno, with data in ls.

  procedure Multitasking_Refinement
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List; 
                 n : in integer32; output : in boolean;
                 epsxa,epsfa : in double_float; maxit : in natural32 );
  procedure Multitasking_Refinement
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List; 
                 n : in integer32; output : in boolean;
                 epsxa,epsfa : in double_float; maxit : in natural32 );

  -- DESCRIPTION :
  --   Given a polynomial system p with solutions in sols,
  --   n threads will be launched to multitask on the solution list.

  -- ON ENTRY :
  --   p         a polynomial system, as many equations as unknowns;
  --   sols      approximations for the solutions of p;
  --   n         number of tasks;
  --   output    flag if output is requested;
  --   epsxa     tolerance on the update for each solution;
  --   epsfa     tolerance on the residual;
  --   maxit     maximum number of iterations.

  -- ON RETURN :
  --   sols      the refined list of solutions.

end Multitasked_DD_QD_Refiners;
