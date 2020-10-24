with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;

package Test_Product_Systems is

-- DESCRIPTION :
--   Development of tools to work with products of polynomials.

  procedure Write_Support ( h : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   The support of h marks the nonzero coefficients in h.

  function Same_Support
             ( h1,h2 : Standard_Complex_Vectors.Vector)
             return boolean;

  -- DESCRIPTION :
  --   Returns true if both h1 and h2 share the same support.

  -- REQUIRED : h1'range = h2'range.

  function Hyper_Multiplicities return Standard_Natural_VecVecs.VecVec;

  -- DESCRIPTION :
  --   The multiplicities of the hyperplanes (or hyper multiplicities) are:
  --   (1) equal to the number of times it occurs with the same support
  --       for the first occurrence of the hyperplane with same support;
  --   (2) zero for all the subsequent occurrences with same support.

  -- REQUIRED :
  --   Standard_Linear_Product_System is initialized properly.

  procedure Write_Multiplicities ( m : in Standard_Natural_VecVecs.VecVec);

  -- DESCRIPTION :
  --   Writes the multiplicities of the linear-product start system.

  procedure Test_Read_and_Write;

  -- DESCRIPTION :
  --   Tests reading and writing of a linear-product system.

  procedure Test_Product_Systems;

  -- DESCRIPTION :
  --   Interactive test on a system of product polynomials.

  function Count_All_Solutions_of_Linear_Product_System
             ( output : in boolean ) return natural32;

  -- DESCRIPTION :
  --   Uses the enumerator in Standard_Linear_Product_System
  --   to count all solutions of a linear-product system.

  procedure Count_All_Solutions_of_Linear_Product_System;

  -- DESCRIPTION :
  --   Applies the count of all solutions to the stored
  --   linear-product system.

  procedure Set_Structure_Linear_Product_Start_System;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and creates a linear-product
  --   start system based on a generated set structure.

  function Residual ( A : Standard_Complex_Matrices.Matrix;
                      b,x : Standard_Complex_Vectors.Vector )
                    return double_float;

  -- DESCRIPTION :
  --   Returns the max norm of b - A*x.

  function Lex_Count_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector ) return natural32;

  -- DESCRIPTON :
  --   Enumerates lexicographically all possible combinations of
  --   the degrees in d.  Solves a linear system for each combination.
  --   There is no intermediate output.

  function Multiplicity
              ( m : Standard_Natural_VecVecs.VecVec;
                a : Standard_Natural_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the multiplicity of the selection a,
  --   using the hyper multiplicities in m.

  function Lex_Count_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector;
                m : Standard_Natural_VecVecs.VecVec ) return natural32;

  -- DESCRIPTON :
  --   Enumerates lexicographically all possible combinations of
  --   the degrees in d.  Solves a linear system for each combination.
  --   To make it go faster, the hyper multiplicities are used.
  --   There is no intermediate output.

  procedure Lex_Enumerate_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Lexicographical enumeration of all solutions
  --   of a linear-product system.

  procedure Enumerate_All_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Enumerates all solutions of a linear-product system.

  procedure Test_Solution_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   With in d the degrees of the polynomials in the start system,
  --   the solution of the system is tested.

  procedure Count_Roots ( q : in Prod_Sys );

  -- DESCRIPTION :
  --   Displays a menu and prompts to select a root counter.

  procedure Solve_Linear_Product_Start_System;

  -- DESCRIPTION :
  --   Prompts for a linear-product system and solves it.

  procedure Interactive_Inverted_Enumerator ( q : in Prod_Sys );

  -- DESCRIPTION :
  --   A step-by-step application of the inverted enumerator
  --   for the solutions of a linear-product system q.

  function Create ( v : Standard_Complex_Vectors.Vector;
                    rcond : double_float ) return Solution;

  -- DESCRIPTION :
  --   Returns the solution data with vector v and rco in rcond.

  procedure Solve_by_Inverted_Enumerator
               ( file : in file_type; q : in Prod_Sys; rc : in natural32 );

  -- DESCRIPTION :
  --   Writes the linear-product system q and its solutions to file.
  --   The root count is given by rc.

  procedure Run_Inverted_Enumerator
              ( file : in out file_type; q : in Prod_Sys );

  -- DESCRIPTION :
  --   Runs the inverted enumerator to count first all the solution of q
  --   and then appends to the file the solutions.

  procedure Test_Inverted_Enumerator;

  -- DESCRIPTION :
  --   Runs an interactive test on the inverted enumerator
  --   to solve a linear-product start system.

  procedure Write_Supports ( nf : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the supporting set structure.

  -- REQUIRED :
  --   Standard_Linear_Product_System has been initialized properly.

  procedure Support_of_Linear_Product_System;

  -- DESCRIPTION :
  --   Prompts for a linear-product system
  --   and then extracts the supporting set structure.

  function Difference ( x,y : Standard_Complex_Matrices.Matrix )
                      return double_float;

  -- DESCRIPTION :
  --   Returns the sum of all differences (in modulus) of all values
  --   in x and y.

  procedure Evaluate_Linear_Product_System;

  -- DESCRIPTION :
  --   Compares different ways to evaluate a linear-product system.

  -- REQUIRED :
  --   Standard_Linear_Product_System contains all coefficients
  --   for a linear-product start system.

  procedure Test_Evaluation;

  -- DESCRIPTION :
  --   Prompts for a linear-product system and tests the evaluation.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Product_Systems;
