with text_io;                            use text_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Exponent_Vectors;                   use Exponent_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Integer_Polyhedral_Continuation is

-- DESCRIPTION :
--   This package implements polyhedral homotopy continuation methods,
--   based on mixed subdivision induced by integer-valued lifting.
--   The continuation is organized in three layers:
--     1. inner normal, tracking of paths for one poly
--     2. mixed cell, recursion is needed when the cell is not fine;
--     3. mixed subdivision, for all cells in the subdivision.
--   Each layer has four versions: all combinations of
--     1. homotopy as polynomial system or coefficient homotopy;
--     2. silent and reporting version.

-- FIRST LAYER : polyhedral continuation for one transformation.

  procedure Mixed_Continuation
                ( p : in Laur_Sys; normal : in Vector; 
                  sols : in out Solution_List );

  procedure Mixed_Continuation
                ( file : in file_type; p : in Laur_Sys;
                  normal : in Vector; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Polyhedral continuation for transformed polyhedral homotopy.

  -- ON ENTRY :
  --   file       file to write intermediate results on;
  --   p          a transformed lifted Laurent polynomial system;
  --   normal     normal to a mixed cell;
  --   sols       start solutions of the subsystem which corresponds
  --              with the mixed cell with given inner normal.

  -- ON RETURN :
  --   sols       the solutions of p, which correspond to one mixed cell.

  procedure Mixed_Continuation
                ( mix : in Standard_Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  normal : in Vector; sols : in out Solution_List );

  procedure Mixed_Continuation
                ( file : in file_type; mix : in Standard_Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  normal : in Vector; sols : in out Solution_List );

  -- DESCRIPTION : polyhedral continuation with coefficient homotopy.

  -- ON ENTRY :
  --   file       file to write intermediate results on;
  --   mix        type of mixture;
  --   lifted     lifted supports of polynomial system, in original order;
  --   h          coefficient homotopy;
  --   c          coefficients of homotopy;
  --   e          the exponent vectors of the unlifted system;
  --   j          coefficient Jacobian matrix;
  --   m          multiplication factors in coefficient Jacobian matrix;
  --   normal     normal to a mixed cell;
  --   sols       start solutions of the subsystem which corresponds
  --              with the mixed cell with given inner normal.

  -- ON RETURN :
  --   sols       the solutions of p, which correspond to one mixed cell.

-- SECOND LAYER : polyhedral continuaton for one mixed cell.

  procedure Mixed_Solve
                ( p : in Laur_Sys; mix : in Vector; 
                  mic : in Mixed_Cell; sols,sols_last : in out Solution_List );

  procedure Mixed_Solve
                ( file : in file_type; p : in Laur_Sys; mix : in Vector;
                  mic : in Mixed_Cell; sols,sols_last : in out Solution_List );

  -- DESCRIPTION :
  --   Polyhedral continuation for one mixed cell.

  -- REQUIRED : polynomials in p must be ordered according to mix.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a lifted Laurent polynomial system;
  --   mix        type of mixture;
  --   mic        a mixed cell;

  -- ON RETURN :
  --   sols       the solution list of p;
  --   sols_last  pointer to last element of the list sols.

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List );

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List );

  -- DESCRIPTION :
  --   Polyhedral coefficient-homotopy continuation for one mixed cell.

  -- REQUIRED : polynomials in p must be ordered according to mix.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a lifted Laurent polynomial system;
  --   lifted     lifted supports, with original order of points;
  --   h          coefficient homotopy;
  --   c          coefficients of homotopy;
  --   e          the exponent vectors of the unlifted system;
  --   j          coefficient Jacobian matrix;
  --   m          multiplication factors in coefficient Jacobian matrix;
  --   mix        type of mixture;
  --   mic        a mixed cell.

  -- ON RETURN :
  --   sols       the solution list of p;
  --   sols_last  pointer to last element of the list sols.

-- THIRD LAYER : polyhedral continuation for a mixed subdivision.

  procedure Mixed_Solve
                ( p : in Laur_Sys;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List );

  procedure Mixed_Solve
                ( file : in file_type; p : in Laur_Sys;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List );

  -- DESCRIPTION :
  --   Polyhedral coefficient-homotopy continuation for a mixed subdivision.

  -- REQUIRED : polynomials in p must be ordered according to mix.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a lifted Laurent polynomial system;
  --   mix        type of mixture;
  --   mixsub     a collection of mixed cells.

  -- ON RETURN :
  --   sols       the solution list of p.

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List );

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Standard_Complex_VecVecs.VecVec;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jaco_Mat; m : in Mult_Factors;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List );

  -- DESCRIPTION :
  --   Polyhedral coefficient-homotopy continuation for a mixed subdivision.

  -- REQUIRED : polynomials in p must be ordered according to mix.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a lifted Laurent polynomial system;
  --   lifted     lifted supports, in original order;
  --   h          coefficient homotopy;
  --   c          coefficients of homotopy;
  --   e          the exponent vectors of the unlifted system;
  --   j          coefficient Jacobian matrix;
  --   m          multiplication factors in coefficient Jacobian matrix;
  --   mix        type of mixture;
  --   mixsub     a collection of mixed cells.

  -- ON RETURN :
  --   sols       the solution list of p.

end Integer_Polyhedral_Continuation;
