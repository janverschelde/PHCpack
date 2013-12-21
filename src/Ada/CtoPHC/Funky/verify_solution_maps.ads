with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Complex_Polynomial_Matrices;        use Complex_Polynomial_Matrices;

procedure Verify_Solution_Maps
            ( file : in file_type; pts : in Vector; planes : in VecMat;
              solmaps : in Array_of_Polynomial_Matrices;
              tol : in double_float; fail : out boolean );

-- DESCRIPTION :
--   Verifies whether the solution maps intersect the planes when
--   evaluated at the points.

-- ON ENTRY :
--   file     for intermediate diagnostics;
--   pts      interpolation points where to evaluate the maps;
--   planes   the input m-planes the solution maps have to meet;
--   solmaps  solution maps producing p-planes of degree q;
--   tol      tolerance to decide whether number is zero.

-- ON RETURN :
--   fail     true if not all determinants are less than tolerance. 
