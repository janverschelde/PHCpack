with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Homotopy_Membership_Tests is

-- DESCRIPTION :
--   This package offers an implementation of the homotopy test to
--   solve the membership problem for a point to be on a component.

  procedure Homotopy_Membership_Test
                 ( file : in file_type;
                   genpts : in Solution_List; x : in Vector;
                   adjsli : in Standard_Complex_VecVecs.VecVec;
                   homtol : in double_float; found : out boolean );

  -- DESCRIPTION :
  --   Implementation of the homotopy membership test procedure.
  --   We invoke the Sample_with_Stop routine from the Sampling Machine.
  --   Since this routine stops when the given vector x is found, 
  --   we compare x against the last vector in the list of new solutions.

  procedure Homotopy_Membership_Test
                 ( file : in file_type;
                   ep : in Poly_Sys; dim : in natural32;
                   sli : in Standard_Complex_VecVecs.VecVec;
                   genpts : in Solution_List; s : in Solution;
                   restol,homtol : in double_float;
                   success,found : out boolean );

  -- DESCRIPTION :
  --   This procedure performs the component test for solution s.

  -- ON ENTRY :
  --   file        for intermediate output and diagnostics;
  --   ep          embedded polynomial system;
  --   dim         dimension of solution components;
  --   sli         slices that cut out the generic points;
  --   genpts      generic points on the components;
  --   s           point up for testing;
  --   restol      tolerance for residual;
  --   homtol      tolerance for homotopy test.

  -- ON RETURN :
  --   success     true when test point satisfies the equations,
  --               false otherwise;
  --   found       true when homotopy test succeeds in finding a matching
  --               generic point, found is false when not success.

  procedure Homotopy_Membership_Test
                 ( file : in file_type;
                   ep : in Poly_Sys; dim : in natural32;
                   genpts,sols : in Solution_List;
                   restol,homtol : in double_float );

  procedure Homotopy_Membership_Test
                 ( file : in file_type;
                   ep : in Poly_Sys; dim : in natural32;
                   genpts,sols : in Solution_List;
                   restol,homtol : in double_float;
                   isols,isols_last : in out Solution_List );

  -- DESCRIPTION :
  --   This procedure performs the component test for the solutions sols.
  --   If provided, the solutions that are not found to be components are
  --   appended to the list isols.

  -- ON ENTRY :
  --   file        for intermediate output and diagnostics;
  --   ep          embedded polynomial system;
  --   dim         number of slices;
  --   genpts      generic points on the slices;
  --   sols        points up for testing, do they lie on the components?
  --   restol      tolerance for residual;
  --   homtol      tolerance for homotopy test.

  -- ON RETURN :
  --   isols       solutions found to be isolated;
  --   isols_last  head of this list is the last element in isols.

end Homotopy_Membership_Tests;
