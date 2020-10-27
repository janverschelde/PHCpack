with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_VecVecs;
with Standard_Floating_VecVecs;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Main_Lifting_Functions is

-- DESCRIPTION :
--   Determines the lifting function depending on user selection.

  function Read_Integer_Lifting
              ( L : Lists_of_Integer_Vectors.List )
              return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION : interactive point-wise lifting of a list.
  --   Returns the list of lifted points.

  function Read_Float_Lifting
              ( L : Lists_of_Floating_Vectors.List )
              return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION : interactive point-wise lifting of a list.
  --   Returns the list of lifted points.

  procedure Integer_Lifting
              ( file : in file_type; p : in Poly_Sys;
                points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lilifu : in out Standard_Integer_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Displays the menu with the available lifting functions and
  --   performs the selected integer lifting function.

  -- NOTE : it is assumed that different supports are submitted.

  -- ON ENTRY :
  --   file     file that must be opened for output;
  --   p        polynomial system;
  --   points   supports of the system p.

  -- ON RETURN :
  --   lifted   the lifted support sets;
  --   lilifu   vectors used for linear lifting, otherwise lilifu = null.

  procedure Floating_Lifting
              ( file : in file_type; p : in Poly_Sys; b : in double_float;
                points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lilifu : in out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Displays the menu with the available lifting functions and
  --   performs the selected floating-point lifting function.

  -- NOTE : it is assumed that different supports are submitted.

  -- ON ENTRY :
  --   file     file that must be opened for output;
  --   p        polynomial system;
  --   b        lifting bound for stable mixed volumes;
  --   points   supports of the system p.

  -- ON RETURN :
  --   lifted   the lifted support sets;
  --   lilifu   vectors used for linear lifting, otherwise lilifu = null.

  procedure Main_Polynomial
              ( file : in file_type; p : in Poly_Sys;
                ipoints : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                fltlif : out boolean; stlb : out double_float;
                fpts : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ilftd : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flftd : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ililifu : in out Standard_Integer_VecVecs.Link_to_VecVec;
                flilifu : in out Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Main_Laurent
              ( file : in file_type; p : in Laur_Sys;
                ipoints : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                fltlif : out boolean; stlb : out double_float;
                fpts : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ilftd : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flftd : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                ililifu : in out Standard_Integer_VecVecs.Link_to_VecVec;
                flilifu : in out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   The user has the choice for integer or floating-point lifting.
  --   On return, output parameter fltlif is true if the user wants
  --   floating-point lifting and false otherwise.
  --   Depending on fltlif, the appropriate parameters are determined.

  -- ON ENTRY :
  --   file     file must be opened for output;
  --   p        a (Laurent) polynomial system;
  --   ipoints  the supports of the polynomial system.

  -- ON RETURN :
  --   fltlif   true if floating-point valued lifting was used;
  --   stlb     equals 0.0 if no lifting for stable mixed volumes,
  --            otherwise equals the bound for stable mixed volumes.

end Main_Lifting_Functions;
