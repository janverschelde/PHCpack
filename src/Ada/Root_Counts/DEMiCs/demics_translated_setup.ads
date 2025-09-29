with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with DEMiCs_Input_Data;
with DEMiCs_MVC;

package DEMiCs_Translated_Setup is

-- DESCRIPTION :
--   The procedures in this package set the stage for applying DEMiCS:
--   (1) extracting the supports from a system,
--   (2) entering the supports into a dataSet,
--   (3) prompting for a user given lifting,
--   (4) adding an artificial origin for stable mixed volumes.

  procedure Make_Supports 
              ( data : in out DEMiCs_Input_Data.class_dataSet.dataSet;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given the supports in sup and type of mixture in mix,
  --   defines the data as a dataSet object.

  procedure Make_Data
              ( res : out DEMiCs_Input_Data.class_dataSet.dataSet;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                storemix : in boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns the data object for input to DEMiCs.
  --   If storemix, then the dimension and mixture type will be
  --   stored in the DEMiCs_Output_Cells for processing the labels
  --   to the cells into a mixed cell configuration.

  function Make_Data ( p : Poly_sys; storemix : boolean;
                       vrblvl : integer32 := 0 ) 
                     return DEMiCs_Input_Data.class_dataSet.dataSet;

  -- DESCRIPTION :
  --   Returns the data object for input to DEMiCs,
  --   for the polynomial system p.

  function Make_Data ( p : Laur_sys; storemix : boolean;
                       vrblvl : integer32 := 0 ) 
                     return DEMiCs_Input_Data.class_dataSet.dataSet;

  -- DESCRIPTION :
  --   Returns the data object for input to DEMiCs,
  --   for the Laurent polynomial system p.

  function Extract_Support
             ( dim : integer32;
               sup : Standard_Floating_Vectors.Link_to_Vector;
               vrblvl : integer32 := 0 )
             return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns the support as a list, as defined in sup.

  function Extract_Supports
             ( ptr2MVC : DEMiCs_MVC.class_mvc.Link_to_mvc;
               vrblvl : integer32 := 0 )
             return Arrays_of_Integer_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Retrieves the original supports stored in class_simplex,
  --   available as ptr2MVC.the_Simplex.

  function Apply_Lifting
             ( sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lft : Standard_Floating_Vectors.Link_to_Vector;
               vrblvl : integer32 := 0 )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Extends the supports sup with the lifting in lft.

end DEMiCs_Translated_Setup;
