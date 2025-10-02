with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --  Makes the data object for input to DEMiCs,
  --  given the supports sup, which may be permuted when the type
  --  of mixture is computed.

  function Make_Data ( p : Poly_sys; storemix : boolean;
                       vrblvl : integer32 := 0 ) 
                     return DEMiCs_Input_Data.class_dataSet.dataSet;

  -- DESCRIPTION :
  --   Returns the data object for input to DEMiCs,
  --   for the polynomial system p.
  --   If storemix, then the dimension and mixture type will be
  --   stored in the DEMiCs_Output_Cells for processing the labels
  --   to the cells into a mixed cell configuration.

  function Make_Data ( p : Laur_sys; storemix : boolean;
                       vrblvl : integer32 := 0 ) 
                     return DEMiCs_Input_Data.class_dataSet.dataSet;

  -- DESCRIPTION :
  --   Returns the data object for input to DEMiCs,
  --   for the Laurent polynomial system p.
  --   If storemix, then the dimension and mixture type will be
  --   stored in the DEMiCs_Output_Cells for processing the labels
  --   to the cells into a mixed cell configuration.

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
  --   Extends the supports sup with the lifting in lft,
  --   assuming every support set appears only once in sup.

  function Apply_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lft : Standard_Floating_Vectors.Link_to_Vector;
               vrblvl : integer32 := 0 )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Extends the supports sup with the lifting in lft,
  --   taking into account the type of mixture,
  --   that is the multiplicity of each support in sup.
  --   The array on return has the range 1..mix'length.

  function User_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Floating_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Prompts the user for a lifting value for each point in sup.
  --   Returns the lifting values.

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               stlb : double_float;
               added : Standard_Integer_Vectors.Vector )
             return Standard_Floating_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns random lifting values for each of the support sets,
  --   assigning the stable lifting bound to the artifial origins,
  --   added to compute the stable mixed volume.

  -- ON ENTRY :
  --   mix     type of mixture;
  --   sup     supports of a polynomial system;
  --   stlb    stable lifting bound, for the stable mixed volume;
  --   added   if added(k), then the artificial origin has been
  --           added to the k-th support set, as returned by the
  --           procedure Add_Artificial_Origins.

  procedure Add_Artificial_Origin
              ( dim : in integer32;
                sup : in out Lists_of_Integer_Vectors.List;
                added : out boolean );

  -- DESCRIPTION :
  --   If the origin does not belong to the list of points sup,
  --   then the origin is appended to sup and added is true on return.
  --   The dimension dim on entry equals the dimension of the points.

  procedure Add_Artificial_Origins
              ( dim : in integer32;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                nbadd : out integer32;
                added : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Add artificial origins to the supports for the stable mixed volume.

  -- ON ENTRY :
  --   dim      dimension of the points in the supports;
  --   sup      supports of a polynomial system.

  -- ON RETURN :
  --   nbadd    number of artificial origins added;
  --   added    vector of sup'range to indicate where the artificial origins
  --            are in sup: sup(k) = 1 then added to k-th support, 0 if not.

end DEMiCs_Translated_Setup;
