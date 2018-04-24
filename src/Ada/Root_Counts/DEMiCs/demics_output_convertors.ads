with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package DEMiCs_Output_Convertors is

-- DESCRIPTION :
--   The output of DEMiCs is converted into a mixed subdivision.

  function Apply_Lifting
             ( sup : Lists_of_Integer_Vectors.List;
               lif : Standard_Floating_Vectors.Vector )
             return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION :
  --   Given in sup is the support of a polynomial
  --   and in lif are the corresponding lifting values.
  --   Returns the support in sup, with the lifting applied
  --   as defined by the lifting values in lif.

  -- REQUIRED : Length_Of(sup) = lif'last.

  function Apply_Lifting
             ( sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lif : Standard_Floating_VecVecs.VecVec )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Given in sup are the supports of a polynomial system
  --   and in lif are the corresponding lifting values.
  --   Returns the supports in sup, with the lifting applied
  --   as defined by the lifting values in lif.

  -- REQUIRED :
  --   for k in lif'range: length_of(sup(k)) = lif(k)'last.

  function Make_Mixed_Cell
             ( dim : integer32;
               lbl : Standard_Integer_Vectors.Vector;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
             return Mixed_Cell;

  -- DESCRIPTION :
  --   Returns the mixed cell corresponding to the labeled points.

  -- ON ENTRY :
  --   dim     dimension before lifting;
  --   lbl     labels to the points that span the mixed cell;
  --   lifsup  lifted supports.

  function Make_Mixed_Cells
             ( dim : integer32;
               labels : Lists_of_Integer_Vectors.List;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
             return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the mixed cells corresponding to the labeled points.

  -- ON ENTRY :
  --   dim     dimension before lifting;
  --   labels  labels to the points that span the mixed cell;
  --   lifsup  lifted supports.

end DEMiCs_Output_Convertors;
