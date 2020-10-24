with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Test_Balance_Lifting is

-- DESCRIPTION :
--   Test balancing of lifting values.

  procedure Read_Subdivision
             ( n : out natural32;
               mix : out Standard_Integer_Vectors.Link_to_Vector;
               mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Prompts for a mixed cell configuration.

  function Inner_Product
             ( point : Standard_Floating_Vectors.Link_to_Vector;
               lifvals,normal : Standard_Floating_Vectors.Link_to_Vector )
             return double_float;

  -- DESCRIPTION :
  --   Returns the inner product of the normal with the point
  --   with lifting values in lifvals.

  procedure Extremal_Inner_Products
             ( points : in Lists_of_Floating_Vectors.List;
               lifvals,normal : in Standard_Floating_Vectors.Link_to_Vector;
               min,max : out double_float );

  -- DESCRIPTION :
  --   Returns in max and min the largest and smalles inner product of 
  --   the normal with the points with lifting values in lifvals.

  function Maximal_Inner_Product
             ( points : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               lifvals,normal : Standard_Floating_Vectors.Link_to_Vector )
             return double_float;

  -- DESCRIPTION :
  --   Returns the largest inner product of the normal with the points
  --   with lifting values in lifvals.

  procedure Balance_Lifting_Values
              ( file : in file_type;
                n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   The lifted points are extracted from the mixed-cell configuration.
  --   Every point has a unique index and the lifting value in every point
  --   is replaced by the index of the point.
  --   The lifting values of the points are stored in one central vector.

  procedure Main;

  -- DESCRIPTION :
  --   Reads a mixed cell configuration and balances the lifting values.

end Test_Balance_Lifting;
