with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_Double_VecVecs_Container is

  procedure Prompt_Dimensions;

  -- DESCRIPTION :
  --   Prompts for the number of arrays and for the number
  --   of vectors in each array.

  procedure Add_Random_Vectors ( dim : in integer32 );

  -- DESCRIPTION :
  --   After the dimensions are set,
  --   fills up the container with random vectors of a fixed size dim.

  procedure Get_Vectors;

  -- DESCRIPTION :
  --   Shows the vectors stored in the container.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimensions, fills the container with
  --   random vectors, and then shows the stored vectors.

end Test_Double_VecVecs_Container;
