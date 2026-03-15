with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_DCMPLX_VecVecs_Container is

-- DESCRIPTION :
--   In the test procedures below, the verbose level is given by vrb.

  procedure Prompt_Dimensions ( vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the number of arrays and for the number
  --   of vectors in each array.

  procedure Add_Random_Vectors
              ( dim : in integer32; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   After the dimensions are set,
  --   fills up the container with random vectors of a fixed size dim.

  procedure Get_Vectors ( vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Shows the vectors stored in the container.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimensions, fills the container with
  --   random vectors, and then shows the stored vectors.

end Test_DCMPLX_VecVecs_Container;
