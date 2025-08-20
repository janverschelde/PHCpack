package Test_DEMiCs_Input is

-- DESCRIPTION :
--   Tests the getting of the input data for DEMiCs.

  procedure read_data_from_file;

  -- DESCRIPTION :
  --   Prompts for the name of the data input file for demics
  --   and then prints the data read from the file.

  procedure interactive_input_data;

  -- DESCRIPTION :
  --   Prompts for input data, dimension, number of distinct supports,
  --   the number of occurrences of each supports, the number of points
  --   in each support, and then prompt for each point in the supports.

  procedure Main;

  -- DESCRIPTION :
  --   The main test runs one of the two test procedures.

end Test_DEMiCs_Input;
