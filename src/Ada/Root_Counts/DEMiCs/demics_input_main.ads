with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with demics_input_data;

package DEMiCs_Input_Main is

-- DESCRIPTION :
--   Exports the main input data procedures.

  procedure read_data_from_file
              ( data : out demics_input_data.class_dataSet.dataSet;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the name of the data input file for demics
  --   and then prints the data read from the file,
  --   if the verbose level vrblvl > 0.

  procedure interactive_input_data
              ( data : out demics_input_data.class_dataSet.dataSet;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for input data, dimension, number of distinct supports,
  --   the number of occurrences of each supports, the number of points
  --   in each support, and then prompt for each point in the supports.

end DEMiCs_Input_Main;
