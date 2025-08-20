with String_Splitters;                   use String_Splitters;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;

package demics_input_data is

-- DESCRIPTION :
--   Originally generated via gcc -c -fdump-ada-spec on inputData.h
--   defines the input data for DEMiCs.

  package class_dataSet is

    type dataSet is record         -- omitted limited before record
      dim : integer32;                 -- inputData.h:32
      supN : integer32;                -- inputData.h:33
      termSumNum : integer32;          -- inputData.h:34
      termMax : integer32;             -- inputData.h:35
      typeMax : integer32;             -- inputData.h:36
      termSet : Standard_Integer_Vectors.Link_to_Vector;   -- inputData.h:38
      termStart : Standard_Integer_Vectors.Link_to_Vector; -- inputData.h:40
      c_type : Standard_Integer_Vectors.Link_to_Vector;    -- inputData.h:43
      support : Standard_Floating_Vectors.Link_to_Vector;  -- inputData.h:46
      coef : Standard_Floating_Vectors.Link_to_Vector;     -- inputData.h:47
      outFile : Link_to_String;        -- inputData.h:49
    end record;

    function new_dataSet return dataSet;  -- inputData.h:29

    procedure delete_dataSet ( this : access dataSet );  -- inputData.h:30

    procedure support_in ( this : access dataSet;
                           rowIdx : integer32; colIdx : integer32;
                           elem : double_float );  -- inputData.h:51

    function support_out ( this : access dataSet;
                           rowIdx : integer32; colIdx : integer32 )
                         return double_float;  -- inputData.h:60

    procedure getInputFile ( this : access dataSet;
                             inputFile : Link_to_String ); -- inputData.h:69

    procedure info_preamble ( this : access dataSet );  -- inputData.h:76

    procedure info_supports ( this : access dataSet );  -- inputData.h:83

  end class_dataSet;

  use class_dataSet;

end demics_input_data;
