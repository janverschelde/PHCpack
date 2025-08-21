with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;

package demics_input_data is

-- DESCRIPTION :
--   Originally generated via gcc -c -fdump-ada-spec on inputData.h
--   defines the input data for DEMiCs.

  package class_dataSet is

    type dataSet is record 
      dim : integer32;           -- dimension of the points
      supN : integer32;          -- number of distinct supports
      termSumNum : integer32;    -- total number of points in all supports  
      termMax : integer32;       -- largest number in the termSet vector
      typeMax : integer32;       -- largest number in the supType vector
      termSet : Standard_Integer_Vectors.Link_to_Vector;
       -- vector of range 1..supN with the number of points
       -- in each distinct support set
      termStart : Standard_Integer_Vectors.Link_to_Vector;
       -- vector of range 1..supN+1 with the index to the first point in 
       -- each support set, the last element in termStart is the total 
       -- number of points
      supType : Standard_Integer_Vectors.Link_to_Vector;
       -- vector of range 1..supN with the number of occurrences
       -- in each distinct support set
      support : Standard_Floating_Vectors.Link_to_Vector; 
       -- coordinates of the points in each support
      coef : Standard_Floating_Vectors.Link_to_Vector;
      outFile : Link_to_String;  -- name of the output file
    end record;

    function new_dataSet return dataSet;

    -- DESCRIPTION :
    --   Initializes an empty data set, setting all values to zero
    --   and all pointers to null.

    procedure delete_dataSet ( this : in out dataSet );

    -- DESCRIPTION :
    --   Deallocates all vectors in the data set.

    procedure support_in ( this : in dataSet;
                           rowIdx : in integer32; colIdx : in integer32;
                           elem : in double_float );

    -- DESCRIPTION :
    --   Assigns the elem to the coordinate of position colIdx
    --   of support with index rowIdx.  All indices start at 1.

    function support_out ( this : in dataSet;
                           rowIdx : in integer32; colIdx : in integer32 )
                         return double_float;
    -- DESCRIPTION :
    --   Returns the coordinate of index colIdx
    --   of the point with index rowIdx.  All indices start at 1.

    function makeOutputFile ( inputFile : Link_to_String )
                            return Link_to_String;

    -- DESCRIPTION :
    --   The string on return has everything before the dot in inputFile,
    --   or entirely copies the input file if there is not dot.
    --   The ending of the file on return is ".out" for use in the outFile
    --   field of the data set.

    procedure parse_preamble ( file : in file_type;
                               this : in out dataSet; fail : out boolean );

    -- DESCRIPTION :
    --   Parses the file, opened for input, for the preamble:
    --   dimension, number of distinct supports, number of elements
    --   in each support, and the occurrences of each support.
    --   If the formatting is wrong, then fail is true on return,
    --   otherwise, the data set has the preamble defined.

    -- EXAMPLE : the preamble of cyclic5 is below
    --
    -- Dim = 5
    -- Support = 5
    -- 
    -- Elem = 2 5 5 5 5 
    -- Type = 1 1 1 1 1 

    procedure parse_supports ( file : in file_type;
                               this : in out dataSet; fail : out boolean );

    -- DESCRIPTION :
    --   Parses the file, opened for input, for the supports,
    --   which can only be done if parse_preamble did not fail.
    --   If the formatting is wrong, then fail is true on return,
    --   otherwise, the data set has the supports defined.

    -- EXAMPLE : the supports of cyclic5 are below
    --
    -- 0 0 0 0 0 
    -- 1 1 1 1 1 
    -- 
    -- 1 1 1 1 0 
    -- 0 1 1 1 1 
    -- 1 0 1 1 1 
    -- 1 1 0 1 1 
    -- 1 1 1 0 1 
    -- 
    -- 1 1 1 0 0 
    -- 0 1 1 1 0 
    -- 0 0 1 1 1 
    -- 1 0 0 1 1 
    -- 1 1 0 0 1 
    -- 
    -- 1 1 0 0 0 
    -- 0 1 1 0 0 
    -- 0 0 1 1 0 
    -- 0 0 0 1 1 
    -- 1 0 0 0 1 
    -- 
    -- 1 0 0 0 0 
    -- 0 1 0 0 0 
    -- 0 0 1 0 0 
    -- 0 0 0 1 0 
    -- 0 0 0 0 1 

    procedure getInputFile ( this : in out dataSet;
                             inputFile : in Link_to_String;
                             fail : out boolean );

    -- DESCRIPTION :
    --   Opens the file with name defined by inputFile, parses the 
    --   information on the file into the data on the support sets.
    --   Stores the ".out" file in the outFile of this,
    --   using everything in front of the dot in the inputFile name.
    --   If the file cannot be opened for input, or the format is wrong,
    --   then fail is true on return.

    procedure info_preamble ( this : in dataSet );

    -- DESCRIPTION :
    --   Writes the dimension, number of distinct supports,
    --   the number of points in each support and
    --   the number of occurrences of each support to screen.

    procedure info_supports ( this : in dataSet );

    -- DESCRIPTION :
    --   Writes each support set to screen.

  end class_dataSet;

end demics_input_data;
