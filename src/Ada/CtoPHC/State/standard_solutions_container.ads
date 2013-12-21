with text_io;                           use text_io;
with Symbol_Table;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_Solutions_Container is

-- DESCRIPTION :
--   This package provides a container for a solution list,
--   designed for the interface with C.

  procedure Initialize ( sols : in Solution_List ); 

  -- DESCRIPTION :
  --   Initializes the container with a solution list.

  function Length return natural32;

  -- DESCRIPTION :
  --   Returns the number of solutions in the container.

  function Dimension return natural32;

  -- DESCRIPTION :
  --   Returns 0 if the container is empty, otherwise the
  --   dimension of the solution vectors is returned.

  function Retrieve return Solution_List;

  -- DESCRIPTION :
  --   Returns the lists of solutions stored in the container.

  procedure Retrieve ( k : in natural32; s : out Solution;
                       fail : out boolean );
  procedure Retrieve ( k : in natural32; s : out Link_to_Solution;
                       fail : out boolean );
  -- DESCRIPTION :
  --   Attemps to retrieve the k-th solution in s.
  --   Reports failure (fail = true) if there are fewer than k
  --   solutions in the container.

  procedure Replace ( k : in natural32; s : in Solution;
                      fail : out boolean );
  procedure Replace ( k : in natural32; s : in Link_to_Solution;
                      fail : out boolean );
  -- DESCRIPTION :
  --   Attemps to replace the k-th solution with the given one.
  --   Reports failure (fail = true) if there are fewer than k
  --   solutions in the container.
	  
  procedure Append ( s : in Solution );
  procedure Append ( s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Appends the solution to the container.

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the solutions in the container.

-- Handling files for incremental read and write of solution lists :

  procedure Silent_Open_Input_File;
  procedure Silent_Open_Input_File ( k : in natural32 );
  procedure Silent_Open_Input_File ( filename : in string );
  procedure Silent_Open_Input_File ( k : in natural32; filename : in string );
  procedure Open_Input_File;
  procedure Open_Input_File ( k : in natural32 );

  -- DESCRIPTION :
  --   Reads the name of a file for input of solutions.
  --   The "Silent_" one just asks for a string of characters,
  --   while the other one prints first a message before prompting.
  --   When the filename is provided, no messages are written.
  --   If k is provided, then it is for reading witness set k.

  procedure Create_Output_File;

  -- DESCRIPTION :
  --   Reads the name of a file for output of solutions.

  function Solution_Input_File return file_type;
  function Solution_Input_File ( k : natural32 ) return file_type;

  function Solution_Output_File return file_type;

  -- DESCRIPTION :
  --   Returns the variable file for input or output.

  procedure Reset_Input_File ( k : in natural32 );

  -- DESCRIPTION :
  --   Resets the file for the k-th witness set for reading.
  --   The number k must be 1 or 2 to have an effect.

  procedure Close_Input_File;
  procedure Close_Input_File ( k : in natural32 );
  procedure Close_Output_File;

  -- DESCRIPTION :
  --   Closes file used for input or output.

-- Management of two symbol tables for diagonal homotopies :

  procedure Store_Symbol_Table
              ( k : in natural32; sbs : Symbol_Table.Array_of_Symbols );

  -- DESCRIPTION :
  --   Stores the symbol table for the k-th witness set, k must be 1 or 2.
  --   The k = 0 can be used for a common symbol table for reading the
  --   solutions in the witness sets from file.

  function Retrieve_Symbol_Table 
             ( k : in natural32 ) return Symbol_Table.Link_to_Array_of_Symbols;

  -- DESCRIPTION :
  --   Returns the symbol table for the k-th witness set, k must be 1 or 2.
  --   With k = 0, one obtains the common symbol table to read solutions.

  procedure Clear_Symbol_Table ( k : in natural32 );

  -- DESCRIPTION :
  --   Clears the symbol table for the k-th witness set, if k = 1 or 2;
  --   or, in case k = 0, then the common solution symbol table is cleared.

end Standard_Solutions_Container;
