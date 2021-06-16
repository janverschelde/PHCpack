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

  procedure Retrieve_Next_Initialize;

  -- DESCRIPTION :
  --   Resets the current pointer to the first solution in the container.

  procedure Retrieve_Next ( s : out Link_to_Solution; k : out natural32 );

  -- DESCRIPTION :
  --   Retrieves the next solution, where the current pointer is,
  --   and returns in k the index of the solution in the list.
  --   If k > 0, then s contains the pointer to the current solutioni
  --   and the current pointer is moved to the next solution.
  --   Otherwise, if the current pointer is null, then k on return is 0.

  procedure Retrieve_Current ( s : out Link_to_Solution; k : out natural32 );

  -- DESCRIPTION :
  --   Retrieves the current solution without advancing the pointer to
  --   the current solution.  If the pointer refers to a valid solution,
  --   then k is the solution number on return, otherwise, k = 0.

  procedure Move_Current ( k : out natural32 );

  -- DESCRIPTION :
  --   Moves the pointer to the current solution to the next solution.
  --   If the pointer is currently null or if there is no next solution,
  --   then k will be zero on return.  Otherwise, k will equal the
  --   value of the cursor, the number of the current solution.

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
