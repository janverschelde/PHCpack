with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Multprec_Complex_Solutions;        use Multprec_Complex_Solutions;

package Multprec_Solutions_Container is

-- DESCRIPTION :
--   This package provides a container for a solution list
--   with multiprecision complex coefficients for the interface with C.

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

end Multprec_Solutions_Container;
