with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with DoblDobl_Complex_Solutions;        use DoblDobl_Complex_Solutions;

package DoblDobl_Solutions_Container is

-- DESCRIPTION :
--   This package provides a container for a solution list
--   with double double complex coefficients for the interface with C.

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

end DoblDobl_Solutions_Container;
