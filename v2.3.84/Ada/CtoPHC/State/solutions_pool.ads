with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Solutions_Pool is

-- DESCRIPTION :
--   The solutions pool manages several lists of solutions.
--   The design of the pool was modeled after Solutions_Container.

-- CREATORS :

  procedure Initialize ( n : in integer32 );

  -- DESCRIPTION :
  --   Initializes the pool for n solution lists.

  procedure Initialize ( k : in integer32; sols : in Solution_List );

  -- DESCRIPTION :
  --   Initializes the k-th solution list with the solutions in sols.

-- SELECTORS :

  function Size return natural32;

  -- DESCRIPTION :
  --   Returns the size of the solutions pool.

  function Length ( k : in integer32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the length of the k-th solution list in the pool.

  function Dimension ( k : in integer32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the k-th solution list in the pool.

  function Retrieve ( k : in integer32 ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the k-th solution list in the pool.

  procedure Retrieve ( k : in integer32; i : in natural32;
                       s : out Solution; fail : out boolean );
  procedure Retrieve ( k : in integer32; i : in natural32;
                       s : out Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Returns in s the i-th solution of the k-th solution list of the pool,
  --   fail on return is false if k > Size or i > Length(k).

  procedure Replace ( k : in integer32; i : in natural32;
                      s : in Solution; fail : out boolean );
  procedure Replace ( k : in integer32; i : in natural32;
                      s : in Link_to_Solution; fail : out boolean );

  -- DESCRIPTION :
  --   Replaces the i-th solution of the k-th solution list of the pool by s,
  --   fail on return is false if k > Size or i > Length(k).

  procedure Append ( k : in integer32; s : in Solution );
  procedure Append ( k : in integer32; s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Appends the solution to the k-th list in the container.

-- DESTRUCTORS :

  procedure Clear ( k : in integer32 );

  -- DESCRIPTION :
  --   Clears the k-th solution list in the pool.

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the entire pool of all solution lists.

end Solutions_Pool;
