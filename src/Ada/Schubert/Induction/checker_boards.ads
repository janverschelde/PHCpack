with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;

package Checker_Boards is

-- DESCRIPTION :
--   This package defines a checker board and provides some
--   operations regarding placing checkers on the board.
--   Black checkers are identified by their row on the board.
--   The k-th black checker is at row k = p(i) and column p'last+1-k.

  type Board is array ( integer32 range <>, integer32 range <> ) of character;

  procedure Initialize ( b : in out Board );

  -- DESCRIPTION :
  --   Initializes the board to spaces.

  function Configuration ( p : Vector ) return Board;
  
  -- DESCRIPTION :
  --   Returns the configuration corresponding to the permutation p.
  --   If a black checker is in row r and column c,
  --   then the permutation p sends n+1-c to r, where n = p'last.

  function Column ( p : Vector; r : natural32 ) return integer32;

  -- DESCRIPTION :
  --   Returns p'last+1-i where p(i) = r.
  --   If no i exists for which p(i) = r, then 0 is returned.

  function Permutation ( b : Board ) return Vector;

  -- DESCRIPTION :
  --   Returns the permutation corresponding to the given board.
  --   If the (i,j)-th position on the board is not blank,
  --   then the permutation p on return will have i at entry n+1-j,
  --   i.e.: p(n+1-j) = i, for a black checker in row i and column j.

  -- NOTE :
  --   The operations Configuration and Permutation are connected:
  --   for any permutation p, Permutation(Configuration(p)) = p; and
  --   for any board b, Configuration(Permutation(b)) = b.

  procedure Place_White ( b : in out Board; r,c : in Vector );

  -- DESCRIPTION :
  --   Places white checkers given by rows r and columns c
  --   on the board b.

  function Configuration ( p,r,c : Vector ) return Board;

  -- DESCRIPTION :
  --   Returns the configuration defined by black checkers in p
  --   and white checkers given by r and c.

end Checker_Boards;
