with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Natural_Vectors;            use Standard_Natural_Vectors;

package Checker_Moves is

-- DESCRIPTION :
--   This package defines the primitive operations to move
--   black and white checkers in a combinatorial game to count
--   the Littlewood-Richardson coefficients.

-- PART I : Moves of the Black Checkers 

  function Identity_Permutation ( n : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the identity permutation of n elements.
  --   This is called the shortest word.

  function Reverse_Permutation ( n : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the first n natural numbers in a permutation,
  --   sorted from high to low, starting at n ending at 1.
  --   This is called the longest word.

  function Descending_Checker ( p : Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns either the index of the descending black checker,
  --   or zero if the permutation is already the longest word.
  --   For use in a specializing moving flag checker game.

  function Falling_Checker ( p : Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns either the index of the descending black checker,
  --   or zero if the permutation is already the shortest word.
  --   For use in a generalizing moving flag checker game.

  function Rising_Checker ( p : Vector; i : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Given the index i to the descending black checker in p,
  --   returns the index to the rising checker.
  --   If i = 0, then 0 is returned.

  function Ascending_Checker ( p : Vector; i : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Given the index i to the falling black checker in p,
  --   returns the index to the ascending checker.
  --   If i = 0, then 0 is returned.

  procedure Move ( p : in out Vector; down,up : in integer32 );

  -- DESCRIPTION :
  --   Performs the move on the permutation p, with the index 
  --   of the descending checker in down, and the index of the
  --   rising checker in up.

  function Number_of_Moves ( n : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of moves made by n black checkers, plus one.
  --   The actual #moves is the total number of swaps to invert a sequence
  --   of n numbers, or n*(n-1)/2, but the #configurations is one more.
  --   The number returned by this function is the size of the poset.

-- PART II : Moves of the White Checkers

  function Critical_Row ( r,c : integer32; rows,cols : Vector )
                        return integer32;

  -- DESCRIPTION :
  --   Returns in the index to entry of rows of the white checker
  --   in the critical row r and column c.  If there are no white 
  --   checkers in the critical row, then the value on return is zero.

  function In_Critical_Diagonal 
              ( dr,dc,n,row,column : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there is a k, dr + k <= n and dc + k <= n, 
  --   such that row = dr + k and column = dc + k;
  --   returns false otherwise.

  -- ON ENTRY :
  --   (dr,dc)  start of the critical diagonal;
  --   n        size of the checkers board;
  --   row      row of a checker;
  --   column   column of a checker.

  function Top_White_Checker
              ( dr,dc,n : integer32; rows,cols : Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of the highest checker in the critical diagonal.
  --   Returns 0 if there is no checker in the critical diagonal.

  function Between ( a,b,x : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x lies between a and b.
  --   Notice that a and b are not necessarily sorted.

  function Blocker ( rows,cols : Vector; t,r : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there is a blocker between the top white checker
  --   indexed by t and the checker indexed by the row r.

  procedure Swap ( rows : in out Vector; i,j : in integer32 );

  -- DESCRIPTION :
  --   If i /= 0, then the rows(i) is swapped with rows(j).

  function Central_Choice
              ( p : Vector; d,r : integer32;
                rows,cols : Vector; cr,cd : integer32; debug : boolean ) 
              return integer32;
  function Central_Choice
              ( file : in file_type; p : Vector; d,r : integer32;
                rows,cols : Vector; cr,cd : integer32; debug : boolean ) 
              return integer32;

  -- DESCRIPTION :
  --   Determines whether to stay, swap, or whether to do both,
  --   depending on the position of the white checker in the
  --   critical row and the top white checker in the critical diagonal.

  -- ON ENTRY :
  --   p        permutation determines black checkers;
  --   d        descending black checker is at (p(d),p'last-d+1);
  --   r        rising black checker is at (p(r),p'last-r+1);
  --   rows     rows of the white checkers;
  --   cols     cols(cols'last-i+1) is column for checker in rows(i);
  --   cr       index to rows for white checker in critical row;
  --   cd       index to cols for top white checker in critical diagonal;
  --   debug    if true, then additional information is printed;
  --            otherwise, there is no additional output.

  -- REQUIRED :
  --   d,r,cr,cd are valid entries in the arrays (all > 0).

  -- ON RETURN :
  --   0        if the white checkers must stay;
  --   1        if the white checkers must swap;
  --   2        if choice between stay or swap.

  function Happy_in_Row ( p : Vector; r,c : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there is a black checker to the left of
  --   the white checker in row r and column c.

  function Happy_in_Column ( p : Vector; r,c : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there is a black checker above the
  --   white checker in row r and column c.

  procedure Make_Happy ( p : in Vector; rows,cols : in out Vector );

  -- DESCRIPTION :
  --   Shifts white checkers to the left or up until they are happy,
  --   or cannot be shifted any further.

  procedure Make_Happy ( p : in Vector; rows,cols : in out Vector;
                         debug : in boolean );

  -- DESCRIPTION :
  --   Same as the other "Make_Happy", prints out extra info if the
  --   flag debug is set to true.

  procedure Check_Happiness
              ( p,rows,cols : in Vector; happy : out boolean );

  -- DESCRIPTION :
  --   Verifies whether the configuration of white checkers is happy.
  --   Writes additional diagnostic information to screen.

  -- ON ENTRY :
  --   p        permutation indicating the black checkers;
  --   rows     rows of the white checkers;
  --   cols     columns of the white checkers.

  -- ON RETURN :
  --   happy    true if all white checkers are happy,
  --            false otherwise.

  function Happy_Checkers ( p,rows,cols : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true for a happy configuration of white checkers,
  --   returns false otherwise, when there is no solution to the
  --   given intersection condition.

  -- ON ENTRY :
  --   p        permutation indicating the black checkers;
  --   rows     rows of the white checkers;
  --   cols     columns of the white checkers.

end Checker_Moves;
