with text_io;                           use text_io;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Standard_Natural_Matrices;         use Standard_Natural_Matrices;
with Checker_Boards;                    use Checker_Boards;

package Checker_Boards_io is

-- DESCRIPTION :
--   This package offers some simple input and output facilities
--   for permutations and checker boards.

  procedure Write ( p : in Vector );
  procedure Write ( file : in file_type; p : in Vector );

  -- DESCRIPTION :
  --   if p'last < 10, then this is the usual put,
  --   otherwise spacings between the entries are adjusted.

  procedure Write ( b : in Board );
  procedure Write ( file : in file_type; b : in Board );

  -- DESCRIPTION :
  --   Writes the board b to standard output or to file.

  procedure Write ( a,b : in Board );
  procedure Write ( file : in file_type; a,b : in Board );

  -- DESCRIPTION :
  --   Writes the boards a and b to standard output or to file.

  procedure Write ( b : in Board; f : in Matrix );
  procedure Write ( file : in file_type; b : in Board; f : in Matrix );

  -- DESCRIPTION :
  --   Writes the board b and moving flag f to standard output or to file.

  -- REQUIRED : b'range(k) = f'range(k), for k in 1..2.

  procedure Write ( b : in Board; f,t : in Matrix );
  procedure Write ( file : in file_type; b : in Board; f,t : in Matrix );

  -- DESCRIPTION :
  --   Writes the board b, moving flag f, the transformation t,
  --   and its product f*t to standard output or to file.

  -- REQUIRED : b'range(k) = f'range(k) = t'range(k), for k in 1..2.

  procedure Write ( b : in Board; f,t,x : in Matrix );
  procedure Write ( file : in file_type; b : in Board; f,t,x : in Matrix );

  -- DESCRIPTION :
  --   Writes board, localization pattern, moving flag, transformation,
  --   and result of the transformation on the moving flag.

  -- ON ENTRY :
  --   b         checker board configuration with black & white checkers;
  --   f         natural representation of the moving flag;
  --   t         natural representation of the transformation;
  --   x         column localization pattern of the solution k-plane.

  procedure Read_Permutation ( p : out Vector );

  -- DESCRIPTION :
  --   Prompts the user to give a permutation vector.

  procedure Write_Permutation ( p : in Vector );
  procedure Write_Permutation ( file : in file_type; p : in Vector );
  procedure Write_Permutation ( p : in Vector; f : in Matrix );
  procedure Write_Permutation ( file : in file_type;
                                p : in Vector; f : in Matrix );
  procedure Write_Permutation ( p : in Vector; f,t : in Matrix );
  procedure Write_Permutation ( file : in file_type;
                                p : in Vector; f,t : in Matrix );
  procedure Write_Permutation ( p : in Vector; f,t,x : in Matrix );
  procedure Write_Permutation ( file : in file_type;
                                p : in Vector; f,t,x : in Matrix );
  procedure Write_Permutation ( p,r,c : in Vector; f : in Matrix );
  procedure Write_Permutation ( file : in file_type;
                                p,r,c : in Vector; f : in Matrix );
  procedure Write_Permutation ( p,r,c : in Vector; f,t : in Matrix );
  procedure Write_Permutation ( file : in file_type;
                                p,r,c : in Vector; f,t : in Matrix );
  procedure Write_Permutation ( p,r,c : in Vector; f,t,x : in Matrix );
  procedure Write_Permutation ( file : in file_type;
                                p,r,c : in Vector; f,t,x : in Matrix );

  -- DESCRIPTION :
  --   Writes the permutation p, the corresponding checker board,
  --   with white checkers at rows r and column c,
  --   the localization pattern for the solution k-plane x,
  --   and also the patterns of the moving flag f and the 
  --   transformation t, if provided.

  procedure Write_Bracket ( v : in Vector );
  procedure Write_Bracket ( file : in file_type; v : in Vector );

  -- DESCRIPTION :
  --   Writes a vector, enclosed by square brackets.
  --   Skips the leading space when writing a vector.

  procedure Write_Coordinates ( rows,cols : in Vector );
  procedure Write_Coordinates ( file : in file_type; rows,cols : in Vector );

  -- DESCRIPTION :
  --   Writes the rows and columns which define the white checkers.

  procedure Write ( b : in Board; p,rows,cols : in Vector );
  procedure Write ( file : in file_type;
                    b : in Board; p,rows,cols : in Vector );

  -- DESCRIPTION :
  --   Writes the board and the vectors defining the positions
  --   of the black and white checkers.

  procedure Write ( b : in Board; f : in Matrix; p,rows,cols : in Vector );
  procedure Write ( file : in file_type;
                    b : in Board; f : in Matrix; p,rows,cols : in Vector );

  -- DESCRIPTION :
  --   Writes the board, the pattern for the moving flag, and the vectors
  --   defining the positions of the black and white checkers.

  procedure Write_Permutations ( p,q,pr,pc,qr,qc : in Vector );
  procedure Write_Permutations
              ( file : in file_type; p,q,pr,pc,qr,qc : in Vector );

  -- DESCRIPTION :
  --   Writes the two boards side by side, defined by permutations p and q,
  --   along with corresponding white checker localizations with rows in
  --   pr and qr, and columns in pc and qc.

  procedure Write_Patterns ( fp,fq,xp,xq : in Matrix );
  procedure Write_Patterns ( file : in file_type; fp,fq,xp,xq : in Matrix );

  -- DESCRIPTION :
  --   Writes two consecutives flags with patterns in fp and fq, with
  --   corresponding column localization patterns of k-planes xp and xq.

  procedure Write_Permutations_and_Patterns
              ( p,q,pr,pc,qr,qc : in Vector; fp,fq,xp,xq : in Matrix );
  procedure Write_Permutations_and_Patterns
              ( file : in file_type;
                p,q,pr,pc,qr,qc : in Vector; fp,fq,xp,xq : in Matrix );

  -- DESCRIPTION :
  --   Writes permutations and patterns to standard output or file.

end Checker_Boards_io;
