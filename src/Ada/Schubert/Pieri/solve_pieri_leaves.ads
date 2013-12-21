with text_io;                            use text_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Brackets;                           use Brackets;

function Solve_Pieri_Leaves
              ( file : file_type; b1,b2 : Bracket; m : Matrix ) return Matrix;

-- DESCRIPTION :
--   Computes a solution at a pair of brackets at the leaves of the Pieri tree.

-- REQUIRED : (b1,b2) satisfy Pieri's condition.

-- ON ENTRY :
--   file       for intermediate output;
--   b1         bracket at leaf of first Pieri tree;
--   b2         bracket at leaf of second Pieri tree;
--   m          m-plane.

-- ON RETURN :
--   The p-plane that intersects m nontrivially, p = b1'length = b2'length.
