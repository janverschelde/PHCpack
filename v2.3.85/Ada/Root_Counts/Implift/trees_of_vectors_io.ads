with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Trees_of_Vectors;                   use Trees_of_Vectors;

package Trees_of_Vectors_io is

-- DESCRIPTION :
--   This package provides i/o operations for Trees_of_Vectors.
--   The formats for input are the same as those for output.

  procedure get ( n : in natural32; tv : in out Tree_of_Vectors );
  procedure get ( file : in file_type;
                  n : in natural32; tv : in out Tree_of_Vectors );

  -- DESCRIPTION :
  --   Reads a tree of vectors from standard input or from file,
  --   n equals the maximum length of a vectors in the tree.

  procedure put ( tv : in Tree_of_Vectors );
  procedure put ( file : in file_type; tv : in Tree_of_Vectors );

  -- DESCRIPTION :
  --   Writes a tree of vectors on standard output or on file.
  --   The left branch in the tree comes first,
  --   for each vector a new line is taken,
  --   the entries in the vector are separated by spaces.
  --   The output of a whole tree is ended by a blank line.

end Trees_of_Vectors_io;
