with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;

package DoblDobl_Complex_Solutions_io is

-- DESCRIPTION :
--   This routines provides routines for input and output of solutions.

-- FOR SOLUTION VECTORS ONLY :

  procedure get_vector ( s : in out Solution );
  procedure get_vector ( file : in file_type; s : in out Solution );

  -- DESCRIPTION :
  --   The input must contain the solution vector.

  procedure put_vector ( v : in DoblDobl_Complex_Vectors.Vector );
  procedure put_vector ( file : in file_type; 
                         v : in DoblDobl_Complex_Vectors.Vector );
  procedure put_vector ( s : in Solution );
  procedure put_vector ( file : in file_type; s : in Solution );

  -- DESCRIPTION :
  --   On the output the solution vector will be written.

  procedure put_diagnostics ( s : in Solution );
  procedure put_diagnostics ( file : in file_type; s : in Solution );

  -- DESCRIPTION :
  --   Writes the diagnostics of the solution to screen or to file.

-- FOR SOLUTIONS :

  procedure get ( s : out Solution );
  procedure get ( n : in natural32; ls : out Link_to_Solution );
  procedure get ( file : in file_type; s : out Solution );
  procedure get ( file : in file_type;
                  n : in natural32; ls : out Link_to_Solution );

  -- DESCRIPTION :
  --   The input must contain the following : s.t, s.m and s.v(i), 
  --   a vector of s.n complex numbers

  procedure put ( s : in Solution );
  procedure put ( file : in file_type; s : in Solution );

  -- DESCRIPTION :
  --   On the output the following will be written :
  --   s.t, s.m and s.v, a vector of s.n complex numbers

-- FOR LISTS OF SOLUTIONS :

  procedure get ( sols : in out Solution_List );
  procedure get ( sols,sols_last : in out Solution_List );
  procedure get ( len,n : in natural32; sols : in out Solution_List );
  procedure get ( len,n : in natural32; sols,sols_last : in out Solution_List );
  procedure get ( file : in file_type; sols : in out Solution_List );
  procedure get ( file : in file_type; sols,sols_last : in out Solution_List );
  procedure get ( file : in file_type; len,n : in natural32;
                  sols : in out Solution_List );
  procedure get ( file : in file_type; len,n : in natural32;
                  sols,sols_last : in out Solution_List );

  -- DESCRIPTION :
  --   A solution list will be read.  If the length len and dimension n
  --   of the list is not supplied, then they will be read first.
  --   If the parameter sols_last is supplied, then this parameter contains
  --   the pointer to the last element of the list on return.
  --   The solutions should be in the appropriate format.

  procedure put ( sols : in Solution_List );
  procedure put ( len,n : in natural32; sols : in Solution_List );
  procedure put ( file : in file_type; sols : in Solution_List );
  procedure put ( file : in file_type; len,n : in natural32;
                  sols : in Solution_List );

  -- DESCRIPTION :
  --   The solutions are written on standard output or on file.
  --   First the length of the list and the dimension of the solutions
  --   will be put on file if they are supplied as parameter.

-- USER-FRIENDLY ROUTINES :

  procedure Display_Format;

  -- DESCRIPTION :
  --   Displays on screen the formatting rules as on-line help facility.

  procedure Read ( sols : in out Solution_List );

  -- DESCRIPTION :
  --   Reads the solution list from file, displays the formatting information
  --   in case of exception and let the user try again.

  procedure Write ( file : in file_type; sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the solution list on file, preceded by length and dimension.
  --   Nothing is written in case the solution list is empty.

-- INCREMENTAL READ and WRITE (for huge solution lists):

  procedure Read_First ( file : in file_type; len,dim : out natural32 );

  -- DESCRIPTION :
  --   Reads the length of the total solution list and the
  --   length of the solution vectors from file.

  -- REQUIRED :
  --   The file must be opened for input and its current position must be
  --   at the beginning of the solution list, starting with the 2 numbers:
  --   length of the list and dimension of the vectors.

  -- ON ENTRY :
  --   file      input file where to get the solutions.

  -- ON RETURN :
  --   len       length of the list;
  --   dim       length of the solution vectors in the list.

  procedure Read_First ( file : in file_type; len,dim : out natural32;
                         s : out Link_to_Solution );

  -- DECRIPTION :
  --   Reads the length of the total solution list, the dimension of
  --   the solution vectors, and the first solution from file.

  -- REQUIRED :
  --   The file must be opened for input and its current position must be
  --   at the beginning of the solution list, starting with the 2 numbers:
  --   length of the list and dimension of the vectors.

  -- ON ENTRY :
  --   file      input file where to get the solutions.

  -- ON RETURN :
  --   len       length of the list;
  --   dim       length of the solution vectors in the list;
  --   s         pointer to the first solution in the list.

  procedure Read_First ( file : in file_type; n : in natural32;
                         len,dim,nbr : out natural32;
                         first,last : in out Solution_List );

  -- DESCRIPTION :
  --   Reads the first n solutions from file.
 
  -- REQUIRED :
  --   The file must be opened for input and its current position must be
  --   at the beginning of the solution list, starting with the 2 numbers:
  --   length of the list and dimension of the vectors.

  -- ON ENTRY :
  --   file      input file where to get the solutions;
  --   n         number of solutions to read from file;
  --   first     current list of solutions (may be empty);
  --   last      pointer to last element in the list first. 

  -- ON RETURN :
  --   len       length of the list;
  --   dim       length of the solution vectors in the list;
  --   nbr       number of solutions read, normally: nbr = min(len,n),
  --             nbr < n if the list on file has fewer than n solutions;
  --   first     the new solutions are appended to the list;
  --   last      updated pointer to last element of first.

  procedure Write_First ( file : in file_type; len,dim : in natural32 );

  -- DESCRIPTION :
  --   Writes the dimensions of the solution list to file.

  -- REQUIRED :
  --   The file must be opened for output.

  -- ON ENTRY :
  --   file      output file where to put the solutions;
  --   len       total number of solutions which will be written to file,
  --             this should not necessarily be the length of the sols;
  --   dim       length of the solution vectors in the solution list.

  procedure Write_First ( file : in file_type; len,dim : in natural32;
                          s : in Solution );
  procedure Write_First ( file : in file_type; len,dim : in natural32;
                          s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Writes the first solution of the list to file.

  -- REQUIRED :
  --   The file must be opened for output.

  -- ON ENTRY :
  --   file      output file where to put the solutions;
  --   len       total number of solutions which will be written to file,
  --             this should not necessarily be the length of the sols;
  --   dim       length of the solution vectors in the solution list;
  --   s         solution to be written.

  procedure Write_First ( file : in file_type; n,len,dim : in natural32;
                          nbr : out natural32; sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the first n solutions of the list to file.

  -- REQUIRED :
  --   The file must be opened for output.

  -- ON ENTRY :
  --   file      output file where to put the solutions;
  --   n         number of solutions to write to file;
  --   len       total number of solutions which will be written to file,
  --             this should not necessarily be the length of the sols;
  --   dim       length of the solution vectors in the solution list;
  --   sols      list of solutions to be written.

  -- ON RETURN :
  --   nbr       number of solutions written to file, this could be less
  --             than n if Length_Of(sols) < n.

  procedure Read_Next ( file : in file_type; s : out Solution );
  procedure Read_Next ( file : in file_type; dim : in natural32;
                        s : out Link_to_Solution );

  -- DESCRIPTION :
  --   Reads the next solution from file.

  -- REQUIRED :
  --   The file must be opened for input and its current position should
  --   be after the length and dimension of the list.

  -- ON ENTRY :
  --   file      input file where to get the solutions;
  --   dim       length of the solution vectors in the list.

  -- ON RETURN :
  --   s         pointer to the next solution in the list.

  procedure Read_Next ( file : in file_type; n,dim : in natural32;
                        nbr : out natural32;
                        first,last : in out Solution_List );

  -- DESCRIPTION :
  --   Reads the next n solutions from file.

  -- REQUIRED :
  --   The file must be opened for input and its current position should
  --   be after the length and dimension of the list.

  -- ON ENTRY :
  --   file      input file where to get the solutions;
  --   n         number of solutions to read from file;
  --   dim       length of the solution vectors in the list;
  --   first     current list of solutions (may be empty);
  --   last      pointer to last element in the list first. 

  -- ON RETURN :
  --   nbr       number of solutions read, nbr < n if the list 
  --             on file has fewer than n solutions left;
  --   first     the new solutions are appended to the list;
  --   last      updated pointer to last element of first.

  procedure Write_Next ( file : in file_type; cnt : in out natural32;
                         s : in Solution );
  procedure Write_Next ( file : in file_type; cnt : in out natural32;
                         s : in Link_to_Solution );

  -- DESCRIPTION :
  --   Writes the next solution to file.

  -- REQUIRED :
  --   The file must be opened for output.  For format consistency,
  --   a successful Write_First must have been executed.

  -- ON ENTRY :
  --   file      file opened for output;
  --   cnt       number of solutions already written to file;
  --   s         solution to be written to file.

  -- ON RETURN :
  --   cnt       updated counter for number of solutions written to file.

  procedure Write_Next ( file : in file_type; n : in natural32;
                         nbr : out natural32; cnt : in out natural32;
                         sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the next n solutions to file.

  -- REQUIRED :
  --   The file must be opened for output.  For format consistency,
  --   a successful Write_First must have been executed.

  -- ON ENTRY :
  --   file      file opened for output;
  --   n         number of solutions to be written to file;
  --   cnt       number of solutions already written to file;
  --   sols      list of solutions to be written to file.

  -- ON RETURN :
  --   nbr       number of solutions written to file;
  --   cnt       updated counter for number of solutions written to file.

end DoblDobl_Complex_Solutions_io;
