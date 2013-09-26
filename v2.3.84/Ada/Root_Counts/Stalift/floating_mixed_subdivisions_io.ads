with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;  

package Floating_Mixed_Subdivisions_io is

-- DESCRIPTION :
--   This package provides some routines for i/o of mixed subdivisions,
--   induced by a floating-point lifting.

-- I/O of coordinate representations :

  procedure get ( n,m : in natural32; mic : out Mixed_Cell );
  procedure get ( file : in file_type;
                  n,m : in natural32; mic : out Mixed_Cell );

  -- DESCRIPTION :
  --   Reads the normal and for each list of points, the length
  --   of the list and the list itself from standard input or from file.

  procedure get ( n,m : out natural32;
                  mixed_type : out Standard_Integer_Vectors.Link_to_Vector;
                  mixsub : out Mixed_Subdivision );
  procedure get ( file : in file_type; n,m : out natural32;
                  mixed_type : out Standard_Integer_Vectors.Link_to_Vector;
                  mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Reads the dimension, the number of different supports,
  --   the type of mixture, the number of mixed cells and the mixed
  --   subdivision from standard input or from file.

  procedure put ( lifvec : in Standard_Floating_Vectors.Vector );
  procedure put ( file : in file_type;
                  lifvec : in Standard_Floating_Vectors.Vector );
  procedure put ( lifsup : in List );
  procedure put ( file : in file_type; lifsup : in List );
  procedure put ( lifsup : in Array_of_Lists );
  procedure put ( file : in file_type; lifsup : in Array_of_Lists );

  -- DESCRIPTION :
  --   Writes the lifted vectors on file or on standard output.
  --   The format uses fixed point format for the integer entries.

  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mic : in Mixed_Cell );
  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mic : in out Mixed_Cell; mv : out natural32;
                  multprec_hermite : in boolean := false );
  procedure put ( file : in file_type;
                  n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mic : in Mixed_Cell );
  procedure put ( file : in file_type;
                  n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mic : in out Mixed_Cell; mv : out natural32;
                  multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   Puts the normal, the length of each point list and the
  --   points belonging to the cell on standard output or on file.
  --   More text banners are provided when the parameter `mv' is supplied.
  --   This `mv' contains the mixed volume of the cell on return.
  --   When the mixed volume is computed, eventually the cell is refined.

  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision );
  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural32;
                  multprec_hermite : in boolean := false );
  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision );
  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural32;
                  multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   Puts the dimension, the type of mixture, the number of mixed
  --   cells and the mixed subdivision on file or on standard output.
  --   More text banners are provided when the parameter `mv' is supplied.
  --   This `mv' contains the mixed volume of the cell on return.
  --   When the mixed volume is computed, eventually the cells are refined.
  --   If multprec_hermite, then the multiprecision Hermite normal form
  --   is used to compute determinants.

  procedure put ( file : in file_type; n : in natural32; b : in double_float;
                  mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; 
                  mv,smv,tmv : out natural32;
                  multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   In addition to writing the mixed subdivision and mixed volume to file,
  --   this routine also computes the stable mixed volume, where b is the
  --   lifting bound used for the artificial origins.

-- INCREMENTAL READ and WRITE (for huge mixed cell configurations) :

  procedure Read_Dimensions
               ( file : in file_type; n,r,m : out natural32;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Reads the dimensions of a mixed cell configuration from file.

  -- ON ENTRY :
  --   file      must be opened for input.

  -- ON RETURN :
  --   n         the ambient dimension before lifting the points;
  --   r         the number of different supports;
  --   m         the number of mixed cells on file;
  --   mix       the type of mixture, #occurrences for each support;
  --   fail      false if no exceptions occurred during reading,
  --             true if some exception occurred.

  procedure Write_Dimensions
               ( file : in file_type; n,r,m : in natural32;
                 mix : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the dimensions of a mixed-cell configuration to file.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   n         the ambient dimension before lifting the points;
  --   r         the number of different supports;
  --   m         the number of mixed cells on file;
  --   mix       the type of mixture, #occurrences for each support.
               
  procedure Read_Next
               ( file : in file_type; n,r : in natural32;
                 mic : out Mixed_Cell; fail : out boolean );

  -- DESCRIPTION :
  --   Reads the next mixed-cell from file.

  -- ON ENTRY :
  --   file      must be opened for input;
  --   n         the ambient dimension, before the lifting;
  --   r         number of different supports.

  -- ON RETURN :
  --   mic       a mixed cell if not fail;
  --   fail      true if some exception occurred during reading;
  --             false if the reading of a mixed cell was successful.

-- I/O of labeled representations :

  procedure get ( file : in file_type; n,r : in natural32;
                  mlb : out Mixed_Labels );

  -- DESCRIPTION :
  --   Reads a labeled representation of a mixed cell in n-space,
  --   for r different support from file.

  procedure get ( file : in file_type; n,r : out natural32;
                  mix : out Standard_Integer_Vectors.Link_to_Vector;
                  sub : out Mixed_Sublabeling );

  -- DESCRIPTION :
  --   Reads a labeled representation of a mixed-cell configuration
  --   from file, which must be opened for input.

  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  mlb : in Mixed_Labels );

  -- DESCRIPTION :
  --   Writes the labeled representation of a mixed cell to file.
  --   The type of mixture is needed in case the cell is not fine,
  --   when its refinement needs to be written to file as well.

  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  sub : in Mixed_Sublabeling );

  -- DESCRIPTION :
  --   Writes the labeled representation of a regular mixed-cell configuration
  --   to file, which must be opened for output.

-- INCREMENTAL READ/WRITE for labeled representations :

  procedure Read_Dimensions
               ( file : in file_type; n,r : out natural32;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Reads the dimensions of a mixed cell configuration from file.

  -- ON ENTRY :
  --   file      must be opened for input.

  -- ON RETURN :
  --   n         the ambient dimension before lifting the points;
  --   r         the number of different supports;
  --   mix       the type of mixture, #occurrences for each support;
  --   fail      false if no exceptions occurred during reading,
  --             true if some exception occurred.

  procedure Write_Dimensions
               ( file : in file_type; n,r : in natural32;
                 mix : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the dimensions of a mixed-cell configuration to file.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   n         the ambient dimension before lifting the points;
  --   r         the number of different supports;
  --   mix       the type of mixture, #occurrences for each support.

  procedure Read_Lifted_Supports
               ( file : in file_type; n,r : in natural32;
                 ls : out Standard_Floating_VecVecs.Array_of_VecVecs;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Reads r lifted supports (vectors of length n+1) from file.

  -- REQUIRED :
  --   The file is open for input and at the proper place;
  --   ls'range = 1..r.

  -- ON ENTRY :
  --   file      input file with at the start r lifted point sets;
  --   n         ambient dimension, before the lifting of the points;
  --   r         number of different supports.

  -- ON RETURN :
  --   ls        lifted supports, must be array of range 1..r;
  --   fail      true if some exception occurred during reading.

  procedure Write_Lifted_Supports
               ( file : in file_type;
                 ls : in Standard_Floating_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Writes the lifted supports to file.

  -- REQUIRED :
  --   The file is open for input.

  -- ON ENTRY :
  --   file      output file to write the lifted point sets;
  --   ls        lifted point sets.

  procedure Read_Next
               ( file : in file_type; n,r,k : in natural32;
                 mlb : out Mixed_Labels; fail : out boolean );

  -- DESCRIPTION :
  --   Reads a mixed cell in labeled representation from file.
  --   This is supposed to be the k-th mixed cell.
  --   The difference with get(file,n,r,mlb) is the suppression of
  --   the raising of an exception.  Instead, when the cell cannot
  --   be read, an error message is printed and fail is true on return.

end Floating_Mixed_Subdivisions_io;
