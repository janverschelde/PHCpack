with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Cell_Stack;                        use Cell_Stack;

package Mixed_Volume is

-- DESCRIPTION :
--   This package contains the main program to compute the mixed volume
--   of a tuple of supports, usually invoked after some preprocessing.

  procedure MixedVol 
               ( nVar,nSpt,CellSize : in integer32;
                 SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 nbCells : out integer32; MCells : out CellStack;
                 MVol : out natural32;
                 multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   Computes the mixed volume of the given lifted supports.

  -- ON ENTRY :
  --   nVar      ambient dimension, number of variables;
  --   nSpt      number of different supports;
  --   CellSize  number of indices in one mixed cell;
  --   SptType   type of each support;
  --   SptIdx    SptIdx[i] indicates the start of i-th support in Spt;
  --   Spt       coordinates of the points in the supports;
  --   lft       lifting of the polynomial system.

  -- ON RETURN :
  --   nbCells   number of mixed cells;
  --   MCells    mixed cells of a regular mixed-cell configuration;
  --   MVol      the mixed volume of the polynomial system;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --             must be used for the Hermite normal form.

  procedure MixedVol_with_Callback
               ( nVar,nSpt,CellSize : in integer32;
                 SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 nbCells : out integer32; MCells : out CellStack;
                 MVol : out natural32;
                 multprec_hermite : in boolean := false;
                 next_cell : access procedure
                   ( idx : Standard_Integer_Vectors.Link_to_Vector ) := null );

  -- DESCRIPTION :
  --   Computes the mixed volume of the given lifted supports.
  --   Each time a new cell is found, next_cell is called.

  -- ON ENTRY :
  --   nVar      ambient dimension, number of variables;
  --   nSpt      number of different supports;
  --   CellSize  number of indices in one mixed cell;
  --   SptType   type of each support;
  --   SptIdx    SptIdx[i] indicates the start of i-th support in Spt;
  --   Spt       coordinates of the points in the supports;
  --   lft       lifting of the polynomial system.

  -- ON RETURN :
  --   nbCells   number of mixed cells;
  --   MCells    mixed cells of a regular mixed-cell configuration;
  --   MVol      the mixed volume of the polynomial system;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --             must be used for the Hermite normal form;
  --   next_cell is a callback procedure, the argument to next_cell are
  --             indices to the points in the cell.

  procedure gcd ( r1,r2 : in integer32; k,l,d : out integer32 );

  -- DESCRIPTION :
  --   Returns in d the gcd of the numbers r1 and r2;
  --   after gcd(r1,r2,k,l,d), d = k*r1 + l*r2 holds on return.

  function cell_size
               ( nSpt : integer32;
                 SptType : Standard_Integer_Vectors.Link_to_Vector )
               return integer32;

  -- DESCRIPTION :
  --   Returns the number of indices it takes to represent a mixed cell.

  -- ON ENTRY :
  --   nSpt      number of different supports;
  --   SptType   #occurrences of each support.

  procedure CellVol
               ( nVar,nSpt : in integer32;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 SptType,Cell : in Standard_Integer_Vectors.Link_to_Vector;
                 Vol : in out natural32;
                 multprec_hermite : in boolean := false );

  -- DESCRIPTION :
  --   Returns the mixed volume of the current cell.

  -- ON ENTRY :
  --   nVar      dimension of the ambient space;
  --   nSpt      number of different supports;
  --   Spt       coordinates of the points in the supports;
  --   Cell      a mixed cell;
  --   Vol       current mixed volume;
  --   multprec_hermite indicates whether multiprecision arithmetic
  --            must be used for the Hermite normal form.

  -- ON RETURN :
  --   Vol       accumulated mixed volume.

  procedure solve_linear_system
               ( n : in integer32;
                 A : in Standard_Floating_Matrices.Matrix;
                 b : in out Standard_Floating_Vectors.Vector;
                 okay : out integer32 );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b.

  -- NOTE :
  --   This routine is used only to compute the inner normal
  --   when writing the mixed-cell configuration to file.
  --
  -- ON ENTRY :   
  --   n         dimension of the linear system;
  --   A         an n-by-n matrix;
  --   b         right hand side vector of size n.

  -- ON RETURN :
  --   b         solution x to A*x = b,
  --   okay      if the integer returned by this function is 1,
  --             otherwise if A is singular, the okay is 0.

  function Inner_Normal
               ( n,r : integer32;
                 m,c : Standard_Integer_Vectors.Link_to_Vector;
                 Spt : Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : Standard_Floating_Vectors.Link_to_Vector )
               return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the inner normal spanned by those points in the supports
  --   determined by the labels.

  -- ON ENTRY :
  --   n         number of variables;
  --   r         number of different supports;
  --   m         type of mixture;
  --   c         labels to the points in Spt which span the mixed cell;
  --   Spt       coordinates of the supports;
  --   lft       values of the lifting on the supports.

  -- ON RETURN :
  --   vector of range 0..n-1, normal to the points in the mixed cell.

  procedure write_cells
               ( file : in string; nVar,nSpt : in integer32;
                 SptType : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nMCells : in integer32; MCells : in out CellStack );

  -- DESCRIPTION :
  --   Writes the regular mixed-cell configuration to file,
  --   in a format suitable for PHCpack.

  -- ON ENTRY :
  --   file      name of the output file for the mixed-cell configuration;
  --   nVar      ambient dimension, number of variables;
  --   nSpt      number of different supports;
  --   SptType   type of each support;
  --   Spt       coordinates of the points in the supports;
  --   lft       lifting of the polynomial system;
  --   CellSize  number of points in a mixed cell;
  --   MCells    mixed cells of a regular mixed-cell configuration.

  -- ON RETURN :
  --   MCells    mixed cells are popped off the stack ...

end Mixed_Volume;
