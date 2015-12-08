with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Integer_VecVecs;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with Cell_Stack;                        use Cell_Stack;

package MixedVol_Algorithm is

-- DESCRIPTION :
--   This package provides an interface to the Ada translation of
--   ACM TOMS Algorithm 846: MixedVol to compute mixed volumes.

  function Flatten ( n,m : integer32; v : Standard_Integer_VecVecs.VecVec )
                   return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns an integer vector of size n*m with the coordinates
  --   of the vectors in v.

  procedure Flatten_Supports
               ( s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 ind : out Standard_Integer_Vectors.Vector;
                 sup : out Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Flattens the array of lists s into one vector.

  -- ON ENTRY :
  --   s         n supports of a polynomial system.

  -- ON RETURN :
  --   ind       indices to the start of each support: ind(i) points to
  --             the start of the i-th support s(i) in the vector sup;
  --   sup       contains points in the support.

  procedure Extract_Supports
               ( n : in integer32; p : in Poly_Sys; m : out integer32;
                 ind,cnt : out Standard_Integer_Vectors.Vector;
                 sup : out Standard_Integer_Vectors.Link_to_Vector );
  procedure Extract_Supports
               ( n : in integer32; p : in Laur_Sys; m : out integer32;
                 ind,cnt : out Standard_Integer_Vectors.Vector;
                 sup : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system,
  --   and converts them into a suitable format.

  -- ON ENTRY :
  --   n         number of variables and equations in p;
  --   p         a polynomial system.

  -- ON RETURN :
  --   m         total number of points in the supports of p;
  --   ind       ind(i) marks the beginning of the i-th support;
  --   cnt       cnt(i) counts the number of elements in the i-th support;
  --   sup       vector of range 1..n*m with coordinates of all supports.

  procedure Extract_Supports
               ( n : in integer32; 
                 s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 m : out integer32;
                 ind,cnt : out Standard_Integer_Vectors.Vector;
                 sup : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Brings the supports of a polynomial systems into a format
  --   suitable for the mixed volume calculation algorithm.

  -- ON ENTRY :
  --   n         number of variables and number of sets in s;
  --   s         supports of a polynomial system.

  -- ON RETURN :
  --   m         total number of points in the supports of p;
  --   ind       ind(i) marks the beginning of the i-th support;
  --   cnt       cnt(i) counts the number of elements in the i-th support;
  --   sup       vector of range 1..n*m with coordinates of all supports.

  procedure Write_Supports 
               ( nSpt : in integer32;
                 Idx : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Writes the supports to screen.

  -- ON ENTRY :
  --   nSpt      number of supports;
  --   Idx       index set to the supports;
  --   Spt       coordinates of the supports.

  procedure Add_Artificial_Origins 
               ( nVar,nSpt : in integer32;
                 Idx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 nbadd : out integer32;
                 added : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Add artificial origins to the supports for the stable mixed volume.

  -- ON ENTRY :
  --   nVar      number of variables;
  --   nSpt      number of supports;
  --   Idx       index set for each support;
  --   Spt       coordinates of the supports.

  -- ON RETURN :
  --   Idx       updated index set for each support;
  --   Spt       updated coordinates of the supports;
  --   nbadd     number of artificial origins added;
  --   added     indices where artificial origins are in Spt,
  --             may be just the zero pointer if nbadd = 0.

  procedure mv ( nVar,nPts : in integer32;
                 ind,cnt,sup : in Standard_Integer_Vectors.Vector;
                 stlb : in double_float; nSpt : out integer32;
                 SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : out Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nbCells : out integer32; cells : out CellStack;
                 mixvol : out natural32;
                 multprec_hermite : in boolean := false );

  -- DESCRITPION :
  --   Computes the mixed volume of the polytopes spanned by given supports.
  --   A regular mixed-cell configuration is available in CellStack.

  -- ON ENTRY :
  --   nVar      ambient dimension, length of the vectors in supports;
  --   nPts      total number of points in the supports;
  --   ind       ind(i) is the start of the i-th support;
  --   cnt       cnt(i) counts the length of the i-th support;
  --   sup       coordinates of the points in the supports;
  --   stlb      lifting bound for stable mixed volumes,
  --             equals 0.0 if no stable mixed volumes are needed.
  --   multprec_hermite indicates whether multiprecision arithmetic
  --             must be used for the Hermite normal form.

  -- ON RETURN :
  --   nSpt      number of different supports;
  --   SptType   type of each support;
  --   perm      permutation of the original supports;
  --   VtxIdx    index vector to the vertex set;
  --   Vtx       vertices of the supports;
  --   lft       lifting values for the vertex points;
  --   CellSize  is the size of a mixed cell;
  --   nbCells   total number of mixed cells;
  --   cells     stack of mixed cells;
  --   mixvol    mixed volume of the polytopes spanned by the supports,
  --             note that this is the total mixed volume, if stlb /= 0.0.

  procedure mv_upto_pre4mv
               ( nVar,nPts : in integer32;
                 ind,cnt,sup : in Standard_Integer_Vectors.Vector;
                 nSpt : out integer32;
                 SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : out Standard_integer_Vectors.Link_to_Vector;
                 Spt : out Standard_Integer_VecVecs.Link_to_VecVec; 
                 NuIdx2OldIdx : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Executes the same code as in mv_with_callback upto the 'pre4mv'
  --   procedure which computes the number of different supports,
  --   the type of mixture, and a permutation of the original supports. 
  --   This information is needed to process the cells returned by
  --   the callback procedure next_cell in mv_with_callback.

  -- ON ENTRY :
  --   nVar      ambient dimension, length of the vectors in supports;
  --   nPts      total number of points in the supports;
  --   ind       ind(i) is the start of the i-th support;
  --   cnt       cnt(i) counts the length of the i-th support;
  --   sup       coordinates of the points in the supports.

  -- ON RETURN :
  --   nSpt      number of different supports;
  --   SptType   type of each support;
  --   perm      permutation of the original supports;
  --   VtxIdx    index vector to the vertex set;
  --   Vtx       vertices of the supports.
  --   SptIdx    index vector to the support set;
  --   Spt       points in the support set;
  --   NuIdx2OldIdx relates new to old indices.

  procedure mv_lift
               ( nVar : in integer32;
                 stlb : in double_float; nSpt : in integer32;
                 VtxIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : out Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Performs a random lifting on the vertex set and if stlb /= 0,
  --   then artificial origins will be added.

  -- ON ENTRY :
  --   nVar      ambient dimension, length of the vectors in supports;
  --   stlb      lifting bound for stable mixed volumes,
  --             equals 0.0 if no stable mixed volumes are needed.
  --   nSpt      number of different supports;
  --   VtxIdx    index vector to the vertex set;
  --   Vtx       vertices of the supports.

  -- ON RETURN :
  --   VtxIdx    index vector to the vertex set,
  --             with possible addition of artificial origins;
  --   Vtx       vertices of the supports, 
  --             with possible addition of artificial origins;
  --   lft       lifting values for the vertex points.
 
  procedure mv_with_callback
               ( nVar,nPts : in integer32;
                 ind,cnt,sup : in Standard_Integer_Vectors.Vector;
                 stlb : in double_float; nSpt : out integer32;
                 SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : out Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nbCells : out integer32; cells : out CellStack;
                 mixvol : out natural32;
                 multprec_hermite : in boolean := false;
                 next_cell : access procedure
                   ( idx : Standard_Integer_Vectors.Link_to_Vector ) := null );

  -- DESCRITPION :
  --   Computes the mixed volume of the polytopes spanned by given supports.
  --   A regular mixed-cell configuration is available in CellStack.
  --   Each time a new cell is found, the procedure next_cell is called.

  -- ON ENTRY :
  --   nVar      ambient dimension, length of the vectors in supports;
  --   nPts      total number of points in the supports;
  --   ind       ind(i) is the start of the i-th support;
  --   cnt       cnt(i) counts the length of the i-th support;
  --   sup       coordinates of the points in the supports;
  --   stlb      lifting bound for stable mixed volumes,
  --             equals 0.0 if no stable mixed volumes are needed.
  --   multprec_hermite indicates whether multiprecision arithmetic
  --             must be used for the Hermite normal form;
  --   next_cell is a callback procedure, the input argument idx are the
  --             indices to the points that span the new cell.

  -- ON RETURN :
  --   nSpt      number of different supports;
  --   SptType   type of each support;
  --   perm      permutation of the original supports;
  --   VtxIdx    index vector to the vertex set;
  --   Vtx       vertices of the supports;
  --   lft       lifting values for the vertex points;
  --   CellSize  is the size of a mixed cell;
  --   nbCells   total number of mixed cells;
  --   cells     stack of mixed cells;
  --   mixvol    mixed volume of the polytopes spanned by the supports,
  --             note that this is the total mixed volume, if stlb /= 0.0.

  procedure uliftmv
               ( nVar,nPts : in integer32;
                 ind,cnt,sup : in Standard_Integer_Vectors.Vector;
                 stlb : in double_float; nSpt : out integer32;
                 SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : out Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nbCells : out integer32; cells : out CellStack;
                 mixvol : out natural32 );

  -- DESCRITPION :
  --   Computes the mixed volume of the polytopes spanned by given supports.
  --   A regular mixed-cell configuration is available in CellStack.
  --   In this procedure the user is prompted for the lifting values.

  -- ON ENTRY :
  --   nVar      ambient dimension, length of the vectors in supports;
  --   nPts      total number of points in the supports;
  --   ind       ind(i) is the start of the i-th support;
  --   cnt       cnt(i) counts the length of the i-th support;
  --   sup       coordinates of the points in the supports;
  --   stlb      lifting bound for stable mixed volumes,
  --             equals 0.0 if no stable mixed volumes are needed.

  -- ON RETURN :
  --   nSpt      number of different supports;
  --   SptType   type of each support;
  --   perm      permutation of the original supports;
  --   VtxIdx    index vector to the vertex set;
  --   Vtx       vertices of the supports;
  --   lft       lifting values for the vertex points;
  --   CellSize  is the size of a mixed cell;
  --   nbCells   total number of mixed cells;
  --   cells     stack of mixed cells;
  --   mixvol    mixed volume of the polytopes spanned by the supports,
  --             note that this is the total mixed volume, if stlb /= 0.0.

  function Supports_of_Mixed_Cell
               ( nVar,nSpt : integer32;
                 SptType,labels : Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : Standard_Floating_Vectors.Link_to_Vector )
               return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the supports of the mixed cell defined by the labels.

  -- ON ENTRY :
  --   nVar      number of variables, dimension before the lifting;
  --   nSpt      number of different supports;
  --   SptType   type of mixture of the supports;
  --   labels    indices to the point in the mixed cell;
  --   Vtx       complete vertex set of all supports;
  --   lft       lifting values for all points in Vtx.

  function Supports_of_Mixed_Cell
               ( nVar,nSpt : integer32;
                 SptType,perm : Standard_Integer_Vectors.Link_to_Vector;
                 labels : Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : Standard_Floating_Vectors.Link_to_Vector )
               return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the supports of the mixed cell defined by the labels,
  --   taking the permutation into account, i.e.: the supports of the
  --   cells follow the original order of the supports.

  -- ON ENTRY :
  --   nVar      number of variables, dimension before the lifting;
  --   nSpt      number of different supports;
  --   SptType   type of mixture of the supports;
  --   perm      permutation of the original support sets;
  --   labels    indices to the point in the mixed cell;
  --   Vtx       complete vertex set of all supports;
  --   lft       lifting values for all points in Vtx.

  function Labels_to_Mixed_Cell
               ( nVar,nSpt : in integer32;
                 SptType : in Standard_Integer_Vectors.Link_to_Vector;
                 labels : in Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector )
               return Mixed_Cell;
  function Labels_to_Mixed_Cell
               ( nVar,nSpt : in integer32;
                 SptType,perm : in Standard_Integer_Vectors.Link_to_Vector;
                 labels : in Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector )
               return Mixed_Cell;

  -- DESCRIPTION :
  --   Returns the mixed cell defined by the labels of indices to the points.

  -- ON ENTRY :
  --   nVar      number of variables, dimension before the lifting;
  --   nSpt      number of different supports;
  --   SptType   type of mixture of the supports;
  --   perm      permutation of the original support sets;
  --   labels    indices to the point in the mixed cell;
  --   Vtx       complete vertex set of all supports;
  --   lft       lifting values for all points in Vtx.

  procedure Create_Mixed_Cell_Configuration
               ( nVar,nSpt,CellSize,nbCells : in integer32;
                 SptType : in Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 cells : in out CellStack; sub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Creates a regular mixed-cell configuration from the stack of cells.

  -- ON ENTRY :
  --   nVar      number of variables, dimension before the lifting;
  --   nSpt      number of different supports;
  --   CellSize  is the number of labels in a mixed cells;
  --   nbCells   total number of mixed cells;
  --   SptType   type of mixture of the supports;
  --   Vtx       complete vertex set of all supports;
  --   lft       lifting values for all points in Vtx;
  --   cells     stack of labels to coordinates of mixed cells.

  -- ON RETURN :
  --   cells     elements have been popped off the stack;
  --   sub       regular mixed-cell configuration.

  procedure Create_Mixed_Cell_Configuration
               ( nVar,nSpt,CellSize,nbCells : in integer32;
                 SptType,perm : in Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 cells : in out CellStack; sub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Creates a regular mixed-cell configuration from the stack of cells,
  --   taking the permutation of the original supports into account,
  --   i.e.: writing supports of the cells in the original order.

  -- ON ENTRY :
  --   nVar      number of variables, dimension before the lifting;
  --   nSpt      number of different supports;
  --   CellSize  is the number of labels in a mixed cells;
  --   nbCells   total number of mixed cells;
  --   SptType   type of mixture of the supports;
  --   perm      permutation of the original supports;
  --   Vtx       complete vertex set of all supports;
  --   lft       lifting values for all points in Vtx;
  --   cells     stack of labels to coordinates of mixed cells.

  -- ON RETURN :
  --   cells     elements have been popped off the stack;
  --   sub       regular mixed-cell configuration.

end MixedVol_Algorithm;
