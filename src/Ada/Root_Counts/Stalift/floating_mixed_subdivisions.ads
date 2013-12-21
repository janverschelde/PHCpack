with Generic_Lists;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;

package Floating_Mixed_Subdivisions is

-- DESCRIPTION :
--   This package enables working with regular mixed subdivisions.

-- DATA STRUCTURES :
--   A mixed cell can be stored in coordinate form, or using labels
--   with respect to a lifted configurations of points.
--   The coordinate form is self contained and perhaps most natural,
--   while the use of labels is more compact and better for global use.

-- coordinate representation :

  type Mixed_Subdivision;              -- list of mixed cells
  type Link_to_Mixed_Subdivision is access Mixed_Subdivision;

  type Mixed_Cell is record
    nor : Standard_Floating_Vectors.Link_to_Vector; -- inner normal
    pts : Link_to_Array_of_Lists;      -- points that span the cell
    sub : Link_to_Mixed_Subdivision;   -- subdivision of the cell
  end record;

  package Lists_of_Mixed_Cells is new Generic_Lists(Mixed_Cell);
  type Mixed_Subdivision is new Lists_of_Mixed_Cells.List;

-- labeled representation :

  type Mixed_Sublabeling;
  type Link_to_Mixed_Sublabeling is access Mixed_Sublabeling;

  type Mixed_Labels is record
    nor : Standard_Floating_Vectors.Link_to_Vector;   -- inner normal
    lab : Standard_Integer_VecVecs.Link_to_VecVec;    -- labels to points
    sub : Link_to_Mixed_Sublabeling;                  -- subdivision of cell
  end record;

  package Lists_of_Mixed_Labels is new Generic_Lists(Mixed_Labels);

  type Mixed_Sublabeling is record
    pts : Standard_Floating_VecVecs.Link_to_Array_of_VecVecs; -- lifted points
    cells : Lists_of_Mixed_Labels.List;
  end record;

-- CREATORS :

  function Compute_Inner_Normal
              ( n : integer32; pts : Array_of_Lists )
              return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the inner normal to the points in pts that span a mixed cell.
  --   The length of the points in pts equals n.

  function Recompute_Inner_Normal
              ( n : integer32; b : double_float; pts : Array_of_Lists )
              return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Recomputes the inner normal for the points in pts of length n, 
  --   multiplying by 10 the lifting value for the artificial origins,
  --   lifted with value equal to b.

  function Create ( pts : Array_of_Lists; 
                    nor : Standard_Floating_Vectors.Vector;
                    tol : double_float ) return Mixed_Cell;
  function Create ( pts : Array_of_Lists; nors : List; tol : double_float )
                  return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Creates the mixed cell(s) of those points whose inner product
  --   with the given normal(s) is minimal.
  --   The parameter tol is the tolerance on the precision.

  function Create ( pts : Array_of_Lists; mixsub : Mixed_Subdivision;
                    tol : double_float ) return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Takes the normals of the cells in the given mixed subdivision
  --   and creates the mixed cells by selecting the points whose inner
  --   product with the normals are minimal.

  procedure Update ( pts : in Array_of_Lists;
                     nor : in Standard_Floating_Vectors.Vector;
                     mixsub,mixsub_last : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Given a tuple of points and a normal, the mixed subdivision will
  --   be updated: either an existing cell will get additional points,
  --   if the normal already occurs in the subdivision, or otherwise,
  --   a new cell will be created and appended to the mixed subdivision.

-- CONVERTORS :

  function Create ( s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Floating_VecVecs.Array_of_VecVecs;

  -- DESCRIPTION :
  --   Turns the array of lists into an array of vectors of vectors,
  --   making a deep copy of the data.

  function Position ( lifpts : Standard_Floating_VecVecs.VecVec;
                      point : Standard_Floating_Vectors.Vector )
                    return integer32;

  -- DESCRIPTION :
  --   Returns the position of the point in the vector of lifted points,
  --   if the point occurs in lifpts, otherwise zero is returned.

  function Create_Labels
              ( pts : Standard_Floating_VecVecs.Array_of_VecVecs;
                mic : Mixed_Cell ) return Mixed_Labels;

  -- DESCRIPTION :
  --   Returns a labeled representation of the mixed cell,
  --   using the positions of the points in a lifted point configuration.

  -- ON ENTRY :
  --   pts      a lifted point configuration;
  --   mic      a mixed cell.

  function Create_Labeled_Subdivision
              ( r : integer32; sub : Mixed_Subdivision )
              return Mixed_Sublabeling;

  -- DESCRIPTION :
  --   Converts a coordinate representation of a mixed-cell configuration
  --   into a labeled representation.

  function Create_Coordinates
              ( pts : Standard_Floating_VecVecs.Array_of_VecVecs;
                mlb : Mixed_Labels ) return Mixed_Cell;

  -- DESCRIPTION :
  --   Returns the coordinate representation of the mixed cell,
  --   using the lifting points in pts.

  function Create_Coordinate_Subdivision
              ( r : integer32; sub : Mixed_Sublabeling )
              return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Converts a labeled representation of a mixed-cell configuration
  --   into a coordinate representation.

-- CONSTRUCTORS :

  procedure Copy ( mic1 : in Mixed_Cell; mic2 : in out Mixed_Cell );
  procedure Copy ( mixsub1 : in Mixed_Subdivision; 
                   mixsub2 : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Makes a deep copy of the cells and the subdivisions.

  procedure Append_Diff ( first,last : in out Mixed_Subdivision;
                          mic : in Mixed_Cell );

  -- DESCRIPTION :
  --   Appends a mixed cell to the list of cells first, where
  --   last points to the last element of the list first.
  --   The suffix _Diff means that only when the cell does not already
  --   belong to the list first, it will be appended.

  procedure Concat_Diff ( first,last : in out Mixed_Subdivision;
                          mixsub : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Concatenates all cells in mixsub to the list of cells first,
  --   last is a pointer to the last cell in first.
  --   The suffix _Diff means that only when those cells that do not already
  --   belong to the list first will be appended.

  procedure Construct ( mixsub : in Mixed_Subdivision;
                        first : in out Mixed_Subdivision );
  procedure Construct_Diff ( mixsub : in Mixed_Subdivision;
                             first : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Constructs all cells in the mixed subdivision to the front of
  --   the list first.
  --   The suffix _Diff means that only when those cells that do not already
  --   belong to the list first will be constructed to first.

-- SELECTORS :

  function Equal ( mic1,mic2 : Mixed_Cell ) return boolean;
  function Equal ( mixsub1,mixsub2 : Mixed_Subdivision ) return boolean;
  function Equal ( mixsub1,mixsub2 : Link_to_Mixed_Subdivision )
                 return boolean;

  -- DESCRIPTION :
  --   Returns true when two mixed cells and mixed subdivisions are equal.

  function Is_In ( mixsub : Mixed_Subdivision;
                   normal : Standard_Floating_Vectors.Vector ) return boolean;
  function Is_In ( mixsub : Mixed_Subdivision; mic : Mixed_Cell )
                 return boolean;

  -- DESCRIPTION :
  --   Returns true if normal or cell belongs to the mixed subdivision,
  --   otherwise false is returned.  When the whole mixed cell is given,
  --   then not only the normal, but also the points will be checked.

  function Is_Original ( mic : Mixed_Cell; b : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if no component of the mixed cell contains
  --   an artificial origin with lifting value equal to b.

  function Is_Stable
               ( nor : Standard_Floating_Vectors.Vector;
                 b : double_float; pts : Array_of_Lists ) return boolean;

  -- DESCRIPTION :
  --   Determines the stability of a mixed cell by recomputing its inner
  --   normal, after lifting the artificial origins to height 10*b.

  -- ON ENTRY :
  --   nor       inner normal to the points in pts;
  --   b         lifting bound used for the stable mixed volume;
  --   pts       points that span the mixed cell.

  -- ON RETURN :
  --   true      if after recomputing the inner normal lifting the
  --             artificial origins to 10*b, no components of the
  --             inner normal became more negative,
  --   false     otherwise.

  function Zero_Type
               ( nor : Standard_Floating_Vectors.Vector;
                 b : double_float; pts : Array_of_Lists )
               return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   By recomputing its inner normal after relifting the artificial
  --   origins, decides which components of the solutions are zero or not.

  -- ON ENTRY :
  --   nor       inner normal to the points in pts;
  --   b         lifting bound used for the stable mixed volume;
  --   pts       points that span the mixed cell.

  -- ON RETURN :
  --   Vector of range nor'first..nor'last-1 of 0, +1, or -1.
  --   if the k-th component is zero then the k-th component of the
  --   solution will be zero, otherwise the k-th component is nonzero.
  --   If the k-th component of the solution will diverge,
  --   then the corresponding component of the vector on return is -1.

  function Is_Stable ( mic : Mixed_Cell; b : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the mixed cell is stable.
  --   A mixed cell is stable if either it contains no artificial origins,
  --   or otherwise no components of its inner normal are less than -1.

  procedure Split_Original_Cells
               ( mixsub : in Mixed_Subdivision; b : in double_float;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 orgcnt,stbcnt : out natural32 );

  -- DESCRIPTION :
  --   Splits the mixed subdivision into original cells (without artificial
  --   origin) and stable mixed cells, using the lifting bound b.
  --   Spurious cells are not in any list on output.

  -- ON ENTRY :
  --   mixsub    a regular mixed-cell configuration;
  --   b         lifting bound for the artificial origins.

  -- ON RETURN :
  --   orgmcc    cells that do not contain an artificial origin;
  --   stbmcc    extra cells with an artificial origin that are stable
  --             and thus contribute to solutions with zero component.
  --   orgcnt    counts the original cells in orgmcc;
  --   stbcnt    counts the extra stable cells.

-- DESTRUCTORS :

  procedure Deep_Clear ( mic : in out Mixed_Cell );
  procedure Deep_Clear ( mixsub : in out Mixed_Subdivision );
  procedure Deep_Clear ( mixsub : in out Link_to_Mixed_Subdivision );
  procedure Shallow_Clear ( mic : in out Mixed_Cell );
  procedure Shallow_Clear ( mixsub : in out Mixed_Subdivision );
  procedure Shallow_Clear ( mixsub : in out Link_to_Mixed_Subdivision );

  -- DESCRIPTION :
  --   The memory space allocated will be freed.
  --   A shallow clear only destroys the list structures,
  --   while with a deep clear, also the contents of the lists are freed.

  procedure Clear ( mlb : in out Mixed_Labels );
  procedure Clear ( mixsub : in out Mixed_Sublabeling );
  procedure Clear ( mixsub : in out Link_to_Mixed_Sublabeling );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the labeled representations.
  --   Since there is no redundancy in these labeled representations,
  --   the Clear operations are always deep.

end Floating_Mixed_Subdivisions;
