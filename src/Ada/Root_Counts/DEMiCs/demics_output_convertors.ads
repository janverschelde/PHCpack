with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package DEMiCs_Output_Convertors is

-- DESCRIPTION :
--   The output of DEMiCs is converted into a mixed subdivision.

  function Apply_Lifting
              ( sup : Lists_of_Integer_Vectors.List;
                lif : Standard_Floating_Vectors.Vector )
              return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION :
  --   Given in sup is the support of a polynomial
  --   and in lif are the corresponding lifting values.
  --   Returns the support in sup, with the lifting applied
  --   as defined by the lifting values in lif.

  -- REQUIRED : Length_Of(sup) = lif'last.

  function Apply_Lifting
              ( sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : Standard_Floating_VecVecs.VecVec )
              return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Given in sup are the supports of a polynomial system
  --   and in lif are the corresponding lifting values.
  --   Returns the supports in sup, with the lifting applied
  --   as defined by the lifting values in lif.
  --   Assumes the supports are fully mixed,
  --   i.e.: every support is distinct from every other support.

  -- REQUIRED :
  --   for k in lif'range: length_of(sup(k)) = lif(k)'last.

  function Apply_Lifting
              ( mix : Standard_Integer_Vectors.Vector;
                sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : Standard_Floating_VecVecs.VecVec )
              return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Given in sup are the supports of a polynomial system,
  --   sorted according to the type of mixture in mix,
  --   and in lif are the corresponding lifting values.
  --   Returns the supports in sup, with the lifting applied
  --   as defined by the lifting values in lif.
  --   Assumes the supports are fully mixed,
  --   i.e.: every support is distinct from every other support.

  -- REQUIRED :
  --   for k in lif'range: length_of(sup(i)) = lif(k)'last,
  --   where i is computed corresponding the type of mixture.

  function Minimum
              ( lifpts : Lists_of_Floating_Vectors.List;
                normal : Standard_Floating_Vectors.Vector )
              return double_float;

  -- DESCRIPTION :
  --   Returns the minimal inner product of the vector in normal
  --   with the lifted points in lifpts.

  function Arguments_of_Minima
              ( lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                normal : Standard_Floating_Vectors.Vector )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the indices in lifsup where the minimum inner product
  --   with the normal is attained.

  function Sort_Labels
              ( labels : Standard_Integer_Vectors.Vector )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Given a vector of labels, the vector on returns contains
  --   the labels where each pair of indices is sorted increasingly.

  procedure Test_Inner_Normal
              ( lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                normal : in Standard_Floating_Vectors.Vector;
                labels : in Standard_Integer_Vectors.Vector;
                fail : out boolean );

  -- DESCRIPTION :
  --   Tests whether the minimum of normal with the points in lifsup
  --   is attained at the indices in labels.
  --   Output is written to screen.
  
  -- ON ENTRY :
  --   lifsup   lifted supports;
  --   normal   computed inner normal in Make_Cell;
  --   labels   labels to the points that span the mixed cell.

  -- ON RETURN :
  --   fail     true if the labels disagree with where the minima are
  --            attained, false otherwise.

  function Make_Mixed_Cell
              ( dim : integer32;
                lbl : Standard_Integer_Vectors.Vector;
                lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                verbose : boolean := true )
              return Mixed_Cell;

  -- DESCRIPTION :
  --   Returns the mixed cell corresponding to the labeled points.

  -- REQUIRED : lifsup'last = dim,
  --   i.e.: all supports are distinct from each other.

  -- ON ENTRY :
  --   dim      dimension before lifting;
  --   lbl      labels to the points that span the mixed cell;
  --   lifsup   lifted supports;
  --   verbose  flag to indicate whether the checks with the inner normal
  --            needs to happen.

  function Make_Mixed_Cell
              ( dim : integer32;
                mix,lbl : Standard_Integer_Vectors.Vector;
                lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                verbose : boolean := true )
              return Mixed_Cell;

  -- DESCRIPTION :
  --   Returns the mixed cell corresponding to the labeled points.

  -- ON ENTRY :
  --   dim      dimension before lifting;
  --   mix      type of mixture;
  --   lbl      labels to the points that span the mixed cell;
  --   lifsup   lifted supports;
  --   verbose  flag to indicate whether the checks with the inner normal
  --            needs to happen.

  function Make_Mixed_Cells
              ( dim : integer32;
                labels : Lists_of_Integer_Vectors.List;
                lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                verbose : boolean := true )
              return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the mixed cells corresponding to the labeled points.

  -- REQUIRED : lifsup'last = dim,
  --   i.e.: all supports are distinct from each other.

  -- ON ENTRY :
  --   dim      dimension before lifting;
  --   labels   labels to the points that span the mixed cell;
  --   lifsup   lifted supports;
  --   verbose  if verbose, then the inner normals are checked
  --            and the outcome of each check is written to screen.

  function Make_Mixed_Cells
              ( dim : integer32;
                mix : Standard_Integer_Vectors.Vector;
                labels : Lists_of_Integer_Vectors.List;
                lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                verbose : boolean := true )
              return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the mixed cells corresponding to the labeled points.

  -- ON ENTRY :
  --   dim      dimension before lifting;
  --   mix      type of mixture;
  --   labels   labels to the points that span the mixed cell;
  --   lifsup   lifted supports;
  --   verbose  if verbose, then the inner normals are checked
  --            and the outcome of each check is written to screen.

end DEMiCs_Output_Convertors;
