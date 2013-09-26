with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Floating_Faces_of_Polytope;         use Floating_Faces_of_Polytope;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Floating_Pruning_Methods is

-- DESCRIPTION :
--   This package contains the creators of a regular mixed subdivision,
--   based on the static lifting algorithm, for computing only those cells
--   of a certain type, in particular the mixed cells.
--   There are facilities for computing only the generating cells and for
--   computing only the stable mixed cells.

  generic

    with procedure Process ( mic : in Mixed_Cell; continue : out boolean );

    -- DESCRIPTION :
    --   This procedure will be invoked after each computation of a new cell.
    --   If the parameter continue is set on false, then the computation will
    --   be stopped, otherwise the creation continues.

  procedure Gen1_Create
               ( n : in integer32; mix : in Vector; fa : in Array_of_Faces;
                 lifted : in Array_of_Lists; tol : in double_float;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 mixsub : out Mixed_Subdivision );

  procedure Create
               ( n : in integer32; mix : in Vector; fa : in Array_of_Faces;
                 lifted : in Array_of_Lists; tol : in double_float;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Creates a mixed subdivision with a criterium to check which
  --   face-face combinations can lead to a cell which contributes to
  --   the mixed volume.

  -- ON ENTRY :
  --   n         dimension before lifting;
  --   mix       type of mixture: indicates how many times each polytope
  --             occurs in the supports;
  --   fa        faces of the lower hull of the lifted point sets:
  --              fa(i) contains the mix(i)-faces of conv(lifted(i));
  --   lifted    the lifted points;
  --   tol       tolerance on the precision;

  -- ON RETURN :
  --   nbsucc    number of times a face-face combination has passed the test,
  --             at each level;
  --   nbfail    number of times a face-face combinations has failed to pass
  --             the test, at each level;
  --   mixsub    collection of cells which contribute to the mixed volume.

end Floating_Pruning_Methods;
