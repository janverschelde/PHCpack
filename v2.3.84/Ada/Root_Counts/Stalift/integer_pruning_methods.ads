with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Integer_Pruning_Methods is

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

  procedure Gen1_Create_CS 
               ( n : in integer32; mix : in Vector; fa : in Array_of_Faces;
                 lifted : in Array_of_Lists;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 mixsub : out Mixed_Subdivision );

  procedure Create_CS 
               ( n : in integer32; mix : in Vector; fa : in Array_of_Faces;
                 lifted : in Array_of_Lists;
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

  -- ON RETURN :
  --   nbsucc    number of times a face-face combination has passed the test,
  --             at each level;
  --   nbfail    number of times a face-face combinations has failed to pass
  --             the test, at each level;
  --   mixsub    collection of cells which contribute to the mixed volume.

  procedure New_Create_CS
               ( n : in integer32; mix : in Vector; fa : in Array_of_Faces;
                 lifted : in Array_of_Lists;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Does the same at the Create_CS, only with a new pruning method.

  procedure Create_CCS
               ( n : in integer32; mix : in Vector; fa : in Array_of_Faces;
                 lifted : in Array_of_Lists;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 normals,suppvals : in out List;
                 mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Construction of the mixed subdivision by constraint complementary
  --   slackness: faces supported by the list of normals pass freely.
  --   This list of normals is updated each time new mixed cells are found.
  --   Following are the descriptions of the new parameters.

  -- ON ENTRY :
  --   normals   list of inner normals on the lower hull of the lifted points;
  --   suppvals  corresponding vectors which contain the supporting values
  --             for each component.

  -- ON RETURN :
  --   normals   contains the inner normals of the mixed cells;
  --   suppvals  corresponding vectors of support values.

end Integer_Pruning_Methods;
