with text_io;                            use text_io;
with Standard_Floating_Vectors;

procedure Pruning_Statistics
              ( file : in file_type;
                nbsucc,nbfail : in Standard_Floating_Vectors.Vector );

-- DESCRIPTION :
--   Writes statistics on the number of face-face combinations.
--   These statistics give the user an idea of the pruning tree.

-- ON ENTRY :
--   file       file must be opened for output;
--   nbsucc     number of successul pruning combinations per level;
--   nbfail     number of failing pruning combinations per level.
