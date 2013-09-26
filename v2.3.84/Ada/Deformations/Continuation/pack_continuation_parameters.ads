with text_io;                            use text_io;
with Standard_Floating_Vectors;

package Pack_Continuation_Parameters is

-- DESCRIPTION :
--   This package provides a compact interface to the settings of the
--   values for the continuation parameters.

  function Get return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..34 with all 34 current values
  --   of the continuation parameters.

  procedure Write ( v : in Standard_Floating_Vectors.Vector );
  procedure Write ( file : in file_type;
                    v : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the vector v to screen or to file, according to the
  --   original format of the values.

  procedure Set ( v : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Given a vector of range 1..34 with values for the continuation
  --   parameters, the values in Continuation_Parameters are changed.

end Pack_Continuation_Parameters;
