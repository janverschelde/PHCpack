with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;

package Pack_Continuation_Parameters is

-- DESCRIPTION :
--   This package provides a compact interface to the settings of the
--   values for the continuation parameters.  In particular, it defines
--   an order on the continuation parameters, mapping every constant
--   defined in the package Continuation_Parameters to a unique number
--   in the range 1 to 34.

  function Get return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..34 with all 34 current values
  --   of the continuation parameters.

  function Get_Value ( k : natural32 ) return double_float;

  -- DESCRIPTION :
  --   If k is in range 1..34, then the value of the k-th parameter
  --   is returned, otherwise -1.0 is returned.

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

  procedure Set_Value ( k : in natural32; v : in double_float );

  -- DESCRIPTION :
  --   Given in k a number in the range 1..34 and in v a corresponding
  --   value for the k-th continuation parameter, the k-th value of the
  --   continuation parameter is set to v.

end Pack_Continuation_Parameters;
