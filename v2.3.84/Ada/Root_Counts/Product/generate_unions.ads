with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

generic

  type boolean_array is array ( natural32 range <> ) of boolean;
  with procedure Process ( ar : in boolean_array; continue : out boolean );

procedure Generate_Unions ( k,first,last : in natural32 );

-- DESCRIPTION :
--   This procedure generates all possible unions of k elements in the
--   range first..last.
