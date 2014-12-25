with String_Splitters;                   use String_Splitters;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;

package Solution_String_Splitters is

-- DESCRIPTION :
--   The package provides operations to split off the coordinates
--   of the string representation of a solution in PHCpack format.

  function Trim_End_to_Newline ( s : string ) return string;

  -- DESCRIPTION :
  --   The string on return is trimmed from the end till its newline.

  procedure Split_Coordinates
              ( s : in string; m : out integer32;
                t : out Standard_Complex_Numbers.Complex_Number;
                c : out Link_to_String; fail : out boolean );

  -- DESCRIPTION :
  --   Given in s the full string representation of a solution,
  --   returns in m its multiplicity flag, in t the current value
  --   of the continuation parameter, and in c the coordinates of
  --   the solution as a string, provided fail is not true.

  function Coordinates ( s : string ) return Link_to_String;

  -- DESCRIPTION :
  --   Given in s the full string representation of a solution,
  --   returns the coordinates of the solution as a string.

end Solution_String_Splitters;
