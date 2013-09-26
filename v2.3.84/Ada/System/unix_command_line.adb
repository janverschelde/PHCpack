with Ada.Command_Line;

package body Unix_Command_Line is

-- IMPLEMENTATION : using Ada.Command_Line

  function Number_of_Arguments return natural is
  begin
    return Ada.Command_Line.Argument_Count;
  end Number_of_Arguments;

  function Argument ( i : natural ) return string is
  begin
    return Ada.Command_Line.Argument(i);
  end Argument;

end Unix_Command_Line;
