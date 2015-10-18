package body QuadDobl_Bracket_Systems is

  procedure Clear ( s : in out Bracket_System ) is
  begin
    for i in s'range loop
      Clear(s(i));
    end loop;
  end Clear;

end QuadDobl_Bracket_Systems;
