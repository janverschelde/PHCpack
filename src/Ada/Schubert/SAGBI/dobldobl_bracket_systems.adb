package body DoblDobl_Bracket_Systems is

  procedure Clear ( s : in out Bracket_System ) is
  begin
    for i in s'range loop
      Clear(s(i));
    end loop;
  end Clear;

end DoblDobl_Bracket_Systems;
