with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Linear_Products_Interface;

function use_roco ( job : integer32;
                    a : C_intarrs.Pointer;
                    b : C_intarrs.Pointer;
                    c : C_dblarrs.Pointer;
                    vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Linear_Products_Interface;

  begin
    case job is
      when 0 => return Linear_Products_Structure_Make(vrblvl);
      when 1 => return Linear_Products_Structure_Write(vrblvl);
      when 2 => return Linear_Products_Structure_Bound(a,vrblvl);
      when 3 => return Linear_Products_System_Make(vrblvl);
      when 4 => return Linear_Products_System_Solve(vrblvl);
      when 5 => return Linear_Products_Clear(vrblvl);
      when 6 => return Linear_Products_Structure_String_Get(a,b,vrblvl);
      when 7 => return Linear_Products_Structure_String_Set(a,b,vrblvl);
      when 8 => return Linear_Products_Structure_Check(a,vrblvl);
      when 10 => return Linear_Products_Partition_Make(a,b,vrblvl);
      when 11 => return Linear_Products_Partition_Bound(a,b,vrblvl);
      when 12 => return Linear_Products_Partition_System(a,b,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  exception
    when others => put("Exception raised in use_roco handling job ");
                   put(job,1); put_line(".  Will try to ignore...");
                   return 1;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_roco;
