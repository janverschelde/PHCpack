with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Scaling_Interface;

function use_scaling ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Scaling_Interface;

  begin
    case job is
      when 1 => return Scale_Standard_System(a,c,vrblvl);
      when 2 => return Scale_DoblDobl_System(a,c,vrblvl);
      when 3 => return Scale_QuadDobl_System(a,c,vrblvl);
      when 4 => return Scale_Multprec_System(a,c,vrblvl);
      when 5 => return Scale_Standard_Solutions(a,b,c,vrblvl);
      when 6 => return Scale_DoblDobl_Solutions(a,b,c,vrblvl);
      when 7 => return Scale_QuadDobl_Solutions(a,b,c,vrblvl);
      when 8 => return Scale_Multprec_Solutions(a,b,c,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_scaling handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_scaling;
