with text_io;                            use text_io;
with Reduction_Interface;

function use_reduction ( job : integer32;
                         a : C_intarrs.Pointer;
                         b : C_intarrs.Pointer;
                         c : C_dblarrs.Pointer;
                         vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Reduction_Interface;

  begin
    case job is
      when 1 => return Reduction_Standard_Linear(a,vrblvl);
      when 2 => return Reduction_DoblDobl_Linear(a,vrblvl);
      when 3 => return Reduction_QuadDobl_Linear(a,vrblvl);
      when 4 => return Reduction_Standard_Nonlinear(a,b,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_reduction.");
                   put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_reduction;
