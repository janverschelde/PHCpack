with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Univariate_Solvers_Interface;

function unisolve ( job : integer32;
                    a : C_intarrs.Pointer;
                    b : C_intarrs.Pointer;
                    c : C_dblarrs.Pointer;
                    vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Univariate_Solvers_Interface;

  begin
    case job is
      when 1 => return Standard_Univariate_Solver(a,b,c,vrblvl);
      when 2 => return DoblDobl_Univariate_Solver(a,b,c,vrblvl);
      when 3 => return QuadDobl_Univariate_Solver(a,b,c,vrblvl);
      when 4 => return Multprec_Univariate_Solver(a,b,c,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in unisolve handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end unisolve;
