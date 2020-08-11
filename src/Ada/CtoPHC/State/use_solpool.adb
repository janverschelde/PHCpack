with text_io;                           use text_io;
with Standard_SolsPool_Interface;

function use_solpool ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Standard_SolsPool_Interface;

  begin
    case job is
      when 0 => return Standard_SolsPool_Initialize(a,vrblvl);
      when 1 => return Standard_SolsPool_Size(a,vrblvl);
      when 2 => return Standard_SolsPool_Length(a,b,vrblvl);
      when 3 => return Standard_SolsPool_Dimension(a,b,vrblvl);
      when 4 => return Standard_SolsPool_Add(a,b,c,vrblvl);
      when 5 => return Standard_SolsPool_Get(a,b,c,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_solpool;
