with text_io;                           use text_io;
with Tableau_Form_Interface;

function use_tabform ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Tableau_Form_Interface;
  
  begin
    case job is
      when 0 => return Tableau_Form_Store(a,b,c,vrblvl);
      when 1 => return Tableau_Form_Dimensions(a,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_tabform;
