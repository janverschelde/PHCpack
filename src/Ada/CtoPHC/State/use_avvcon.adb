with Ada.Text_IO;                       use Ada.Text_IO;
with Double_VecVecs_Interface;

function use_avvcon ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Double_VecVecs_Interface;

  begin
    case job is
      when 0 => return Double_VecVecs_Initialize(a,b,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_avvcon;
