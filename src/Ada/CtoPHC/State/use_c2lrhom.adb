with text_io;                           use text_io;
with Schubert_Interface;

function use_c2lrhom ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Schubert_Interface;

  begin
    case job is
      when 0 => return Schubert_Intersection_Conditions(a,b,c,vrblvl);
      when 1 => return Standard_LR_Homotopies(a,b,c,vrblvl);
      when 2 => return DoblDobl_LR_Homotopies(a,b,c,vrblvl);
      when 3 => return QuadDobl_LR_Homotopies(a,b,c,vrblvl);
      when others => put_line("  Sorry.  Invalid operation in use_c2lrhom.");
                     return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2lrhom;
