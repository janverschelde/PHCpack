with text_io;                           use text_io;
with Member_Interface;

function use_c2mbt ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Member_Interface;

  begin
    case job is
      when 0 => return Member_Standard_Polynomial_Numbers(a,b,c,vrblvl);
      when 1 => return Member_DoblDobl_Polynomial_Numbers(a,b,c,vrblvl);
      when 2 => return Member_QuadDobl_Polynomial_Numbers(a,b,c,vrblvl);
      when 3 => return Member_Standard_Laurent_Numbers(a,b,c,vrblvl);
      when 4 => return Member_DoblDobl_Laurent_Numbers(a,b,c,vrblvl);
      when 5 => return Member_QuadDobl_Laurent_Numbers(a,b,c,vrblvl);
      when 6 => return Member_Standard_Polynomial_String(a,b,c,vrblvl);
      when 7 => return Member_DoblDobl_Polynomial_String(a,b,c,vrblvl);
      when 8 => return Member_QuadDobl_Polynomial_String(a,b,c,vrblvl);
      when 9 => return Member_Standard_Laurent_String(a,b,c,vrblvl);
      when 10 => return Member_DoblDobl_Laurent_String(a,b,c,vrblvl);
      when 11 => return Member_QuadDobl_Laurent_String(a,b,c,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2mbt;
