with text_io;                           use text_io;
with Monomial_Maps_Interface;

function use_mapcon ( job : integer32;
                      a : C_intarrs.Pointer;
		      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Monomial_Maps_Interface;

  begin
    case job is
      when 0 => return Monomial_Maps_Solve(a,vrblvl);
      when 1 => return Monomial_Maps_Write(vrblvl);
      when 2 => return Monomial_Maps_Clear(vrblvl);
      when 3 => return Monomial_Maps_Top_Dimension(a,vrblvl);
      when 4 => return Monomial_Maps_Size(a,b,vrblvl);
      when 5 => return Monomial_Maps_Degree(a,b,vrblvl);
      when 6 => return Monomial_Maps_Coefficients(a,c,vrblvl);
      when 7 => return Monomial_Maps_Exponents(a,b,vrblvl);
      when 8 => return Monomial_Maps_Data(a,b,c,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_mapcon;
