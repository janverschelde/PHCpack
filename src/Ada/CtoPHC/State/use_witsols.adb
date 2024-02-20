with Irreducible_Components_Interface;

function use_witsols ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function do_jobs return integer32 is

    use Irreducible_Components_Interface;

  begin
    case job is
      when  0 => return Standard_Polynomial_Solver(a,b,vrblvl);
      when  1 => return Standard_Laurent_Solver(a,b,vrblvl);
      when  2 => return DoblDobl_Polynomial_Solver(a,b,vrblvl);
      when  3 => return DoblDobl_Laurent_Solver(a,b,vrblvl);
      when  4 => return QuadDobl_Polynomial_Solver(a,b,vrblvl);
      when  5 => return QuadDobl_Laurent_Solver(a,b,vrblvl);
      when  6 => return Standard_Polynomial_WitSet_Copy(a,vrblvl);
      when  7 => return Standard_Laurent_WitSet_Copy(a,vrblvl);
      when  8 => return DoblDobl_Polynomial_WitSet_Copy(a,vrblvl);
      when  9 => return DoblDobl_Laurent_WitSet_Copy(a,vrblvl);
      when 10 => return QuadDobl_Polynomial_WitSet_Copy(a,vrblvl);
      when 11 => return QuadDobl_Laurent_WitSet_Copy(a,vrblvl);
      when 12 => return Standard_WitSet_Clear(vrblvl);
      when 13 => return DoblDobl_WitSet_Clear(vrblvl);
      when 14 => return QuadDobl_WitSet_Clear(vrblvl);
      when 15 => return Irreducible_Factor_String(a,b,vrblvl);
      when others => return -1;
    end case;
  end do_jobs;

begin
  return do_jobs;
end use_witsols;
