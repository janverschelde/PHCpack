with Power_Series_Interface;

function use_series ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function do_jobs return integer32 is

    use Power_Series_Interface;

  begin
    case job is
      when 1 => return Series_Standard_Newton_at_Point(a,b,vrblvl);
      when 2 => return Series_DoblDobl_Newton_at_Point(a,b,vrblvl);
      when 3 => return Series_QuadDobl_Newton_at_Point(a,b,vrblvl);
      when 4 => return Series_Standard_Newton_at_Series(a,b,vrblvl);
      when 5 => return Series_DoblDobl_Newton_at_Series(a,b,vrblvl);
      when 6 => return Series_QuadDobl_Newton_at_Series(a,b,vrblvl);
      when 7 => return Series_Standard_Pade(a,b,vrblvl);
      when 8 => return Series_DoblDobl_Pade(a,b,vrblvl);
      when 9 => return Series_QuadDobl_Pade(a,b,vrblvl);
      when others => return -1;
    end case;
  end do_jobs;

begin
  return do_jobs;
end use_series;
