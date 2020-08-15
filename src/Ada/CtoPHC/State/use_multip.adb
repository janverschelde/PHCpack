with text_io;                           use text_io;
with Deflation_Interface;

function use_multip ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Deflation_Interface;

  begin
    case job is
      when 0 => return Deflation_Standard_Multiplicity(a,b,c,vrblvl-1);
      when 1 => return Deflation_DoblDobl_Multiplicity(a,b,c,vrblvl-1);
      when 2 => return Deflation_QuadDobl_Multiplicity(a,b,c,vrblvl-1);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_multip;
