with Ada.Text_IO;                       use Ada.Text_IO;
with Double_VecVecs_Interface;
with DCMPLX_VecVecs_Interface;

function use_avvcon ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Double_VecVecs_Interface;
    use DCMPLX_VecVecs_Interface;

  begin
    case job is
      when 0 => return Double_VecVecs_Initialize(a,b,vrblvl);
      when 1 => return DCMPLX_VecVecs_Initialize(a,b,vrblvl);
      when 2 => return Double_VecVecs_Set(a,b,c,vrblvl);
      when 3 => return DCMPLX_VecVecs_Set(a,b,c,vrblvl);
      when 4 => return Double_VecVecs_Get_Dimension(a,vrblvl);
      when 5 => return DCMPLX_VecVecs_Get_Dimension(a,vrblvl);
      when 6 => return Double_VecVecs_Get_Size_Array(a,b,vrblvl);
      when 7 => return DCMPLX_VecVecs_Get_Size_Array(a,b,vrblvl);
      when 8 => return Double_VecVecs_Get_Size_Vector(a,b,vrblvl);
      when 9 => return DCMPLX_VecVecs_Get_Size_Vector(a,b,vrblvl);
      when 10 => return Double_VecVecs_Get(a,b,c,vrblvl);
      when 11 => return DCMPLX_VecVecs_Get(a,b,c,vrblvl);
      when 12 => return Double_VecVecs_Clear(vrblvl);
      when 13 => return DCMPLX_VecVecs_Clear(vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_avvcon;
