with text_io;                           use text_io;
with Standard_SysPool_Interface;
with DoblDobl_SysPool_Interface;
with QuadDobl_SysPool_Interface;
with Newton_Interface;

function use_syspool ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Standard_SysPool_Interface;
    use DoblDobl_SysPool_Interface;
    use QuadDobl_SysPool_Interface;
    use Newton_Interface;

  begin
    case job is
      when  0 => return Standard_SysPool_Initialize(a,vrblvl-1);
      when  1 => return Standard_SysPool_Size(a,vrblvl-1);
      when  2 => return Standard_SysPool_Read(a,vrblvl-1);
      when  3 => return Standard_SysPool_Write(a,vrblvl-1);
      when  4 => return Standard_SysPool_from_Container(a,vrblvl-1);
      when  5 => return Newton_Standard_SysPool_Refine(a,vrblvl-1);
      when  6 => return Standard_SysPool_into_Container(a,vrblvl-1);
      when  7 => return DoblDobl_SysPool_into_Container(a,vrblvl-1);
      when  8 => return QuadDobl_SysPool_into_Container(a,vrblvl-1);
      when  9 => return DoblDobl_SysPool_Size(a,vrblvl-1);
      when 10 => return QuadDobl_SysPool_Size(a,vrblvl-1);
      when 11 => return DoblDobl_SysPool_Initialize(a,vrblvl-1);
      when 12 => return QuadDobl_SysPool_Initialize(a,vrblvl-1);
      when 13 => return Standard_SysPool_Clear(vrblvl-1);
      when 14 => return DoblDobl_SysPool_Clear(vrblvl-1);
      when 15 => return QuadDobl_SysPool_Clear(vrblvl-1);
      when 16 => return DoblDobl_SysPool_from_Container(a,vrblvl-1);
      when 17 => return QuadDobl_SysPool_from_Container(a,vrblvl-1);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_syspool;
