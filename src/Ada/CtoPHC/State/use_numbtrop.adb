with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Numerical_Tropisms_Interface;

function use_numbtrop ( job : integer32;
                        a : C_intarrs.Pointer;
                        b : C_intarrs.Pointer;
                        c : C_dblarrs.Pointer;
                        vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Numerical_Tropisms_Interface;

  begin
    case job is
      when 1 => return Standard_Initialize(a,b,c,vrblvl);
      when 2 => return DoblDobl_Initialize(a,b,c,vrblvl);
      when 3 => return QuadDobl_Initialize(a,b,c,vrblvl);
      when 4 => return Store_Standard_Tropism(a,b,c,vrblvl);
      when 5 => return Store_Dobldobl_Tropism(a,b,c,vrblvl);
      when 6 => return Store_Quaddobl_Tropism(a,b,c,vrblvl);
      when 7 => return Standard_Retrieve_All_Tropisms(a,b,c,vrblvl);
      when 8 => return DoblDobl_Retrieve_All_Tropisms(a,b,c,vrblvl);
      when 9 => return QuadDobl_Retrieve_All_Tropisms(a,b,c,vrblvl);
      when 10 => return Standard_Size(a,vrblvl);
      when 11 => return DoblDobl_Size(a,vrblvl);
      when 12 => return QuadDobl_Size(a,vrblvl);
      when 13 => return Standard_Retrieve_One_Tropism(a,b,c,vrblvl);
      when 14 => return DoblDobl_Retrieve_One_Tropism(a,b,c,vrblvl);
      when 15 => return QuadDobl_Retrieve_One_Tropism(a,b,c,vrblvl);
      when 16 => return Standard_Clear(vrblvl);
      when 17 => return DoblDobl_Clear(vrblvl);
      when 18 => return QuadDobl_Clear(vrblvl);
      when 19 => return Standard_Dimension(a,vrblvl);
      when 20 => return DoblDobl_Dimension(a,vrblvl);
      when 21 => return QuadDobl_Dimension(a,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  exception
    when others => put("Exception raised in use_numbtrop handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_numbtrop;
