with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Giftwrap_Interface;

function use_giftwrap ( job : integer32;
                        a : C_intarrs.Pointer;
                        b : C_intarrs.Pointer;
                        c : C_dblarrs.Pointer;
                        vrblvl : integer32 := 0 ) return integer32 is

  function Handle_Jobs return integer32 is

    use Giftwrap_Interface;

  begin
    case job is
      when  1 => return Giftwrap_Planar_Hull(a,b,vrblvl);
      when  2 => return Giftwrap_Spatial_Hull(a,b,vrblvl);
      when  3 => return Giftwrap_Number_of_Facets(a,b,vrblvl);
      when  4 => return Giftwrap_String_of_Facet(a,b,vrblvl);
      when  5 => return Giftwrap_3d_Clear(vrblvl);
      when  6 => return Giftwrap_4d_Clear(vrblvl);
      when  7 => return Giftwrap_String_Size(a,vrblvl);
      when  8 => return Giftwrap_String_of_Support(b,vrblvl);
      when  9 => return Giftwrap_String_Clear(vrblvl);
      when 10 => return Giftwrap_Laurent_Initial_Form(a,b,vrblvl);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_giftwrap handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_giftwrap;
