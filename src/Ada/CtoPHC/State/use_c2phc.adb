with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with use_outdata; -- output data for DEMiCs
with use_c2phc4c;

function use_c2phc ( job : integer32;
                     a : C_intarrs.Pointer;
		     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32 is

  pragma Unreserve_All_Interrupts;

  function Handle_Jobs return integer32 is
  begin
    case job is
     -- extract DEMiCs output data
      when 834 => return use_outdata(0,a,b,c); -- allocate memory for lifting
      when 835 => return use_outdata(1,a,b,c); -- assign a lifting value
      when 836 => return use_outdata(2,a,b,c); -- retrieve a lifting value
      when 837 => return use_outdata(3,a,b,c); -- clear lifting values
      when 838 => return use_outdata(4,a,b,c); -- append cell indices
      when 839 => return use_outdata(5,a,b,c); -- retrieve cell indices
      when 840 => return use_outdata(6,a,b,c); -- clear cell indices
      when 841 => return use_outdata(7,a,b,c); -- store mixed volume
      when 842 => return use_outdata(8,a,b,c); -- retrieve mixed volume
      when 843 => return use_outdata(9,a,b,c); -- call DEMiCs for mixed volume
      when 844 => return use_outdata(10,a,b,c); -- stable mv by DEMiCs
      when others => return use_c2phc4c(job,a,b,c,vrblvl);
    end case;
  exception
    when others => put("Exception raised in use_c2phc handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 raise; -- return job;
end use_c2phc;
