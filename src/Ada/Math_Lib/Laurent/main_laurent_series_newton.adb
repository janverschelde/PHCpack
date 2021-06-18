with text_io;                            use text_io;

package body Main_Laurent_Series_Newton is

  procedure Run_Laurent_Series_Newton
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is
  begin
    if vrb > 0 then
      put("-> in main_laurent_series_newton.");
      put_line("Run_Laurent_Series_Newton ...");
    end if;
  end Run_Laurent_Series_Newton;

end Main_Laurent_Series_Newton;
