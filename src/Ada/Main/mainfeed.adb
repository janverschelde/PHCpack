with text_io;                           use text_io;

procedure mainfeed ( infilename,outfilename : in string ) is

  function main_feedback ( input,output : string ) return integer;
  pragma Import(C,main_feedback,"main_feedback");

  procedure Main is

    res : integer := 0;

  begin
    if infilename = "" or else outfilename = "" then
      new_line;
      put_line("Usage: phc -k input_file output_file");
      new_line;
    else
      put_line("reading from " & infilename);
      put_line("writing to " & outfilename);
      put_line("Calling feedback...");
      res := main_feedback(infilename,outfilename);
    end if;
  end Main;

begin
  Main;
end mainfeed;
