with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Ada.Characters.Latin_1;

procedure mainfeed ( infilename,outfilename : in string;
                     verbose : in integer32 := 0 ) is

  function main_feedback ( input,output : string ) return integer;
  pragma Import(C,main_feedback,"main_feedback");

  procedure Main is

  -- DESCRIPTION :
  --   Checks whether the infilename and outfilename are not empty,
  --   appends the '\0' end of string character to the names,
  --   before calling the C function main_feedback.

    NUL : constant character := Ada.Characters.Latin_1.NUL;
    res : integer;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in mainfeed.Main ...");
    end if;
    if infilename = "" or else outfilename = "" then
      new_line;
      put_line("Usage: phc -k input_file output_file");
      new_line;
    else
      if verbose > 0 then
        put_line("reading from " & infilename);
        put_line("writing to " & outfilename);
        put_line("Calling feedback...");
      end if;
      res := main_feedback(infilename & NUL,outfilename & NUL);
    end if;
  end Main;

begin
  Main;
end mainfeed;
