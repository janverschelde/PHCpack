with text_io;                           use text_io;
with Ada.Characters.Latin_1;

procedure mainfeed ( infilename,outfilename : in string ) is

  function main_feedback ( input,output : string ) return integer;
  pragma Import(C,main_feedback,"main_feedback");

  procedure Main is

  -- DESCRIPTION :
  --   Checks whether the infilename and outfilename are not empty,
  --   appends the '\0' end of string character to the names,
  --   before calling the C function main_feedback.

    res : integer := 0;
    NUL : constant character := Ada.Characters.Latin_1.NUL;

  begin
    if infilename = "" or else outfilename = "" then
      new_line;
      put_line("Usage: phc -k input_file output_file");
      new_line;
    else
     -- put_line("reading from " & infilename);
     -- put_line("writing to " & outfilename);
     -- put_line("Calling feedback...");
      res := main_feedback(infilename & NUL,outfilename & NUL);
    end if;
  end Main;

begin
  Main;
end mainfeed;
