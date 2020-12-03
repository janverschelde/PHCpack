with Ada.Calendar;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Characters_and_Numbers;
with Standard_Parse_Numbers;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_Strings;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions; 

package body Parse_Strings_to_Polynomials is

  procedure First_Line_Format ( file : in file_type ) is

  -- DESCRIPTION :
  --   Writes the format for the first line to file.

  begin
    put_line(file,
             "The first line on input should start with a positive integer,");
    put_line(file,
             "equal to the number of polynomials in the system.");
    put_line(file,
             "The optional second positive number on the first line equals");
    put_line(file,
             "the number of variables in the system.");
    put_line(file,
             "If the second number on the first input line is omitted,");
    put_line(file,
             "then we assume the system has as many variables as equations.");
    put_line(file,
             "If the system is not square (over or underdetermined),");
    put_line(file,
             "then the second number on the first line is not optional.");
  end First_Line_Format;

  procedure Parse_Dimensions
              ( outfile : in file_type; to_file : in boolean;
                s : in string; nq,nv : out natural32;
                fail : out boolean ) is

    p : integer := s'first;
    i : integer32 := 0;
    ni : natural32 := 0;
    sign : character := '+';

  begin
    fail := false;
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if p > s'last then
      if to_file then
        new_line(outfile);
        put_line(outfile,"The first line contains only spaces.");
      else
        new_line;
        put_line("The first line contains only spaces.");
      end if;
      fail := true;
    elsif Characters_and_Numbers.Convert(s(p)) >= 10 then
      if to_file then
        new_line(outfile);
        put_line(outfile,
                 "The first nonspace on the first line is not a number.");
        put_line(outfile,"The character is " & s(p) & ".");
      else
        new_line;
        put_line("The first nonspace on the first line is not a number.");
        put_line("The character is " & s(p) & ".");
      end if;
      fail := true;
    else
      Standard_Parse_Numbers.Parse(s,p,i,ni,sign);
      if not to_file then
        new_line;
        put_line("Parsing the first line of the input...");
        put("The first number parsed : "); put(i,1); put_line(".");
      end if;
      nq := natural32(i);
      if p > s'last then
        if not to_file then
          put_line("At the end of the line, assuming the system is square.");
        end if;
        nv := nq;
      else
        Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
        if p > s'last then
          if not to_file then
            put_line("No second number found, assuming the system is square.");
          end if;
          nv := nq;
        elsif Characters_and_Numbers.Convert(s(p)) >= 10 then
          if to_file then
            new_line(outfile);
            put_line(outfile,
                     "The second nonspace on the first line is not a number.");
          else
            new_line;
            put_line("The second nonspace on the first line is not a number.");
          end if;
          fail := true;
        else
          Standard_Parse_Numbers.Parse(s,p,i,ni,sign);
          if not to_file
           then put("The second number parsed : "); put(i,1); put_line(".");
          end if;
          nv := natural32(i);
        end if;
      end if;
    end if;
  end Parse_Dimensions;

  procedure Read_First_Line
              ( infile,outfile : in file_type; to_file : in boolean;
                nq,nv : out natural32; fail : out boolean ) is

    tmp : string(1..256);
    cnt : integer := 0;

  begin
    get_line(infile,tmp,cnt);
   -- put("Read "); put(cnt,1); put_line(" characters on the first line.");
   -- put_line("The first line follows :"); 
   -- put_line(tmp(1..cnt));
    Parse_Dimensions(outfile,to_file,tmp(1..cnt),nq,nv,fail);
    if fail then
      if to_file
       then new_line(outfile); First_Line_Format(outfile);
       else new_line; First_Line_Format(standard_output);
      end if;
    end if;
  end Read_First_Line;

  procedure Read_Polynomials
             ( infile,outfile : in file_type; to_file : in boolean;
               n : in natural32;
               la : out String_Splitters.Link_to_Array_of_Strings;
               fail : out boolean ) is

    use String_Splitters;
    a : Array_of_Strings(1..integer(n));

  begin
    a := Read_till_Delimiter(infile,integer(n),';');
    if not to_file then
      new_line;
      put_line("The strings read from file : ");
    end if;
    fail := false;
    for i in a'range loop
      if a(i) = null then
        if to_file then
          put(outfile,"There is no polynomial ");
          put(outfile,integer32(i),1);
          put_line(outfile,", missing ';' symbol at end of polynomial?");
        else
          put("There is no polynomial "); put(integer32(i),1);
          put_line(", missing ';' symbol at end of polynomial?");
        end if;
        fail := true;
      end if;
      exit when fail;
      if not to_file
       then put_line(a(i).all);
      end if;
    end loop;
    if not fail
     then la := new Array_of_Strings'(a);
    end if;
  exception
    when others =>
      if to_file then
        put(outfile,"There are no "); put(outfile,n,1);
        put_line(outfile," ';' on file.");
      else
        put("There are no "); put(n,1); put_line(" ';' on file.");
      end if;
      fail := true;
  end Read_Polynomials;

  procedure Parse_Polynomials
              ( file : in file_type; to_file : in boolean;
                n : in natural32;
                a : in String_Splitters.Link_to_Array_of_Strings;
                p : in Link_to_Laur_Sys; fail : out integer32 ) is

    i : integer32 := integer32(a'first);

    procedure Write_Error ( index : in integer32 ) is

    -- DESCRIPTION :
    --   Writes a generic error message when the error occurs
    --   in the string with the given index.

    begin
      if to_file then
        put(file,"Error raised in string "); put(file,index,1);
        put_line(file,", the wrong input is on the next line");
        put_line(file,a(integer(index)).all);
      else
        put("Error raised in string "); put(index,1);
        put_line(", the wrong input is on next line");
        put_line(a(integer(index)).all);
      end if;
    end Write_Error;

  begin
    fail := 0; 
    while i <= integer32(a'last) loop
      exit when (fail > 0);
      declare
        q : Poly;
        s : constant string := a(integer(i)).all;
      begin
        q := Standard_Complex_Laur_Strings.Parse(n,s);
        p(i) := q;
        i := i + 1;
      exception
        when ILLEGAL_CHARACTER =>
          fail := i; Write_Error(i);
          if to_file then
            put(file,"An invalid character was found at polynomial ");
            put(file,i,1); put_line(file,".");
          else
            put("An invalid character was found at polynomial ");
            put(i,1); put_line(".");
          end if;
          raise;
        when ILLEGAL_OPERATION =>
          fail := i; Write_Error(i);
          if to_file then
            put(file,"An invalid operation was found at polynomial ");
            put(file,i,1); put_line(file,".");
          else
            put("An invalid operation was found at polynomial ");
            put(i,1); put_line(".");
          end if;
          raise;
        when OVERFLOW_OF_UNKNOWNS =>
          fail := i; Write_Error(i);
          if to_file then
            put(file,"Too many unknowns at polynomial ");
            put(file,i,1); put_line(file,".");
            put(file,"The current symbol table : ");
            Symbol_Table_io.Write(file); new_line(file);
          else
           -- put("Symbol_Table.Number : "); put(Symbol_Table.Number,1);
            new_line;
            put("Too many unknowns at polynomial ");
            put(fail,1); put_line(".");
           -- put("The current symbol table : ");
           -- Symbol_Table_io.Write; new_line;
          end if;
          raise;
        when BAD_BRACKET =>
          fail := i; Write_Error(i);
          if to_file then
            put(file,"Misplaced bracket at polynomial ");
            put(file,i,1); put_line(file,".");
          else
            put("Misplaced bracket at polynomial ");
            put(i,1); put_line(".");
          end if;
          raise;
        when others =>
          if to_file then
            put(file,"Some error occurred at polynomial ");
            put(file,i,1); put_line(file,".");
          else
            put("Some error occurred at polynomial ");
            put(i,1); put_line(".");
          end if;
          fail := i; Write_Error(i);
          raise;
      end;
    end loop;
   --  exception
   -- when OVERFLOW_OF_UNKNOWNS
   --   => put("fail = "); put(fail,1); new_line;
   --      put_line("trying to stay silent ..."); return;
   -- a simple return does not work with -O3 inlining option ...
   -- when others => raise;
  end Parse_Polynomials;

  procedure Create_Output_File
              ( file : in out file_type; name : in string;
                no_output_file : out boolean ) is
  begin
    if name = "" then
      declare
        ans : character;
      begin
        new_line;
        put("Do you want the output to file ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'n' then
          no_output_file := true;
        else
          no_output_file := false;
          Read_Name_and_Create_File(file);
        end if;
      exception
        when others =>
          new_line;
          put_line("Could not create the output file.");
      end;
    else
      no_output_file := false;
     -- new_line;
     -- put_line("Creating the file with name '" & name & "'");
      declare
      begin
        Create_Output_File(file,name);
      exception
        when others =>
          new_line;
          put("Could not create the output file with name '");
          put_line(name & "'.");
      end;
    end if;
  end Create_Output_File;

  procedure Write_Results
              ( file : in file_type;
                nq,nv : in natural32;
                a : in String_Splitters.Link_to_Array_of_Strings;
                p : in Link_to_Laur_Sys; fail_index : in integer32 ) is

    now : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  begin
    if fail_index = 0 then
      put(file,nq,1);
      if nv /= nq 
       then put(file," "); put(file,nv,1);
      end if;
      new_line(file);
      for i in a'range loop
        put_line(file,a(i).all);
      end loop;
      new_line(file);
      put(file,"PHC parsed input successfully on ");
      Time_Stamps.Write_Time_Stamp(file,now); put_line(file,".");
      new_line(file);
      put_line(file,"The symbols of the variables :");
      Symbol_Table_io.Write(file);
      new_line(file);
      put_line(file,"The polynomials parsed with standard doubles :");
      put(file,p.all); new_line(file);
    else
      new_line(file);
      put(file,"PHC failed to parse input on ");
      Time_Stamps.Write_Time_Stamp(file,now); put_line(file,".");
    end if;
  end Write_Results;

  procedure Read_from_File ( inname,outname : in string ) is

    now : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    infile,outfile : file_type;
    to_file : boolean := false;
    nq,nv : natural32;
    fail,nofile : boolean;
    las : String_Splitters.Link_to_Array_of_Strings;
    fail_index : integer32 := 0;
    lq : Link_to_Laur_Sys;

  begin
    Open(infile,in_file,inname);
    Create_Output_File(outfile,outname,nofile);
    to_file := not nofile;
    Read_First_Line(infile,outfile,to_file,nq,nv,fail);
    if fail then
      if to_file then
        new_line(outfile);
        put(outfile,"PHC failed to parse input on file '" & inname);
        put_line(outfile,"',");
        put(outfile,"on ");
        Time_Stamps.Write_Time_Stamp(outfile,now); put_line(outfile,".");
      else
        new_line;
        put("PHC failed to parse input on file '" & inname);
        put_line("',");
        put("on ");
        Time_Stamps.Write_Time_Stamp(standard_output,now); put_line(".");
      end if;
    else
      Read_Polynomials(infile,outfile,to_file,nq,las,fail);
      if not fail then
        lq := new Laur_Sys(1..integer32(nq));
        if not to_file then
          new_line;
          put("Parsing "); put(nq,1); put(" strings in ");
          put(nv,1); put_line(" variables...");
        end if;
       -- if to_file
       --  then put_line("Writing error messages to file.");
       --  else put_line("Writing error messages to screen.");
       -- end if;
        Symbol_Table.Init(nv);
        Parse_Polynomials(outfile,to_file,nv,las,lq,fail_index);
        if to_file
         then Write_Results(outfile,nq,nv,las,lq,fail_index);
         else Write_Results(standard_output,nq,nv,las,lq,fail_index);
        end if;
      end if;
    end if;
  exception
    when others =>
      if to_file then
        new_line(outfile);
        put(outfile,"Parse error in the file with name '");
        put(outfile,inname); put_line(outfile,"'."); raise;
      else
        new_line;
        put("Parse error in the file with name '"); put(inname);
        put_line("'."); raise;
      end if;
  end Read_from_File;

  procedure Read_Input_File_Name is

    name : constant string := String_Splitters.Read_String;

  begin
    Read_from_File(name,"");
 -- exception
 --   when others => raise; -- put_line("suppressing exception..."); return;
  end Read_Input_File_Name;

  procedure Main ( infilename,outfilename : in string ) is
  begin
    if infilename /= "" then
      Read_from_File(infilename,outfilename);
    else
      new_line;
      put_line("Reading the name of the input file ...");
      Read_Input_File_Name; 
    end if;
  end Main;

end Parse_Strings_to_Polynomials;
