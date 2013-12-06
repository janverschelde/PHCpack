with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Numbers_io;                        use Numbers_io;
with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Multprec_Parse_Numbers;            use Multprec_Parse_Numbers;
with Symbol_Table_io;
with Multprec_Complex_Laurentials;      use Multprec_Complex_Laurentials;
with Multprec_Complex_Laurentials_io;   use Multprec_Complex_Laurentials_io;
with Parse_Polynomial_Exceptions;       use Parse_Polynomial_Exceptions;

package body Multprec_Complex_Laur_Systems_io is

-- EXCEPTION HANDLER :

  procedure Write_Symbol_Table is

  -- DESCRIPTION :
  --   Writes the current list of symbols on one line on standard output.

  begin
    put("Current symbols : ");
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(Symbol_Table.Get(i)); put(" ");
    end loop;
    new_line;
  end Write_Symbol_Table;

  procedure Handler ( k : in integer32 ) is
  begin
    put(" raised while reading polynomial ");
    put(k,1);
    put_line(".");
  end Handler;

-- INPUT ROUTINES :

  procedure get ( p : in out Laur_Sys ) is
  begin
    get(standard_input,p);
  end get;

  procedure get ( file : in file_type; p : in out Laur_Sys ) is

    i : integer32 := p'first;

  begin
    while i <= p'last loop
      get(file,p(i));
      i := i+1;
    end loop;
  exception
    when ILLEGAL_CHARACTER    => put("ILLEGAL_CHARACTER");
                                 Handler(i); raise;
    when INVALID_SYMBOL       => put("INVALID_SYMBOL");
                                 Handler(i); raise;
    when ILLEGAL_OPERATION    => put("ILLEGAL_OPERATION");
                                 Handler(i); raise;
    when INFINITE_NUMBER      => put("INFINITE_NUMBER");
                                 Handler(i); raise;
    when OVERFLOW_OF_UNKNOWNS => put("OVERFLOW_OF_UNKNOWNS");
                                 Handler(i);
                                 Write_Symbol_Table; raise;
    when BAD_BRACKET          => put("BAD_BRACKET");
                                 Handler(i); raise;
    when others => put("UNKNOWN ERROR"); Handler(i); raise;
  end get;

  procedure get ( n : in out natural32; p : in out Laur_Sys ) is
  begin
    get(standard_input,n,p);
  end get;

  procedure get ( file : in file_type;
                  n : in out natural32; p : in out Laur_Sys ) is

    i : integer32 := p'first;

  begin
    get(file,n);
    while i <= p'first+integer32(n)-1 loop
      get(file,p(i));
      i := i+1;
    end loop;
  exception
    when ILLEGAL_CHARACTER    => put("ILLEGAL_CHARACTER");
                                 Handler(i); raise;
    when INVALID_SYMBOL       => put("INVALID_SYMBOL");
                                 Handler(i); raise;
    when ILLEGAL_OPERATION    => put("ILLEGAL_OPERATION");
                                 Handler(i); raise;
    when INFINITE_NUMBER      => put("INFINITE_NUMBER");
                                 Handler(i); raise;
    when OVERFLOW_OF_UNKNOWNS => put("OVERFLOW_OF_UNKNOWNS");
                                 Handler(i);
                                 Write_Symbol_Table; raise;
    when BAD_BRACKET          => put("BAD_BRACKET");
                                 Handler(i); raise;
    when others => put("UNKNOWN ERROR"); Handler(i); raise;
  end get;

-- MORE USER FRIENDLY INPUT OPERATIONS :

  procedure get ( lp : out Link_to_Laur_Sys ) is

    inpt : file_type;
    n,m : natural32;
    onfile : character;

  begin
    loop
      put("Is the system on a file ? (y/n/i=info) ");
      Ask_Alternative(onfile,"yni");
      if onfile = 'i'
       then new_line; Display_Format; new_line;
      end if;
      exit when onfile /= 'i';
    end loop;
    new_line;
    if onfile = 'y' then
      declare
        procedure Read is
        begin	 
          put_line("Reading the name of the input file.");
          Read_Name_and_Open_File(inpt);
          get(inpt,n);
        end Read;
      begin
        Read;
      exception
        when others => put_line("The data on the file is not correct.");
                       put_line("A natural number is expected first.");
                       put_line("Supply another file name."); Close(inpt);
                       Read;
      end;
    else
      put("Give the number of polynomials : "); Read_Natural(n);
    end if;
    lp := new Laur_Sys(1..integer32(n));
    declare
      procedure Read is
      begin 
        if onfile = 'y' then
          m := natural32(Scan_Line_for_Number(inpt));
          if Symbol_Table.Empty then
            if m > 0
             then Symbol_Table.Init(m);
             else Symbol_Table.Init(n);
            end if;
          end if;
          get(inpt,lp.all);
          Close(inpt);
        else
          put("Give the number of unknowns : "); Read_Natural(m);
          if Symbol_Table.Empty
           then Symbol_Table.Init(m);
          end if;
          put("Give "); put(n,2);
          if n = 1
           then put_line(" polynomial : ");
           else put_line(" polynomials : ");
          end if;
          get(lp.all);
          skip_line;  -- skip end_of_line symbol
        end if;
      exception
        when others => if onfile = 'y' then Close(inpt); end if;
                       put_line("Polynomial system read : "); put(lp.all);
                       Clear(lp); Symbol_Table.Clear; raise;
      end Read;
    begin
      Read;
    exception
      when others =>
             if onfile = 'y'
              then put_line("The polynomials on the file are incorrect."
                          & " Try again...");
              else put_line("The polynomials are incorrect. Try again...");
             end if;
             Clear(lp); Symbol_Table.Clear;
             get(lp);
    end;
  end get;

  procedure get ( file : in file_type; lp : out Link_to_Laur_Sys ) is

    n,m : natural32 := 0;

  begin
    get(file,n);
    declare
      p : Laur_Sys(1..integer32(n));
    begin
      m := natural32(Scan_Line_for_Number(file));
      if Symbol_Table.Empty then
        if m > 0
         then Symbol_Table.Init(m);
         else Symbol_Table.Init(n);
        end if;
      end if;
      get(file,p);
      lp := new Laur_Sys'(p);
    end;
   -- skip_line(file);    -- skip end_of_line symbol
   -- this skip causes exception "END_ERROR" to be raised
   -- if the system is not followed by anything else!
  end get;

-- OUTPUT ROUTINES :

  procedure put ( p : in Laur_Sys ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( n : in natural32; p : in Laur_Sys ) is
  begin
    put(standard_output,n,p);
  end put;

  procedure put ( file : in file_type; p : in Laur_Sys ) is
  begin
    for i in p'range loop
      put(file,p(i));
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type; n : in natural32; p : in Laur_Sys ) is
  begin
    put(file,n,1);
    new_line(file);
    for i in p'range loop
      put(file,p(i));
      new_line(file);
    end loop;
  end put;

  procedure put ( p : in Laur_Sys; s : in Array_of_Symbols ) is
  begin
    put(standard_output,p,s);
  end put;

  procedure put ( n : in natural32;
                  p : in Laur_Sys; s : in Array_of_Symbols ) is
  begin
    put(standard_output,n,p,s);
  end put;

  procedure put ( file : in file_type;
                  p : in Laur_Sys; s : in Array_of_Symbols ) is
  begin
    for i in p'range loop
      put(file,p(i),s);
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type; n : in natural32;
                  p : in Laur_Sys; s : in Array_of_Symbols ) is
  begin
    put(file,n,1);
    new_line(file);
    put(file,p,s);
  end put;

  procedure put_line ( p : in Laur_Sys ) is
  begin
    put_line(standard_output,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Laur_Sys ) is

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    nq : constant natural32 := natural32(p'length);

  begin
    put(file,nq,2);
    if n /= nq
     then put(file," "); put(file,n,1);
    end if;
    new_line(file);
    for i in p'range loop
      put_line(file,p(i));
      new_line(file);
    end loop;
  end put_line;

  procedure put_line ( p : in Laur_Sys; s : in Array_of_Symbols ) is
  begin
    put_line(standard_output,p,s);
  end put_line;

  procedure put_line ( file : in file_type;
                       p : in Laur_Sys; s : in Array_of_Symbols ) is

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    nq : constant natural32 := natural32(p'length);

  begin
    put(file,nq,2);
    if n /= nq
     then put(file," "); put(file,n,1);
    end if;
    new_line(file);
    for i in p'range loop
      put_line(file,p(i),s);
      new_line(file);
    end loop;
  end put_line;

end Multprec_Complex_Laur_Systems_io;
