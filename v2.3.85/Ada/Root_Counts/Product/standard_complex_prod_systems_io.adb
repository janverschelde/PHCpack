with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Symbol_Table;
with Standard_Complex_Poly_Lists_io;     use Standard_Complex_Poly_Lists_io;

package body Standard_Complex_Prod_Systems_io is

  procedure get ( lp : out Link_to_Prod_Sys ) is

    ans : character;
    file : file_type;
    n,m : natural32 := 0;

  begin
    put("Is the product system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the input file.");
      Read_Name_and_Open_File(file);
      get(file,lp);
    else
      put("Give the number of equations : "); get(n);
      put("Give the number of unknowns : "); get(m);
      if Symbol_Table.Empty
       then Symbol_Table.Init(m);
      end if;
      put("Give "); put(n,1); put(" product polynomials in ");
      put(m,1); put_line(" unknowns, terminate by semicolon :");
      declare
        p : Prod_Sys(1..integer32(n));
      begin
        get(p);
        lp := new Prod_Sys'(p);
      end;
    end if;
  end get;

  procedure get ( file : in file_type; lp : out Link_to_Prod_Sys ) is

    n,m : natural32 := 0;

  begin
    get(file,n);
    m := natural32(Scan_Line_for_Number(file));
    if Symbol_Table.Empty then
      if m > 0 
       then Symbol_Table.Init(m);
       else Symbol_Table.Init(n);
      end if;
    end if;
    declare
      p : Prod_Sys(1..integer32(n));
    begin
      get(file,p);
      lp := new Prod_Sys'(p);
    end;
  end get;

  procedure get ( p : out Prod_Sys ) is
  begin
    get(Standard_Input,p);
  end get;

  procedure get ( file : in file_type; p : out Prod_Sys ) is
  begin
    for i in p'range loop
      get(file,p(i));
    end loop;
  end get;

  procedure put ( p : in Prod_Sys ) is
  begin
    put(Standard_Output,p);
  end put;

  procedure put ( n : in natural32; p : in Prod_Sys ) is
  begin
    put(Standard_Output,n,p);
  end put;

  procedure put ( n,m : in natural32; p : in Prod_Sys ) is
  begin
    put(Standard_Output,n,m,p);
  end put;

  procedure put ( file : in file_type; p : in Prod_Sys ) is
  begin
    for i in p'range loop
      put(file,p(i)); new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  n : in natural32; p : in Prod_Sys ) is
  begin
    put(file,n,1); new_line(file);
    for i in p'range loop
      put(file,p(i)); new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  n,m : in natural32; p : in Prod_Sys ) is
  begin
    put(file,n,1); put(file,"  "); put(file,m,1); new_line(file);
    for i in p'range loop
      put(file,p(i)); new_line(file);
    end loop;
  end put;

  procedure put_line ( p : in Prod_Sys ) is
  begin
    put_line(Standard_Output,p);
  end put_line;

  procedure put_line ( n : in natural32; p : in Prod_Sys ) is
  begin
    put_line(Standard_Output,n,p);
  end put_line;

  procedure put_line ( n,m : in natural32; p : in Prod_Sys ) is
  begin
    put_line(Standard_Output,n,m,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Prod_Sys ) is
  begin
    for i in p'range loop
      put_line(file,p(i)); new_line(file);
    end loop;
  end put_line;

  procedure put_line ( file : in file_type;
                       n : in natural32; p : in Prod_Sys ) is
  begin
    put(file,n,1); new_line(file);
    for i in p'range loop
      put_line(file,p(i)); new_line(file);
    end loop;
  end put_line;

  procedure put_line ( file : in file_type;
                       n,m : in natural32; p : in Prod_Sys ) is
  begin
    put(file,n,1); put(file,"  "); put(file,m,1); new_line(file);
    for i in p'range loop
      put_line(file,p(i)); new_line(file);
    end loop;
  end put_line;

end Standard_Complex_Prod_Systems_io;
