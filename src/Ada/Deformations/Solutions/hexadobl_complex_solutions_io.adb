with File_Scanning;                      use File_Scanning;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Solutions_io;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with HexaDobl_Complex_Numbers_io;        use HexaDobl_Complex_Numbers_io;
with Symbol_Table;                       use Symbol_Table;
with Symbols_io;

package body HexaDobl_Complex_Solutions_io is

-- INPUT OF A SOLUTION VECTOR :

  procedure get_vector ( s : in out Solution ) is
  begin
    get_vector(Standard_Input,s);
  end get_vector;

  procedure get_vector ( file : in file_type; s : in out Solution ) is

    ind : integer32;

  begin
    if Symbol_Table.Number < natural32(s.n) then
      Symbol_Table.Clear;
      Symbol_Table.Init(natural32(s.n));
      for i in s.v'range loop
        declare
          sb : constant Symbol := Symbols_io.Read_Symbol(file);
        begin
          Symbol_Table.Add(sb);
          Symbols_io.Skip_Symbol(file); get(file,s.v(i));
        exception
          when others
            => put("exception raised at component "); put(i,1); new_line;
               raise;
        end;
      end loop;
    else
      for i in s.v'range loop
        declare
        begin
          ind := integer32(Symbols_io.Get_Symbol(file));
          Symbols_io.Skip_Symbol(file); get(file,s.v(ind));
        exception
          when others
            => put("exception raised at component "); put(i,1); new_line;
               raise;
        end;
      end loop;
    end if; 
  end get_vector;

-- INPUT THE DIAGNOSTICS BANNER :

  procedure Scan_Diagnostics
              ( file : in file_type; err,rco,res : out hexa_Double ) is

  -- DESCRIPTION :
  --   Scans the file for the three diagnostic fields produced by the
  --   root refiner.

    found : boolean;
    d_err,d_rco,d_res : double_float := 0.0;

  begin
    Scan_Line(file,"err :",found);
    if not found then
      err := Create(integer(0));
      rco := Create(integer(0));
      res := Create(integer(0));
    else
      get(file,d_err); err := Create(d_err);
      Scan_Line(file,"rco :",found);
      if not found then
        rco := Create(integer(0));
        res := Create(integer(0));
      else
        get(file,d_rco); rco := Create(d_rco);
        Scan_Line(file,"res :",found);
        if not found
         then res := Create(integer(0));
         else get(file,d_res); res := Create(d_res);
        end if;
      end if;
    end if;
  end Scan_Diagnostics;

-- INPUT OF A SOLUTION :

  procedure get ( s : out Solution ) is
  begin
    get(Standard_Input,s);
  end get;

  procedure get ( n : in natural32; ls : out Link_to_Solution ) is
  begin
    get(Standard_Input,n,ls);
  end get;

  procedure get ( file : in file_type; s : out Solution ) is

    c : character;

  begin
    get(file,c); Standard_Complex_Solutions_io.write_warning('t',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(':',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
    get(file,s.t); skip_line(file);
    declare
    begin
      get(file,c); Standard_Complex_Solutions_io.write_warning('m',c);
      get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
      get(file,c); Standard_Complex_Solutions_io.write_warning(':',c);
      get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
      get(file,s.m);
    exception
      when others => new_line; 
                     put_line("exception in reading m"); raise;
    end;
    if not End_of_Line(file)
     then get(file,c); Skip_line(file);  -- skip information on this line
    end if;
    get(file,c); skip_line(file);
    get_vector(file,s);
    skip_line(file);
    Scan_Diagnostics(file,s.err,s.rco,s.res);
  exception
    when others => new_line; put_line("exception in get of solution"); raise;
  end get;

  procedure get ( file : in file_type;
                  n : in natural32; ls : out Link_to_Solution ) is

    s : Solution(integer32(n));
    c : character;

  begin
    get(file,c); Standard_Complex_Solutions_io.write_warning('t',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(':',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
    get(file,s.t); skip_line(file);
    get(file,c); Standard_Complex_Solutions_io.write_warning('m',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(':',c);
    get(file,c); Standard_Complex_Solutions_io.write_warning(' ',c);
    get(file,s.m);
    if not End_of_Line(file)
     then get(file,c); Skip_line(file);  -- skip information on this line
    end if;
    get(file,c); skip_line(file);
    get_vector(file,s);
    ls := new Solution'(s);
  end get;

-- OUTPUT OF A SOLUTION VECTOR :

  procedure put_vector ( v : in HexaDobl_Complex_Vectors.Vector ) is
  begin
    put_vector(Standard_Output,v);
  end put_vector;

  procedure put_vector ( file : in file_type;
                         v : in HexaDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := v'last;

  begin
    if Symbol_Table.Number < natural32(n) then
      for i in v'range loop
        put(file," x"); put(file,i,1); put(file," : ");
        put(file,v(i)); new_line(file);
      end loop;
    else
      for i in v'range loop
        put(file,' '); 
        Symbols_io.put_symbol(file,natural32(i)); put(file," : ");
        put(file,v(i)); new_line(file);
      end loop;
    end if;
  end put_vector;

  procedure put_vector ( s : in Solution ) is
  begin
    put_vector(Standard_Output,s);
  end put_vector;

  procedure put_vector ( file : in file_type; s : in Solution ) is
  begin
    put_vector(file,s.v);
  end put_vector;

  procedure put_diagnostics ( s : in Solution ) is
  begin
    put_diagnostics(Standard_Output,s);
  end put_diagnostics;

  procedure put_diagnostics ( file : in file_type; s : in Solution ) is
  begin
    put(file,"==");
    put(file," err : "); put(file,s.err,3); put(file," =");
    put(file," rco : "); put(file,s.rco,3); put(file," =");
    put(file," res : "); put(file,s.res,3); put(file," =");
  end put_diagnostics;

-- OUTPUT OF A SOLUTION :

  procedure put ( s : in Solution ) is
  begin
    put(Standard_Output,s);
  end put;

  procedure put ( file : in file_type; s : in Solution ) is
  begin
    put(file,"t : "); put(file,s.t); new_line(file);
    put(file,"m : "); put(file,s.m,1); new_line(file);
    put_line(file,"the solution for t :");
    put_vector(file,s);
    put(file,"==");
    put(file," err : "); put(file,s.err,3); put(file," =");
    put(file," rco : "); put(file,s.rco,3); put(file," =");
    put(file," res : "); put(file,s.res,3); put(file," =");
  end put;

-- INPUT OF A LIST OF SOLUTIONS :

  procedure get ( len,n : in natural32;
                  sols,sols_last : in out Solution_List ) is
  begin
    get(Standard_Input,sols,sols_last);
  end get;

  procedure get ( len,n : in natural32; sols : in out Solution_List ) is
  begin
    get(Standard_Input,len,n,sols);
  end get;

  procedure get ( sols,sols_last : in out Solution_List ) is
  begin
    get(Standard_Input,sols,sols_last);
  end get;

  procedure get ( sols : in out Solution_List ) is
  begin
    get(Standard_Input,sols);
  end get;

  procedure get ( file : in file_type; len,n : in natural32;
                  sols,sols_last : in out Solution_List ) is

    c : character;
    i : natural32 := 1;

  begin
    while i <= len loop
      get(file,c); skip_line(file);    -- skip opening bar
      get(file,c); skip_line(file);    -- skip line with solution number
      declare
        s : Solution(integer32(n));
      begin
        get(file,s);
        Append(sols,sols_last,s);
      exception
        when others => new_line;
                       put_line("exception raised in inner get"); raise;
      end;
      i := i+1;
    end loop;
    get(file,c); skip_line(file);     -- skip closing bar
  exception
    when others
      => new_line;
         put("Exception raised while reading solution ");
         put(i,1); put_line(".");
         raise;
  end get;

  procedure get ( file : in file_type; len,n : in natural32;
                  sols : in out Solution_List ) is

    sols_last : Solution_List;

  begin
    get(file,len,n,sols,sols_last);
  end get;

  procedure get ( file : in file_type;
                  sols,sols_last : in out Solution_List ) is

    len,n : natural32 := 0;

  begin
    get(file,len); get(file,n);
    get(file,len,n,sols,sols_last);
  exception
    when others
      => new_line;
         put("Exception raised when reading solutions,");
         put(" length = "); put(len,1); 
         put(", dimension = "); put(n,1); put_line(".");
         raise;
  end get;

  procedure get ( file : in file_type; sols : in out Solution_List ) is

    sols_last : Solution_List;

  begin
    get(file,sols,sols_last);
  end get;

-- OUTPUT OF A LIST OF SOLUTIONS :

  procedure put ( sols : in Solution_List ) is
  begin
    put(Standard_Output,sols);
  end put;

  procedure put ( len,n : in natural32; sols : in Solution_List ) is
  begin
    put(Standard_Output,len,n,sols);
  end put;

  procedure put ( file : in file_type; sols : in Solution_List ) is
  begin
    if not Is_Null(sols) then
      declare
        count : natural32 := 1;
        temp : Solution_List := sols;
      begin
        Standard_Complex_Solutions_io.put_bar(file);
        while not Is_Null(temp) loop
          put(file,"solution "); put(file,count,1);
          put(file," :"); new_line(file);
          put(file,Head_Of(temp).all);
          put_line(file,"="); 
          temp := Tail_Of(temp);
          count := count + 1;
        end loop;
      end;
    end if;
  end put;

  procedure put ( file : in file_type; len,n : in natural32;
                  sols : in Solution_List ) is
  begin
    put(file,len,1); put(file," "); put(file,n,1); new_line(file);
    put(file,sols);
  end put;

  procedure Display_Format is

    s : array(1..24) of string(1..65);

  begin
    s( 1):="  A solution list of a complex polynomial system  is  denoted  by";
    s( 2):="the  number of solutions and the dimension, followed by a list of";
    s( 3):="solutions.   The  solutions  are  separated  by  a  banner  line,";
    s( 4):="followed by their position in the list.                          ";
    s( 5):="  A solution consists of the current value  of  the  continuation";
    s( 6):="parameter  t,  its  multiplicity  (or  winding number) m, and the";
    s( 7):="solution vector.                                                 ";
    s( 8):="  A solution vector contains as many lines as the dimension.  The";
    s( 9):="i-th  line  starts  with  the  symbol  that  represents  the i-th";
    s(10):="unknown, followed by the colon `:' and two floating-point numbers";
    s(11):="representing  respectively  the  real  and  imaginary part of the";
    s(12):="solution component.                                              ";
    s(13):="  As example we list the solution  list of  the  regular solution";
    s(14):="(1,i) of a 2-dimensional system in the unknowns x and y at t=1.  ";
    s(15):="                                                                 ";
    s(16):="1 2                                                              ";
    s(17):="=================================================================";
    s(18):="solution 1 :                                                     ";
    s(19):="t :  1.00000000000000E+00  0.00000000000000E+00                  ";
    s(20):="m : 1                                                            ";
    s(21):="the solution for t :                                             ";
    s(22):=" x : 1.00000000000000E+00  0.00000000000000E+00                  ";
    s(23):=" y : 0.00000000000000E+00  1.00000000000000E+00                  ";
    s(24):="=================================================================";
    for i in s'range loop
      put_line(s(i));
    end loop;
  end Display_Format;

  procedure Read ( sols : in out Solution_List ) is

    file : file_type;

  begin
    put_line("Reading the name of the file for the solutions.");
    Read_Name_and_Open_File(file);
    get(file,sols);
    Close(file);
  exception
    when others => Close(file); Clear(sols);
                   put_line("INCORRECT FORMAT OF SOLUTION LIST");
                   Display_Format; new_line;
                   Read(sols);
  end Read;

  procedure Write ( file : in file_type; sols : in Solution_List ) is
  begin
    if not Is_Null(sols)
     then put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write;

-- INCREMENTAL READ and WRITE (for huge solution lists):

  procedure Read_First ( file : in file_type; len,dim : out natural32 ) is
  begin
    len := 0; dim := 0;
    get(file,len); get(file,dim);
  exception
    when others => new_line;
                   put("Exception raised when reading dimensions.");
                   raise;
  end Read_First;

  procedure Read_First ( file : in file_type; len,dim : out natural32;
                         s : out Link_to_Solution ) is

    c : character;
    
  begin
    Read_First(file,len,dim);
    get(file,c); skip_line(file);    -- skip opening bar
    get(file,c); skip_line(file);    -- skip line with solution number
    declare
      sol : Solution(integer32(dim));
    begin
      get(file,sol);
      s := new Solution'(sol);
    exception 
      when others
        => new_line;
           put_line("Exception raised while reading first solution.");
           raise;
    end;
  end Read_First;

  procedure Read_First ( file : in file_type; n : in natural32;
                         len,dim,nbr : out natural32;
                         first,last : in out Solution_List ) is

    c : character;
    
  begin
    nbr := 0;
    Read_First(file,len,dim);
    while (nbr <= len) and (nbr < n) loop
      nbr := nbr+1;
      get(file,c); skip_line(file);    -- skip opening bar
      get(file,c); skip_line(file);    -- skip line with solution number
      declare
        s : Solution(integer32(dim));
      begin
        get(file,s);
        Append(first,last,s);
      end;
    end loop;
  exception
    when others
      => new_line;
         put("Exception raised while reading solution ");
         put(nbr,1); put_line(".");
         raise;
  end Read_First;

  procedure Write_First ( file : in file_type; len,dim : in natural32 ) is
  begin
    put(file,len,1); put(file," "); put(file,dim,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
  end Write_First;

  procedure Write_First ( file : in file_type; n,len,dim : in natural32;
                          nbr : out natural32; sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    Write_First(file,len,dim);
    nbr := 0;
    while not Is_Null(tmp) and (nbr < n) loop
      ls := Head_Of(tmp);
      nbr := nbr + 1;
      put(file,"solution "); put(file,nbr,1);
      put(file," :"); new_line(file);
      put(file,ls.all);
      put_line(file,"=");
      tmp := Tail_Of(tmp);
    end loop;
  end Write_First;

  procedure Write_First ( file : in file_type; len,dim : in natural32;
                          s : in Solution ) is
  begin
    Write_First(file,len,dim);
    put_line(file,"solution 1 :");
    put(file,s);
    put_line(file,"=");
  end Write_First;

  procedure Write_First ( file : in file_type; len,dim : in natural32;
                          s : in Link_to_Solution ) is
  begin
    Write_First(file,len,dim,s.all);
  end Write_First;

  procedure Read_Next ( file : in file_type; s : out Solution ) is

    c : character;

  begin
    get(file,c); skip_line(file);    -- skip opening bar
    get(file,c); skip_line(file);    -- skip line with solution number
    get(file,s);
  exception
    when others
      => new_line;
         put_line("Exception raised while reading next solution.");
         raise;
  end Read_Next;

  procedure Read_Next ( file : in file_type; dim : in natural32;
                        s : out Link_to_Solution ) is

    sol : Solution(integer32(dim));

  begin
    Read_Next(file,sol);
    s := new Solution'(sol);
  end Read_Next;

  procedure Read_Next ( file : in file_type; n,dim : in natural32;
                        nbr : out natural32;
                        first,last : in out Solution_List ) is

    c : character;

  begin
    nbr := 0;
    while (nbr < n) loop
      nbr := nbr+1;
      get(file,c); skip_line(file);    -- skip opening bar
      get(file,c); skip_line(file);    -- skip line with solution number
      declare
        s : Solution(integer32(dim));
      begin
        get(file,s);
        Append(first,last,s);
      end;
    end loop;
  exception
    when others
      => new_line;
         put("Exception raised while reading solution ");
         put(nbr,1); put_line(".");
         raise;
  end Read_Next;

  procedure Write_Next ( file : in file_type; cnt : in out natural32;
                         s : in Solution ) is

  begin
    cnt := cnt + 1;
    put(file,"solution "); put(file,cnt,1);
    put(file," :"); new_line(file);
    put(file,s);
    put_line(file,"="); 
  end Write_Next;

  procedure Write_Next ( file : in file_type; cnt : in out natural32;
                         s : in Link_to_Solution ) is
  begin
    Write_Next(file,cnt,s.all);
  end Write_Next;

  procedure Write_Next ( file : in file_type; n : in natural32;
                         nbr : out natural32; cnt : in out natural32;
                         sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    nbr := 0;
    while not Is_Null(tmp) and (nbr < n) loop
      ls := Head_Of(tmp);
      cnt := cnt + 1;
      put(file,"solution "); put(file,cnt,1);
      put(file," :"); new_line(file);
      put(file,ls.all);
      put_line(file,"=");
      tmp := Tail_Of(tmp);
      nbr := nbr + 1;
    end loop;
  end Write_Next;

end HexaDobl_Complex_Solutions_io;
