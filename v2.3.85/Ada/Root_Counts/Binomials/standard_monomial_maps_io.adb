with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Characters_and_Numbers;
with Standard_Write_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Vectors;
with Symbol_Table_io;
with Standard_Complex_Laurentials_io;   use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Binomial_Varieties_io;
with Standard_Monomial_Map_Filters;

package body Standard_Monomial_Maps_io is

-- MANAGEMENT OF SYMBOL TABLE :

  procedure Insert_Parameter_Symbols
               ( d : in natural32;
                 symbols : in Symbol_Table.Array_of_Symbols ) is

    n : constant natural32 := natural32(symbols'last);

  begin
    Symbol_Table.Clear;
    Symbol_Table.Init(n+d);
    for i in 1..d loop
      declare
        sb : Symbol_Table.Symbol;
        nb : constant string := Characters_and_Numbers.nConvert(i);
      begin
        sb(sb'range) := (sb'range => ' ');
        sb(sb'first) := 't';
        for k in nb'range loop
          sb(sb'first+k) := nb(k);
        end loop;
        Symbol_Table.Add(sb);
      end;
    end loop;
    for i in symbols'range loop
      Symbol_Table.Add(symbols(i));
    end loop;
  end Insert_Parameter_Symbols;

  procedure Insert_Parameter_Symbols ( d : in natural32 ) is

    n : constant natural32 := Symbol_Table.Number;
    symbols : constant Symbol_Table.Array_of_Symbols := Symbol_Table.Content;

  begin
    Symbol_Table.Clear;
    Symbol_Table.Init(n+d);
    Insert_Parameter_Symbols(d,symbols);
  end Insert_Parameter_Symbols;

  procedure Check_Parameter_Symbols ( n,d : in natural32 ) is

    nb : constant natural32 := Symbol_Table.Number;

  begin
    if nb = n then
      Insert_Parameter_Symbols(d);
    elsif nb /= n + d then
      declare
        symbols : constant Symbol_Table.Array_of_Symbols
                := Symbol_Table.Content;
        ambient : Symbol_Table.Array_of_Symbols(1..integer32(n));
      begin
        for i in 1..n loop
          ambient(integer32(i)) := symbols(symbols'last-integer32(n-i));
        end loop;
        Insert_Parameter_Symbols(d,ambient);
      end;
    end if;
  end Check_Parameter_Symbols;

  procedure Remove_Parameter_Symbols ( n : in natural32 ) is

    d : constant natural32 := Symbol_Table.Number - n;

  begin
    if d > 0 then
      declare
        symbols : constant Symbol_Table.Array_of_Symbols
                := Symbol_Table.Content;
      begin
        Symbol_Table.Clear;
        Symbol_Table.Init(n);
        for i in 1..n loop
          Symbol_Table.Add(symbols(integer32(d+i)));
        end loop;
      end;
    end if;
  end Remove_Parameter_Symbols;

-- INPUT ROUTINES :

  procedure get ( map : out Link_to_Monomial_Map ) is
  begin
    get(standard_input,map);
  end get;

  procedure get ( file : in file_type; map : out Link_to_Monomial_Map ) is

    p : Link_to_Laur_Sys;
    T : Standard_Integer_Matrices.Link_to_Matrix;
    c : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Symbol_Table.Clear;
    get(file,p);
    Standard_Binomial_Varieties_io.Parse_Binomial_Variety(p.all,T,c);
    declare
      m : Monomial_Map(T'last(1));
      v : Standard_Integer_Vectors.Vector(T'range(2));
    begin
      m.d := T'last(2);
      m.c := c.all;
      for i in m.v'range loop
        for j in T'range(2) loop
          v(j) := T(i,j);
        end loop;
        m.v(i) := new Standard_Integer_Vectors.Vector'(v);
      end loop;
      map := new Monomial_Map'(m);
    end;
    Standard_Integer_Matrices.Clear(T);
    Standard_Complex_Vectors.Clear(c);
  end get;

  procedure get ( maps : out Monomial_Map_Array ) is
  begin
    get(standard_input,maps);
  end get;

  procedure get ( file : in file_type; maps : out Monomial_Map_Array ) is
  begin
    for i in maps'range loop
      declare
      begin
        get(file,maps(i));
      exception
        when others
          => put("exception raised when reading map "); put(i,1); new_line;
             raise;
      end;
    end loop;
  end get;

  procedure get ( maps : out Link_to_Monomial_Map_Array ) is

    n : integer32 := 0;

  begin
    get(n);
    declare
      the_maps : Monomial_Map_Array(1..n);
    begin
      get(the_maps);
      maps := new Monomial_Map_Array'(the_maps);
    end;
  end get;

  procedure get ( file : in file_type;
                  maps : out Link_to_Monomial_Map_Array ) is

    n : integer32 := 0;

  begin
    get(file,n); skip_line(file);
    declare
      the_maps : Monomial_Map_Array(1..n);
    begin
      get(file,the_maps);
      maps := new Monomial_Map_Array'(the_maps);
    end;
  end get;

  procedure get ( maps : out Monomial_Map_List ) is

    n : integer32 := 0;
    maps_last : Monomial_Map_List;

  begin
    get(n);
    for i in 1..n loop
      declare
        link_to_map : Link_to_Monomial_Map;
      begin
        get(link_to_map);
        Append(maps,maps_last,link_to_map.all);
      exception
        when others
          => put("something is wrong with map "); put(i,1); new_line;
             raise;
      end;
    end loop;
  end get;

  procedure get ( file : in file_type; maps : out Monomial_Map_List ) is

    n : integer32 := 0;
    maps_last : Monomial_Map_List;

  begin
    get(file,n); skip_line(file);
    for i in 1..n loop
      declare
        link_to_map : Link_to_Monomial_Map;
      begin
        get(file,link_to_map);
        Append(maps,maps_last,link_to_map.all);
      end;
    end loop;
  end get;

-- OUTPUT ROUTINES :

  procedure put ( map : in Monomial_Map ) is
  begin
    put(standard_output,map);
  end put;

  procedure put ( file : in file_type; map : in Monomial_Map ) is

    cnt : natural32 := 0;
    exp : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Check_Parameter_Symbols(natural32(map.n),natural32(map.d));
    Standard_Binomial_Varieties_io.Write_Header
      (file,natural32(map.n),natural32(map.d));
    for j in 1..map.n loop
      declare
        s : constant Symbol_Table.Symbol
          := Symbol_Table.get(natural32(map.d+j)); 
        first : boolean := true;
        isone : constant boolean := Standard_Monomial_Maps.Is_One(map.c(j));
        bracket : boolean;
      begin
        Symbol_Table_io.put(file,s);
        if not Standard_Monomial_Maps.Is_Zero(map.c(j)) then
          put(file," - ");
          if not isone then
            bracket := Standard_Write_Numbers.Is_Real(map.c(j))
                    or Standard_Write_Numbers.Is_Imag(map.c(j));
            if bracket then put(file,"("); end if;
            Standard_Write_Numbers.Write_Number(file,map.c(j),cnt);
            if bracket then put(file,")"); end if;
            first := false;
          end if;
          exp := map.v(j);
          for i in 1..map.d loop
            if exp(i) /= 0 then
              if first
               then first := false;
               else put(file,"*");
              end if;
              put(file,"t"); put(file,i,1);
              if exp(i) /= 1
               then put(file,"^"); put(file,exp(i),1);
              end if;
            end if;
          end loop;
          if first then
            if isone
             then put(file,"1");
            end if;
          end if;
        end if;
        put_line(file,";");
      end;
    end loop;
  end put;

  procedure put ( map : in Link_to_Monomial_Map ) is
  begin
    put(standard_output,map);
  end put;

  procedure put ( file : in file_type; map : in Link_to_Monomial_Map ) is
  begin
    if map /= null
     then put(file,map.all);
    end if;
  end put;

  procedure put ( maps : in Monomial_Map_Array ) is
  begin
    put(standard_output,maps);
  end put;

  procedure put ( file : in file_type; maps : in Monomial_Map_Array ) is
  begin
    for i in maps'range loop
      put(file,maps(i).all);
    end loop;
  end put;

  procedure put ( maps : in Monomial_Map_List ) is
  begin
    put(standard_output,maps);
  end put;

  procedure put ( file : in file_type; maps : in Monomial_Map_List ) is

    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    put(file,Length_Of(maps),1); new_line(file);
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      put(file,link_to_map);
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( maps : in Array_of_Monomial_Map_Lists ) is
  begin
    put(standard_output,maps);
  end put;

  procedure put ( file : in file_type;
                  maps : in Array_of_Monomial_Map_Lists ) is

    cnt : natural32 := 0;
    tmp : Monomial_Map_List;
    link_to_map : Link_to_Monomial_Map;

  begin
    for i in reverse maps'range loop
      cnt := cnt + Length_Of(maps(i));
    end loop;
    put(file,cnt,1); new_line(file);
    for i in reverse maps'range loop
      tmp := maps(i);
      while not Is_Null(tmp) loop
        link_to_map := Head_Of(tmp);
        put(file,link_to_map);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end put;

-- IDEAL REPRESENTATIONS :

  procedure Show_Ideal ( p : in Laur_Sys; map : in Monomial_Map ) is
  begin
    Show_Ideal(standard_output,p,map);
  end Show_Ideal;

  procedure Show_Ideal ( file : in file_type;
                         p : in Laur_Sys; map : in Monomial_Map ) is

    tol : constant double_float := 1.0E-8;
    first : boolean := true;
    f : Link_to_Laur_Sys
      := Standard_Monomial_Map_Filters.Filter(p,map);

  begin
    Remove_Parameter_Symbols(natural32(map.n));
    put(file,"< ");
    for i in map.c'range loop
      if AbsVal(map.c(i)) < tol then
        if first
         then first := false;
         else put(file," , ");
        end if;
        declare
          sb : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(i));
        begin
          Symbol_Table_io.put(file,sb);
        end;
      end if;
    end loop;
    for i in f'range loop
      if first
       then first := false;
       else put(file," , ");
      end if;
      put(file,f(i));
    end loop;
    put_line(file," >");
    Clear(f);
  end Show_Ideal;

  procedure Show_Ideals ( p : in Laur_Sys; maps : in Monomial_Map_List ) is
  begin
    Show_Ideals(standard_output,p,maps);
  end Show_Ideals;

  procedure Show_Ideals ( file : in file_type;
                          p : in Laur_Sys; maps : in Monomial_Map_List ) is

    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    while not Is_Null(tmp) loop
      link_to_map := Head_Of(tmp);
      Show_Ideal(file,p,link_to_map.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Show_Ideals;

  procedure Show_Ideals ( p : in Laur_Sys;
                          maps : in Array_of_Monomial_Map_Lists ) is
  begin
    Show_Ideals(standard_output,p,maps);
  end Show_Ideals;

  procedure Show_Ideals ( file : in file_type;
                          p : in Laur_Sys;
                          maps : in Array_of_Monomial_Map_Lists ) is
  begin
    for i in reverse maps'range loop
      Show_Ideals(file,p,maps(i));
    end loop;
  end Show_Ideals;

-- APPEND MAPS TO A FILE :

  procedure Append ( name : in string; maps : in Monomial_Map_List ) is

    file : file_type;

  begin
    if not Is_Null(maps) then
      Open_Append_File(file,name);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,maps);
    end if;
  end Append;

  procedure Append ( name : in string;
                     maps : in Array_of_Monomial_Map_Lists ) is

    file : file_type;

  begin
    Open_Append_File(file,name);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,maps);
  end Append;

-- STATISTICS :

  procedure Write_Lengths ( c : in Array_of_Monomial_Map_Lists ) is
  begin
    Write_Lengths(standard_output,c);
  end Write_Lengths;

  procedure Write_Lengths
              ( file : in file_type; c : in Array_of_Monomial_Map_Lists ) is

    len : constant Standard_Natural_Vectors.Vector := Lengths(c);
    sum : constant natural32 := Standard_Natural_Vectors.Sum(len);

  begin
    put(file,"#components per dimension :"); put(file,len);
    put(file," sum : "); put(file,sum,1); new_line(file);
  end Write_Lengths;

  procedure Show_Degrees ( maps : in Monomial_Map_List ) is
  begin
    Show_Degrees(standard_output,maps);
  end Show_Degrees;

  procedure Show_Degrees ( maps : in Monomial_Map_List; sum : out natural32 ) is
  begin
    Show_Degrees(standard_output,maps,sum);
  end Show_Degrees;

  procedure Show_Degrees ( p : in Laur_Sys; maps : in Monomial_Map_List ) is
  begin
    Show_Degrees(standard_output,p,maps);
  end Show_Degrees;

  procedure Show_Degrees ( p : in Laur_Sys;
                           maps : in Monomial_Map_List; sum : out natural32 ) is
  begin
    Show_Degrees(standard_output,p,maps,sum);
  end Show_Degrees;

  procedure Show_Degrees ( file : in file_type;
                           maps : in Monomial_Map_List ) is

    deg : constant Standard_Natural_Vectors.Vector := Degrees(maps);

  begin
    for i in deg'range loop
      put(file,"-> map "); put(file,i,1); 
      put(file," has degree : "); put(file,deg(i),1); new_line(file);
    end loop;
  end Show_Degrees;

  procedure Show_Degrees ( file : in file_type;
                           maps : in Monomial_Map_List; sum : out natural32 ) is

    deg : constant Standard_Natural_Vectors.Vector := Degrees(maps);

  begin
    for i in deg'range loop
      put(file,"-> map "); put(file,i,1); 
      put(file," has degree : "); put(file,deg(i),1); new_line(file);
    end loop;
    sum := Standard_Natural_Vectors.Sum(deg);
    put(file,"Sum of the degrees of the maps : ");
    put(file,sum,1); put_line(file,".");
  end Show_Degrees;

  procedure Show_Degrees ( file : in file_type; p : in Laur_Sys;
                           maps : in Monomial_Map_List ) is

    deg : constant Standard_Natural_Vectors.Vector := Degrees(maps);
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    for i in deg'range loop
      link_to_map := Head_Of(tmp);
      put(file,"-> map "); put(file,i,1); put_line(file," is defined by ");
      Show_Ideal(file,p,link_to_map.all);
      put(file,"   and has degree : "); put(file,deg(i),1); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Show_Degrees;

  procedure Show_Degrees ( file : in file_type; p : in Laur_Sys;
                           maps : in Monomial_Map_List; sum : out natural32 ) is

    deg : constant Standard_Natural_Vectors.Vector := Degrees(maps);
    tmp : Monomial_Map_List := maps;
    link_to_map : Link_to_Monomial_Map;

  begin
    for i in deg'range loop
      link_to_map := Head_Of(tmp);
      put(file,"-> map "); put(file,i,1); put_line(file," is defined by ");
      Show_Ideal(file,p,link_to_map.all);
      put(file,"   and has degree : "); put(file,deg(i),1); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    sum := Standard_Natural_Vectors.Sum(deg);
    put(file,"Sum of the degrees of the maps : ");
    put(file,sum,1); put_line(file,".");
  end Show_Degrees;

  procedure Show_Degrees ( maps : in Array_of_Monomial_Map_Lists ) is
  begin
    Show_Degrees(standard_output,maps);
  end Show_Degrees;

  procedure Show_Degrees ( file : in file_type;
                           maps : in Array_of_Monomial_Map_Lists ) is

    sum,total_degree : natural32 := 0;

  begin
    for i in reverse maps'range loop
      if Length_Of(maps(i)) > 0 then
        put(file,"Degrees of maps of dimension ");
        put(file,i,1); put_line(file," :");
        Show_Degrees(file,maps(i),sum);
        total_degree := total_degree + sum;
      end if;
    end loop;
    if total_degree /= sum then
      put(file,"Sum of all degrees over all dimensions : ");
      put(file,total_degree,1); new_line(file);
    end if;
  end Show_Degrees;

  procedure Show_Degrees ( p : in Laur_Sys;
                           maps : in Array_of_Monomial_Map_Lists ) is
  begin
    Show_Degrees(standard_output,p,maps);
  end Show_Degrees;

  procedure Show_Degrees ( file : in file_type; p : in Laur_Sys;
                           maps : in Array_of_Monomial_Map_Lists ) is

    sum,total_degree : natural32 := 0;

  begin
    for i in reverse maps'range loop
      if Length_Of(maps(i)) > 0 then
        put(file,"Degrees of maps of dimension ");
        put(file,i,1); put_line(file," :");
        Show_Degrees(file,p,maps(i),sum);
        total_degree := total_degree + sum;
      end if;
    end loop;
    if total_degree /= sum then
      put(file,"Sum of all degrees over all dimensions : ");
      put(file,total_degree,1); new_line(file);
    end if;
  end Show_Degrees;

-- READ SYSTEM AND SOLUTIONS :

  procedure Read_System_and_Maps
               ( p : out Link_to_Laur_Sys;
                 maps : out Monomial_Map_List ) is

    file : file_type;

  begin
    put_line("Reading a binomial system ...");
    Read_Name_and_Open_File(file);
    Read_System_and_Maps(file,p,maps);
    close(file);
  end Read_System_and_Maps;

  procedure Read_System_and_Maps
               ( file : in file_type; p : out Link_to_Laur_Sys;
                 maps : out Monomial_Map_List ) is

    found : boolean;

  begin
    get(file,p);
    Scan_and_Skip(file,"THE SOLUTIONS",found);
    if found
     then get(file,maps);
    end if;
  end Read_System_and_Maps;

end Standard_Monomial_Maps_io;
