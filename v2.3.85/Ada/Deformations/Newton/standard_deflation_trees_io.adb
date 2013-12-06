with Characters_and_Numbers;             use Characters_and_Numbers;
with Symbol_Table;                       use Symbol_Table;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;

package body Standard_Deflation_Trees_io is

  procedure Append_Number ( sb : in out Symbol; ind : in out natural32;
                            nb : in natural32 ) is

  -- DESCRIPTION :
  --   Adds the character representation of the number nb to the
  --   symbol sb, starting at the position given by ind.
  --   The position index ind is updated by this append.

    nb_str : constant string := Convert(integer32(nb));

  begin
    for i in nb_str'range loop
      ind := ind + 1;
      sb(integer(ind)) := nb_str(i);
    end loop;
  end Append_Number;

  procedure Add_Multiplier_Symbols ( k,n : in natural32 ) is

    sb : Symbol;
    ind : natural32;

  begin
    sb := (sb'range => ' ');
    sb(1) := 'l';
    sb(2) := 'm';
    sb(3) := '[';
    ind := 3;
    Append_Number(sb,ind,k);
    ind := ind + 1;
    sb(integer(ind)) := ',';
    Symbol_Table.Enlarge(n);
    for i in 1..n loop
      declare
        new_sb : Symbol;
        new_ind : integer;
      begin
        for j in new_sb'range loop
          new_sb(j) := sb(j);
        end loop;
        new_ind := integer(ind);
        Append_Number(new_sb,natural32(new_ind),i);
        new_ind := new_ind + 1;
        new_sb(new_ind) := ']';
        Symbol_Table.Add(new_sb);
      end;
    end loop;
  end Add_Multiplier_Symbols;

  procedure Write_Deflated_System
              ( file : in file_type; p,dp : in Poly_Sys ) is

    ne : constant natural32 := natural32(dp'last);
    nv : constant natural32 := Number_of_Unknowns(dp(dp'first));

  begin
    put(file,ne,nv,dp(p'range));
    for i in p'last+1..dp'last loop
      put_line(file,dp(i));
    end loop;
  end Write_Deflated_System;

  procedure Write ( file : in file_type; t : in Node ) is

    nb : natural32;

  begin
    put(file,"Level "); put(file,t.d,1); put_line(file," :");
    put(file,t.s);
    for i in t.children'range loop
      if t.children(i) /= null then
        for j in 1..(t.d+1)*2 loop 
          put(file," ");
        end loop;
        put(file,"Child "); put(file,i,1); put_line(file," :");
        Add_Multiplier_Symbols(t.d+1,natural32(i));
        Write(file,t.children(i).all);
        nb := Symbol_Table.Number;
        for j in nb-natural32(i)+1..nb loop
          Symbol_Table.Remove(j);
        end loop;
      end if;
    end loop;
  end Write;

  procedure Write_Nodes ( file : in file_type; name : in string;
                          nd : in Node; p : in Poly_Sys ) is

    nb : natural32;
    ft : file_type;

  begin
    if not Is_Null(nd.sols) then
      if nd.d = 0 then
        declare
          fn : constant string := name & "_d0";
        begin
          put_line(file,"See the file " & fn);
          create(ft,out_file,fn);
          put(ft,natural32(nd.s'last),nd.s);
        end;
      else
        put_line(file,"See the file " & name);
        create(ft,out_file,name);
        Write_Deflated_System(ft,p,nd.s);
      end if;
      new_line(ft);
      put(ft,"THE SOLUTIONS OF DEFLATION LEVEL ");
      put(ft,nd.d,1); put_line(ft," :");
      put(ft,Length_Of(nd.sols),natural32(Head_Of(nd.sols).n),nd.sols);
    end if;
    for i in nd.children'range loop
      if nd.children(i) /= null then
        Add_Multiplier_Symbols(nd.d+1,natural32(i));
        declare
          fn : constant string
             := name & "_d" & Convert(integer32(nd.d+1)) & "R"  & Convert(i-1);
        begin
          Write_Nodes(file,fn,nd.children(i).all,p);
        end;
        nb := Symbol_Table.Number;
        for j in reverse nb-natural32(i)+1..nb loop
          Symbol_Table.Remove(j);
        end loop;
      end if;
    end loop;
  end Write_Nodes;

  procedure Write ( file : in file_type; name : in string; t : in Node ) is
  begin
    new_line(file);
    Write_Nodes(file,name,t,t.s);
  end Write;

end Standard_Deflation_Trees_io;
