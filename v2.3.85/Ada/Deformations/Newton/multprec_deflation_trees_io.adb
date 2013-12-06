with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Symbol_Table;                       use Symbol_Table;
with Multprec_Complex_Polynomials;       use Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Standard_Deflation_Trees_io;

package body Multprec_Deflation_Trees_io is

  procedure Write_Deflated_System
              ( file : in file_type; p,dp : in Poly_Sys ) is

    ne : constant integer32 := dp'last - dp'first + 1;
    nv : constant integer32 := integer32(Number_of_Unknowns(dp(dp'first)));

  begin
    put(file,ne,1); put(file,"  "); put(file,nv,1); new_line(file);
    put(file,dp(p'range));
    for i in p'last+1..dp'last loop
      put_line(file,dp(i));
    end loop;
    new_line(file);
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
        Standard_Deflation_Trees_io.Add_Multiplier_Symbols(t.d+1,natural32(i));
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
          put(ft,nd.s'last,1); new_line(ft);
          put(ft,nd.s);
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
        Standard_Deflation_Trees_io.Add_Multiplier_Symbols(nd.d+1,natural32(i));
        declare
          fn : constant string
             := name & "_d" & Convert(integer32(nd.d)+1)
                     & "R"  & Convert(integer32(i)-1);
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

end Multprec_Deflation_Trees_io;
