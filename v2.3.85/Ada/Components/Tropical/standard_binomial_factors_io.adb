with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Symbol_Table;
with Standard_Power_Transformations;     use Standard_Power_Transformations;

package body Standard_Binomial_Factors_io is

  procedure Initialize_Symbol_Table ( s : in string ) is

    sb : Symbol_Table.Symbol;

  begin
    sb := (sb'range => ' ');
    sb(s'range) := s;
    Symbol_Table.Clear;
    Symbol_Table.Init(1);
    Symbol_Table.Add(sb);
  end Initialize_Symbol_Table;

  procedure Initialize_Symbol_Table ( s1,s2 : in string ) is

    sb1,sb2 : Symbol_Table.Symbol;

  begin
    sb1 := (sb1'range => ' '); sb1(s1'range) := s1;
    sb2 := (sb2'range => ' '); sb2(s2'range) := s2;
    Symbol_Table.Clear;
    Symbol_Table.Init(2);
    Symbol_Table.Add(sb1); Symbol_Table.Add(sb2);
  end Initialize_Symbol_Table;

  procedure put ( t : in Term ) is
  begin
    put(standard_output,t);
  end put;

  procedure put ( file : in file_type; t : in Term ) is

    k : constant integer32 := Pivot(t.dg.all);
    m : constant Matrix(1..2,1..2) := Eliminate(t.dg.all,k);
    cx,cy : Complex_Number;

  begin
    if k = 2 then
      put(file,"[ x = "); put(file,t.cf);
      put(file,", y = t^"); put(file,t.dg(2),1); put_line(file," ]");
    else
      if m(2,1) = 0 then
        put(file,"[ x = t^"); put(file,t.dg(1),1);
      else
        cx := t.cf**integer(m(2,1));
        put(file,"[ x = ("); put(file,cx); put(file,"*i)*t^");
        put(file,t.dg(1),1); 
      end if;
      put(file,", ");
      if  m(2,2) = 0 then
        put(file,"y = t^");
      else
        if m(2,1) /= 0
         then new_line(file);
        end if;
        cy := t.cf**integer(m(2,2));
        put(file,"y = ("); put(file,cy); put(file,"*i)*t^");
      end if;
      put(file,t.dg(2),1); put_line(file," ]");
    end if;
  end put;

  procedure put ( t : in List_of_Terms ) is
  begin
    put(standard_output,t);
  end put;

  procedure put ( file : in file_type; t : in List_of_Terms ) is

    p : List_of_Terms := t;

  begin
    while not Is_Null(p) loop
      put(file,Head_Of(p));
      p := Tail_Of(p);
    end loop;
  end put;

end Standard_Binomial_Factors_io;
