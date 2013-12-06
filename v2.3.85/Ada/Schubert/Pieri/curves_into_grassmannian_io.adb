with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Symbol_Table;                       use Symbol_Table;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;

package body Curves_into_Grassmannian_io is

  procedure Set_up_Symbol_Table
              ( m,p,q : in natural32; top,bottom : in Bracket;
                addmix : in natural32 ) is

    n : constant natural32 := Number_of_Variables(top,bottom) + 2 + addmix;
    rws : constant natural32 := (m+p)*(q+1);
    row,deg : natural32;
    sb : Symbol;

  begin
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Symbol_Table.Init(n);             -- initialization with #variables
    sb := (sb'range => ' ');
    sb(1) := 'x';
    sb(4) := 's';
    for j in 1..integer32(p) loop     -- adding the rest columnwise
      row := 0; deg := 0;
      for i in 1..rws loop
        row := row + 1;
        if i >= top(j) and i <= bottom(j) then
          sb(2) := Convert_Hexadecimal(row);
          sb(3) := Convert_Hexadecimal(natural32(j));
          sb(5) := Convert_Hexadecimal(deg);
          Symbol_Table.Add(sb);
        end if;
        if i mod (m+p) = 0
         then row := 0; deg := deg+1;
        end if;
      end loop;
    end loop;
    sb := (sb'range => ' ');
    sb(1) := 's';
    if addmix = 0 then
      Symbol_Table.Add(sb);       -- adding "s"
    else
      sb(2) := '1';
      Symbol_Table.Add(sb);       -- adding "s1"
      sb(2) := '2';
      Symbol_Table.Add(sb);       -- adding "s2"
      sb(2) := ' ';
    end if;
    sb(1) := 't';
    Symbol_Table.Add(sb);             -- adding "t"
  end Set_up_Symbol_Table;

  procedure One_Set_up_Symbol_Table
              ( m,p,q : in natural32; top,bottom : in Bracket ) is
  begin
    Set_up_Symbol_Table(m,p,q,top,bottom,0);
  end One_Set_up_Symbol_Table;

  procedure Two_Set_up_Symbol_Table
              ( m,p,q : in natural32; top,bottom : in Bracket ) is
  begin
    Set_up_Symbol_Table(m,p,q,top,bottom,1);
  end Two_Set_up_Symbol_Table;

  procedure Reduce_Symbols ( top,bottom : in Bracket; locmap : in Matrix ) is

    ind : natural32 := 0;

  begin
    for j in locmap'range(2) loop
      for i in top(j)..bottom(j) loop
        ind := ind + 1;
        if locmap(integer32(i),j) = 1
         then Symbol_Table.Remove(ind-natural32(j)+1);
        end if;
      end loop;
    end loop;
  end Reduce_Symbols;

end Curves_into_Grassmannian_io;
