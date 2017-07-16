with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Symbol_Table;
with Standard_Complex_Laur_Strings;      use Standard_Complex_Laur_Strings;

package body Parse_Dimensions is

  p : Poly;
  s : Link_to_Laur_Sys;

  function Dim ( maxvar : natural32; strpol : string ) return natural32 is

    res : natural32;

  begin
    Symbol_Table.Init(maxvar);
    p := Parse(maxvar,strpol);
    res := Symbol_Table.Number; -- faster than Size_of_Support(p)
    return res;
  end Dim;

  function Dim ( maxvar : natural32; strsys : Array_of_Strings )
               return natural32 is

    res : natural32;

  begin
    s := new Laur_Sys(integer32(strsys'first)..integer32(strsys'last));
    Symbol_Table.Init(maxvar);
    for i in strsys'range loop
      s(integer32(i)) := Parse(maxvar,strsys(i).all);
    end loop;
    res := Symbol_Table.Number;
    return res;
  end Dim;

  function Get return Poly is
  begin
    return p;
  end Get;

  function Get return Link_to_Laur_Sys is
  begin
    return s;
  end Get;

  procedure Clear is
  begin
   -- Symbol_Table.Clear;
    Clear(p);
    Clear(s);
  end Clear;

end Parse_Dimensions;
