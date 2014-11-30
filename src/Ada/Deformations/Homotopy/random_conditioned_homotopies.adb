with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Strings;
with Standard_Complex_Poly_Randomizers;

package body Random_Conditioned_Homotopies is

  function Strip_Semicolon ( f : string ) return string is
  begin
    if f(f'last) = ';'
     then return f(f'first..f'last-1);
     else return f;
    end if;
  end Strip_Semicolon;

  function Random_Coefficient_System
             ( f : String_Splitters.Array_of_Strings )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    dim : constant integer32 := integer32(f'last);
    nvr : constant natural32 := natural32(dim);
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..dim);
    pls : Standard_Complex_Poly_Systems.Poly_Sys(1..dim);

  begin
    if Symbol_Table.Number < nvr then
      Symbol_Table.Clear;
      Symbol_Table.Init(nvr);
    end if;
    pls := Standard_Complex_Poly_Strings.Parse(nvr,f);
    res := Standard_Complex_Poly_Randomizers.Complex_Randomize1(pls);
    Standard_Complex_Poly_Systems.Clear(pls);
    return res;
   end Random_Coefficient_System;

  function Conditioned_Homotopy
             ( f : String_Splitters.Array_of_Strings )
             return String_Splitters.Array_of_Strings is

    dim : constant integer32 := integer32(f'last);
    res : String_Splitters.Array_of_Strings(f'range);
    target : Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
           := Random_Coefficient_System(f);
    f_one : String_Splitters.Array_of_Strings(f'range)
          := Standard_Complex_Poly_Strings.Write(target);
    starts : Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
           := Random_Coefficient_System(f);
    f_zero : String_Splitters.Array_of_Strings(f'range)
           := Standard_Complex_Poly_Strings.Write(starts);
    c_zero : constant string := "2*t^2 - 3*t + 1";
    c_half : constant string := "-4*t^2 + 4*t";
    c_one : constant string := "2*t^2 - t";

  begin
    for i in res'range loop
      declare
        s : constant string 
          := "(" & c_zero & ")*(" & Strip_Semicolon(f_zero(i).all) & ") + "
           & "(" & c_half & ")*(" & Strip_Semicolon(f(i).all) & ") + "
           & "(" & c_one & ")*(" & Strip_Semicolon(f_one(i).all) & ");";
      begin
        res(i) := new string'(s);
      end;
    end loop;
    Standard_Complex_Poly_Systems.Clear(target);
    Standard_Complex_Poly_Systems.Clear(starts);
    String_Splitters.Clear(f_zero);
    String_Splitters.Clear(f_one);
    return res;
  end Conditioned_Homotopy;

end Random_Conditioned_Homotopies;
