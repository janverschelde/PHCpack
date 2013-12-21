with text_io;                           use text_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Matrices;         use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Integer_Circuits;         use Standard_Integer_Circuits;
with Standard_Complex_Laur_Functions;   use Standard_Complex_Laur_Functions;

package body Standard_Monomial_Map_Circuits is

  function Equation ( map : Monomial_Map; c : Vector ) return Poly is

    res : Poly;
    lead,et : Term;

  begin
    lead.dg := new Standard_Integer_Vectors.Vector'(c);
    lead.cf := Create(1.0);
    et.dg := new Standard_Integer_Vectors.Vector'(c'range => 0);
    et.cf := Eval(lead,map.c);
    for i in lead.dg'range loop
      if lead.dg(i) < 0 then
        et.dg(i) := -lead.dg(i);   -- clear negative powers
        lead.dg(i) := 0;
      end if;
    end loop;
    res := Create(lead);
    Sub(res,et);
    Clear(lead);
    Clear(et);
    return res;
  end Equation;

  function Equations ( map : Monomial_Map; c : List ) return Laur_Sys is

    res : Laur_Sys(1..integer32(Length_Of(c)));
    tmp : List := c;
    lv : Link_to_Vector;

  begin
    for i in res'range loop
      lv := Head_Of(tmp);
      res(i) := Equation(map,lv.all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Equations;

  function Circuits ( map : Monomial_Map; d : integer32 ) return List is

    A : constant Matrix := Tropism_Configuration(map);
    res : List;

  begin
    put_line("The tropism configuration : "); put(A);
    res := Circuits(A,d);
    return res;
  end Circuits;

end Standard_Monomial_Map_Circuits;
