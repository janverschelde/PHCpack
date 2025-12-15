with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Complex_Laur_Functions;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Laurentials_IO;   use Standard_Complex_Laurentials_IO;
with Standard_Complex_Laur_Systems_IO;  use Standard_Complex_Laur_Systems_IO;

package body Laurent_Homotopy_Derivatives is

  function Eval ( dim : integer32;
                  cff : Standard_Complex_Vectors.Vector;
                  ctp : Standard_Floating_Vectors.Vector;
                  deg : Standard_Integer_VecVecs.VecVec;
                  t : double_float )
                return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;
    trm : Term;
    ldg : Standard_Integer_Vectors.Link_to_Vector;
    pwt : double_float;

  begin
    trm.dg := new Standard_Integer_Vectors.Vector(1..dim);
    for i in cff'range loop
      pwt := t**ctp(i);
      trm.cf := cff(i)*pwt;
      ldg := deg(i);
      for j in ldg'range loop
        trm.dg(j) := ldg(j);
      end loop;
      Add(res,trm);
    end loop;
    Clear(trm);
    return res;
  end Eval;

  function Diff ( p : Standard_Complex_Laurentials.Poly;
                  idx : Standard_Integer_Vectors.Vector;
                  vrblvl : integer32 := 0 )
                return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res,wrk : Poly;
    done : boolean := false;

  begin
    if vrblvl > 0 then
      put_line("-> in Laurent_Homotopy_Derivatives.diff, p = ");
      put_line(p);
      put("idx :"); put(idx); new_line;
    end if;
    Copy(p,wrk);
    for i in idx'range loop
      if vrblvl > 0 then
        put("i = "); put(i,1); put(", idx("); put(i,1); put(") = ");
        put(idx(i),1); new_line;
      end if;
      for j in 1..idx(i) loop
        res := Diff(wrk,i);
        if vrblvl > 0 then
          put_line("p :"); put_line(wrk);
          put("The derivative of p w.r.t. i = "); put(i,1); put(" :");
          put_line(res);
        end if;
        if Equal(res,Null_Poly)
         then done := true;
        end if;
        exit when done;
        Copy(res,wrk);
      end loop;
      exit when done; -- only exit when a derivative became zero
    end loop;
    if vrblvl > 0 then
      put_line("the derivative :"); put_line(res);
    end if;
    Clear(wrk);
    return res;
  end Diff;

  function Eval ( cff : Standard_Complex_VecVecs.VecVec;
                  ctp : Standard_Floating_VecVecs.VecVec;
                  deg : Standard_Integer_VecVecs.Array_of_VecVecs;
                  t : double_float )
                return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(cff'range);

  begin
    for i in res'range loop
      res(i) := Eval(cff'last,cff(i).all,ctp(i).all,deg(i).all,t);
    end loop;
    return res;
  end Eval;

  function Diff ( cff : Standard_Complex_VecVecs.VecVec;
                  ctp : Standard_Floating_VecVecs.VecVec;
                  deg : Standard_Integer_VecVecs.Array_of_VecVecs;
                  idx : Standard_Integer_Vectors.Vector;
                  z : Standard_Complex_Vectors.Vector; t : double_float;
                  vrblvl : integer32 := 0 )
                return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;

    res : Standard_Complex_Vectors.Vector(cff'range);
    htp : Laur_Sys(cff'range) := Eval(cff,ctp,deg,t);
    dfp : Poly;

  begin
    if vrblvl > 0 then
      put_line("-> in Laurent_Homotopy_Derivatives.diff ...");
      put("Homotopy evaluated at t ="); put(t); put_line(" :");
      put_line(htp);
    end if;
    for i in htp'range loop
      dfp := Diff(htp(i),idx,vrblvl-1);
      if Equal(dfp,Null_Poly) 
       then res(i) := create(0.0);
       else res(i) := Standard_Complex_Laur_Functions.Eval(dfp,z);
      end if;
      Clear(dfp);
    end loop;
    Clear(htp);
    return res;
  end Diff;

end Laurent_Homotopy_Derivatives;
