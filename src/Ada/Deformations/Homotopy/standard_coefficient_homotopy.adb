with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;

package body Standard_Coefficient_Homotopy is

-- PART I : for one pair of polynomials 

  function Coefficients ( p : Poly ) return Standard_Complex_Vectors.Vector is

    nbt : constant integer32 := integer32(Number_of_Terms(p));
    res : Standard_Complex_Vectors.Vector(1..nbt);
    idx : integer32 := 0;

    procedure Visit ( t : in Term; continue : out boolean ) is
    begin
      idx := idx + 1;
      res(idx) := t.cf;
      continue := true;
    end Visit;
    procedure Visit_Terms is new Visiting_Iterator(Visit);

  begin
    Visit_Terms(p);
    return res;
  end Coefficients;

  function Labeled_Coefficients ( p : Poly; real : boolean ) return Poly is

    res : Poly := Null_Poly;
    cnt : integer32 := 0;

    procedure Visit ( t : in Term; continue : out boolean ) is

      rt : Term;

    begin
      cnt := cnt + 1;
      rt.dg := t.dg;
      if real
       then rt.cf := create(cnt);
       else rt.cf := create(0.0,double_float(cnt));
      end if;
      Add(res,rt);
      continue := true;
    end Visit;

    procedure Visit_Terms is new Visiting_Iterator(Visit);

  begin
    Visit_Terms(p);
    return res;
  end Labeled_Coefficients;

  function Index_of_Labels
             ( c : Standard_Complex_Vectors.Vector; real : boolean )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(c'range);
    cnt : integer32 := 0;
    idx : integer32;

  begin
    for i in c'range loop
      if real
       then idx := integer32(REAL_PART(c(i)));
       else idx := integer32(IMAG_PART(c(i)));
      end if;
      if idx > 0 then
        cnt := cnt + 1;
        res(cnt) := i;
      end if;
    end loop;
    return res(1..cnt);
  end Index_of_Labels;

  procedure Evaluated_Coefficients
             ( cff : in out Standard_Complex_Vectors.Vector;
               cp,cq : in Standard_Complex_Vectors.Vector;
               ip,iq : in Standard_Integer_Vectors.Vector;
               t : in double_float ) is

    idx : integer32;

  begin
    for i in cp'range loop
      cff(ip(i)) := (1.0-t)*cp(i);
    end loop;
    for i in cq'range loop
      idx := iq(i);
      cff(idx) := cff(idx) + t*cq(i);
    end loop;
  end Evaluated_Coefficients;

  procedure Evaluated_Coefficients
             ( cff : in out Standard_Complex_Vectors.Link_to_Vector;
               cp,cq : in Standard_Complex_Vectors.Link_to_Vector;
               ip,iq : in Standard_Integer_Vectors.Link_to_Vector;
               t : in double_float ) is

    idx : integer32;

  begin
    for i in cp'range loop
      cff(ip(i)) := (1.0-t)*cp(i);
    end loop;
    for i in cq'range loop
      idx := iq(i);
      cff(idx) := cff(idx) + t*cq(i);
    end loop;
  end Evaluated_Coefficients;

  function Evaluated_Coefficients
             ( nbcff : integer32;
               cp,cq : Standard_Complex_Vectors.Vector;
               ip,iq : Standard_Integer_Vectors.Vector;
               t : double_float ) return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(1..nbcff);

  begin
    res := (res'range => Create(0.0));
    Evaluated_Coefficients(res,cp,cq,ip,iq,t);
    return res;
  end Evaluated_Coefficients;

-- PART II : for a pair of systems

  function Coefficients
             ( p : Poly_Sys ) return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(p'range);

  begin
    for i in p'range loop
      declare
        c : constant Standard_Complex_Vectors.Vector
          := Coefficients(p(i));
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(c);
      end;
    end loop;
    return res;
  end Coefficients;

  function Labeled_Coefficients
             ( p : Poly_Sys; real : boolean ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Labeled_Coefficients(p(i),real);
    end loop;
    return res;
  end Labeled_Coefficients;

  function Index_of_Labels
             ( c : Standard_Complex_VecVecs.VecVec; real : boolean )
             return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(c'range);

  begin
    for i in c'range loop
      declare
        e : constant Standard_Integer_Vectors.Vector
          := Index_of_Labels(c(i).all,real);
      begin
        res(i) := new Standard_Integer_Vectors.Vector'(e);
      end;
    end loop;
    return res;
  end Index_of_Labels;

  procedure Evaluated_Coefficients
             ( cff : in out Standard_Complex_VecVecs.VecVec;
               cp,cq : in Standard_Complex_VecVecs.VecVec;
               ip,iq : in Standard_Integer_VecVecs.VecVec;
               t : in double_float ) is

    cffk : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in cff'range loop
      cffk := cff(i);
      for j in cffk'range loop
        cffk(j) := Create(0.0);
      end loop;
      Evaluated_Coefficients(cffk,cp(i),cq(i),ip(i),iq(i),t);
    end loop;
  end Evaluated_Coefficients;

end Standard_Coefficient_Homotopy;
