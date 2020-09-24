with unchecked_deallocation;
with TripDobl_Complex_Poly_SysFun;      use TripDobl_Complex_Poly_SysFun;
with TripDobl_Complex_Jaco_Matrices;    use TripDobl_Complex_Jaco_Matrices;

package body TripDobl_Coefficient_Homotopy is

-- INTERNAL DATA STRUCTURES :

  type homdata ( n : integer32 ) is record
    k : natural32;
    gamma : TripDobl_Complex_Vectors.Vector(1..n);
    htcff : TripDobl_Complex_VecVecs.VecVec(1..n);
    cp,cq : TripDobl_Complex_VecVecs.VecVec(1..n);
    ip,iq : Standard_Integer_VecVecs.VecVec(1..n);
    hf : Eval_Coeff_Poly_Sys(1..n);
    jf : Eval_Coeff_Jaco_Mat(1..n,1..n);
    mf : Mult_Factors(1..n,1..n);
  end record;

  type homtp is access homdata;
  hom : homtp;

-- PART I : for one pair of polynomials 

  function Coefficients ( p : Poly ) return TripDobl_Complex_Vectors.Vector is

    nbt : constant integer32 := integer32(Number_of_Terms(p));
    res : TripDobl_Complex_Vectors.Vector(1..nbt);
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
    zero : constant triple_double := create(integer(0));

    procedure Visit ( t : in Term; continue : out boolean ) is

      rt : Term;
      imcnt : triple_double;

    begin
      cnt := cnt + 1;
      rt.dg := t.dg;
      if real then
        rt.cf := create(cnt);
      else
        imcnt := create(cnt);
        rt.cf := create(zero,imcnt);
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
             ( c : TripDobl_Complex_Vectors.Vector; real : boolean )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(c'range);
    cnt : integer32 := 0;
    idx : integer32;

  begin
    for i in c'range loop
      if real
       then idx := integer32(to_double(REAL_PART(c(i))));
       else idx := integer32(to_double(IMAG_PART(c(i))));
      end if;
      if idx > 0 then
        cnt := cnt + 1;
        res(cnt) := i;
      end if;
    end loop;
    return res(1..cnt);
  end Index_of_Labels;

  procedure Evaluated_Coefficients
             ( cff : in out TripDobl_Complex_Vectors.Vector;
               cp,cq : in TripDobl_Complex_Vectors.Vector;
               ip,iq : in Standard_Integer_Vectors.Vector;
               t : in triple_double ) is

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
             ( cff : in out TripDobl_Complex_Vectors.Vector;
               cp,cq : in TripDobl_Complex_Vectors.Vector;
               ip,iq : in Standard_Integer_Vectors.Vector;
               k : in natural32;
               gamma : in TripDobl_Complex_Vectors.Vector;
               t : in Complex_Number ) is

    idx : integer32;
    pfac : Complex_Number := gamma(gamma'first);
    one : constant triple_double := create(integer(1));
    onemint : constant Complex_Number := one - t;
    powtk1 : Complex_Number := onemint;
    powtk2 : Complex_Number := t;

  begin
    if k > 1 then
      for i in 2..k loop
        powtk1 := onemint*powtk1;
        powtk2 := t*powtk2;
      end loop;
    end if;
    pfac := powtk1*pfac;
    for i in cp'range loop
      cff(ip(i)) := pfac*cp(i);
    end loop;
    for i in cq'range loop
      idx := iq(i);
      cff(idx) := cff(idx) + powtk2*cq(i);
    end loop;
  end Evaluated_Coefficients;

  procedure Evaluated_Coefficients
             ( cff : in TripDobl_Complex_Vectors.Link_to_Vector;
               cp,cq : in TripDobl_Complex_Vectors.Link_to_Vector;
               ip,iq : in Standard_Integer_Vectors.Link_to_Vector;
               t : in triple_double ) is

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
             ( cff : in TripDobl_Complex_Vectors.Link_to_Vector;
               cp,cq : in TripDobl_Complex_Vectors.Link_to_Vector;
               ip,iq : in Standard_Integer_Vectors.Link_to_Vector;
               k : in natural32;
               gamma : in TripDobl_Complex_Vectors.Vector;
               t : in Complex_Number ) is

    idx : integer32;
    pfac : Complex_Number := gamma(gamma'first);
    one : constant triple_double := create(integer(1));
    onemint : constant Complex_Number := one - t;
    powtk1 : Complex_Number := onemint;
    powtk2 : Complex_Number := t;

  begin
    if k > 1 then
      for i in 2..k loop
        powtk1 := onemint*powtk1;
        powtk2 := t*powtk2;
      end loop;
    end if;
    pfac := powtk1*pfac;
    for i in cp'range loop
      cff(ip(i)) := pfac*cp(i);
    end loop;
    for i in cq'range loop
      idx := iq(i);
      cff(idx) := cff(idx) + powtk2*cq(i);
    end loop;
  end Evaluated_Coefficients;

  function Evaluated_Coefficients
             ( nbcff : integer32;
               cp,cq : TripDobl_Complex_Vectors.Vector;
               ip,iq : Standard_Integer_Vectors.Vector;
               t : triple_double ) return TripDobl_Complex_Vectors.Vector is

    res : TripDobl_Complex_Vectors.Vector(1..nbcff);

  begin
    res := (res'range => Create(integer(0)));
    Evaluated_Coefficients(res,cp,cq,ip,iq,t);
    return res;
  end Evaluated_Coefficients;

-- PART II : for a pair of systems

  function Coefficients
             ( p : Poly_Sys ) return TripDobl_Complex_VecVecs.VecVec is

    res : TripDobl_Complex_VecVecs.VecVec(p'range);

  begin
    for i in p'range loop
      declare
        c : constant TripDobl_Complex_Vectors.Vector
          := Coefficients(p(i));
      begin
        res(i) := new TripDobl_Complex_Vectors.Vector'(c);
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
             ( c : TripDobl_Complex_VecVecs.VecVec; real : boolean )
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
             ( cff : in out TripDobl_Complex_VecVecs.VecVec;
               cp,cq : in TripDobl_Complex_VecVecs.VecVec;
               ip,iq : in Standard_Integer_VecVecs.VecVec;
               t : in triple_double ) is

    cffk : TripDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in cff'range loop
      cffk := cff(i);
      for j in cffk'range loop
        cffk(j) := Create(integer(0));
      end loop;
      Evaluated_Coefficients(cffk,cp(i),cq(i),ip(i),iq(i),t);
    end loop;
  end Evaluated_Coefficients;

  procedure Evaluated_Coefficients
             ( cff : in out TripDobl_Complex_VecVecs.VecVec;
               cp,cq : in TripDobl_Complex_VecVecs.VecVec;
               ip,iq : in Standard_Integer_VecVecs.VecVec;
               k : in natural32;
               gamma : in TripDobl_Complex_Vectors.Vector;
               t : in Complex_Number ) is

    cffk : TripDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in cff'range loop
      cffk := cff(i);
      for j in cffk'range loop
        cffk(j) := Create(integer(0));
      end loop;
      Evaluated_Coefficients(cffk,cp(i),cq(i),ip(i),iq(i),k,gamma,t);
    end loop;
  end Evaluated_Coefficients;

-- PART III : encapsulation

  procedure Create ( p,q : in Poly_Sys; k : in natural32;
                     a : in Complex_Number ) is

    n : constant integer32 := p'last;
    lp,lq,lh : Poly_Sys(p'range);
    hd : homdata(n);

  begin
    hd.k := k;
    for i in 1..n loop
      hd.gamma(i) := a;
    end loop;
    hd.cp := Coefficients(p);
    hd.cq := Coefficients(q);
    lp := Labeled_Coefficients(p,true);
    lq := Labeled_Coefficients(q,false);
    lh := lp + lq;
    hd.htcff := Coefficients(lh);
    hd.ip := Index_of_Labels(hd.htcff,true);
    hd.iq := Index_of_Labels(hd.htcff,false);
    hd.hf := Create(lh); 
    TripDobl_Complex_Jaco_Matrices.Create(lh,hd.jf,hd.mf);
    hom := new homdata'(hd);
    Clear(lp); Clear(lq); Clear(lh);
  end Create;

  function Number_of_Equations return integer32 is
  begin
    if hom = null
     then return -1;
     else return hom.n;
    end if;
  end Number_of_Equations;

  function All_Start_Coefficients return TripDobl_Complex_VecVecs.VecVec is
  begin
    return hom.cp; -- p is the start system
  end All_Start_Coefficients;

  function Start_Coefficients
             ( idx : integer32 )
             return TripDobl_Complex_Vectors.Link_to_Vector is
  begin
    return hom.cp(idx);
  end Start_Coefficients;

  function All_Target_Coefficients return TripDobl_Complex_VecVecs.VecVec is
  begin
    return hom.cq; -- q is the target system
  end All_Target_Coefficients;

  function Target_Coefficients
             ( idx : integer32 )
             return TripDobl_Complex_Vectors.Link_to_Vector is
  begin
    return hom.cq(idx);
  end Target_Coefficients;

  function Eval ( x : TripDobl_Complex_Vectors.Vector;
                  t : Complex_Number )
                return TripDobl_Complex_Vectors.Vector is

    res : TripDobl_Complex_Vectors.Vector(x'range);

  begin
    Evaluated_Coefficients
     -- (hom.htcff,hom.cp,hom.cq,hom.ip,hom.iq,REAL_PART(t));
      (hom.htcff,hom.cp,hom.cq,hom.ip,hom.iq,hom.k,hom.gamma,t);
    res := Eval(hom.hf,hom.htcff,x);
    return res;
  end Eval;

  function Diff ( x : TripDobl_Complex_Vectors.Vector;
                  t : Complex_Number ) return Matrix is

    res : Matrix(x'range,x'range);
 
  begin
    Evaluated_Coefficients
     -- (hom.htcff,hom.cp,hom.cq,hom.ip,hom.iq,REAL_PART(t));
      (hom.htcff,hom.cp,hom.cq,hom.ip,hom.iq,hom.k,hom.gamma,t);
    res := Eval(hom.jf,hom.mf,hom.htcff,x);
    return res;
  end Diff;

  procedure free is new unchecked_deallocation (homdata,homtp);

  procedure Clear is
  begin
    if hom /= null then
      Clear(hom.hf);
      Clear(hom.jf);
      Clear(hom.mf);
      TripDobl_Complex_VecVecs.Clear(hom.htcff);
      TripDobl_Complex_VecVecs.Clear(hom.cp);
      TripDobl_Complex_VecVecs.Clear(hom.cq);
      Standard_Integer_VecVecs.Clear(hom.ip);
      Standard_Integer_VecVecs.Clear(hom.iq);
      free(hom);
    end if;
  end Clear;

end TripDobl_Coefficient_Homotopy;
