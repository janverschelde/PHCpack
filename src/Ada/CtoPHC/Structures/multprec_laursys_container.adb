with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Standard_Integer_Vectors;          use Standard_Integer_Vectors;

package body Multprec_LaurSys_Container is

-- INTERNAL DATA :

  lp : Link_to_Laur_Sys;
  ep : Link_to_Eval_Laur_Sys;
  jm : Link_to_Jaco_Mat;
  ej : Link_to_Eval_Jaco_Mat;

-- CREATORS :

  procedure Initialize ( p : in Laur_Sys ) is
  begin
    lp := new Laur_Sys(p'range);
    for i in p'range loop
      Copy(p(i),lp(i));
    end loop;
  end Initialize;

  procedure Initialize ( n : in integer32 ) is
  begin
    lp := new Laur_Sys'(1..n => Null_Poly);
  end Initialize;

  procedure Create_Evaluator is
  begin
    if lp /= null
     then ep := new Eval_Laur_Sys'(Create(lp.all));
    end if;
  end Create_Evaluator;

  procedure Create_Jacobian_Matrix is
  begin
    if lp /= null
     then jm := new Jaco_Mat'(Create(lp.all));
    end if;
  end Create_Jacobian_Matrix;

  procedure Create_Jacobian_Evaluator is
  begin
    if jm = null
     then Create_Jacobian_Matrix;
    end if;
    if jm /= null
     then ej := new Eval_Jaco_Mat'(Create(jm.all));
    end if;
  end Create_Jacobian_Evaluator;

-- CONSTRUCTORS :

  procedure Add_Term ( k : in integer32; t : in Term ) is
  begin
    Add(lp(k),t);
  end Add_Term;

  procedure Add_Poly ( k : in integer32; p : in Poly ) is
  begin
    Add(lp(k),p);
  end Add_Poly;

-- SELECTORS :

  function Dimension return natural32 is
  begin
    if lp = null 
     then return 0;
     else return natural32(lp'last);
    end if;
  end Dimension;

  function Degree ( k : integer32 ) return integer32 is
  begin
    if lp = null then
      return -1;
    elsif lp(k) = Null_Poly then
      return -1;
    else
      return Degree(lp(k));
    end if;
  end Degree;

  function Number_of_Terms ( k : integer32 ) return natural32 is
  begin
    if lp = null
     then return 0;
     else return Number_of_Terms(lp(k));
    end if;
  end Number_of_Terms;

  function Retrieve_Term ( k : integer32; i : natural32 ) return Term is

    res : Term;
    cnt : natural32 := 0;

    procedure Find_Term ( t : in Term; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      if cnt = i then
        Copy(t.cf,res.cf);
        res.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
        continue := false;
      else
        continue := true;
      end if;
    end Find_Term;
    procedure Scan_Terms is new Visiting_Iterator(Find_Term);

  begin
    res.cf := Multprec_Complex_Numbers.Create(natural32(0));
    if lp = null then
      return res;
    elsif i = 0 or i > Number_of_Terms(k) then
      return res;
    else
      Scan_Terms(lp(k));
    end if;
    return res;
  end Retrieve_Term;

  function Retrieve_Poly ( k : integer32 ) return Poly is
  begin
    if lp = null then
      return Null_Poly;
    elsif k = 0 or k > lp'last then
      return Null_Poly;
    else
      return lp(k);
    end if;
  end Retrieve_Poly;

  function Retrieve return Link_to_Laur_Sys is
  begin
    return lp;
  end Retrieve;

  function Evaluator return Link_to_Eval_Laur_Sys is
  begin
    return ep;
  end Evaluator;

  function Jacobian_Matrix return Link_to_Jaco_Mat is
  begin
    return jm;
  end Jacobian_Matrix;

  function Jacobian_Evaluator return Link_to_Eval_Jaco_Mat is
  begin
    return ej;
  end Jacobian_Evaluator;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Clear(lp);
    Clear(ep);
    Clear(jm);
    Clear(ej);
  end Clear;

begin
  lp := null;
  ep := null;
  jm := null;
  ej := null;
end Multprec_LaurSys_Container;
