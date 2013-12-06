with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
--with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
--with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Power_Transformations;
with Standard_Initial_Forms;
with Standard_Binomial_Factors_io;       use Standard_Binomial_Factors_io;

package body Standard_Puiseux_Certificates is

  function Other_Index ( pivot : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   This is an auxiliary routine for inner normals with pivot = 2.
  --   The routine returns the other index, depending on the pivot.
  
  begin
    if pivot = 1
     then return 2;
     else return 1;
    end if;
  end Other_Index;

  function Second_Power ( p : Poly; k,e : integer32 ) return integer32 is

    res : integer32 := 0;
    first_time : boolean := true;
    ind : constant integer32 := Other_Index(k);

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      if t.dg(k) = e then
        if first_time then
          res := t.dg(ind);
          first_time := false;
        elsif t.dg(ind) < res then 
          res := t.dg(ind);
        end if;
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Second_Power;

  function Shift ( p : Poly; k : in integer32 ) return Poly is

    res : Poly;
    e1 : constant integer32 := Minimal_Degree(p,k);
    e2 : constant integer32 := Second_Power(p,k,e1);
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Integer_Vectors.Vector(1..2);
    if k = 1
     then t.dg(1) := -e1; t.dg(2) := -e2;
     else t.dg(2) := -e1; t.dg(1) := -e2;
    end if;
    res := p*t;
    Clear(t);
    return res;
  end Shift;

  function Number_of_Initial_Roots
             ( p : Poly; k : integer32 ) return integer32 is

    res : integer32;
    first_time : boolean := true;
    ind : constant integer32 := Other_Index(k);

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      if t.dg(k) = 0 and t.dg(ind) > 0 then
        if first_time then
          res := t.dg(ind);
          first_time := false;
        elsif t.dg(ind) > res then 
          res := t.dg(ind);
        end if;
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Number_of_Initial_Roots;

  function Initial_Coefficients
              ( p : Poly; k : integer32 )
              return Standard_Complex_Vectors.Vector is

    n : constant integer32 := Number_of_Initial_Roots(p,k);
    res : Standard_Complex_Vectors.Vector(0..n) := (0..n => Create(0.0));
    ind : constant integer32 := Other_Index(k);

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      if t.dg(k) = 0
       then res(t.dg(ind)) := t.cf;
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Initial_Coefficients;

  function Derivative ( p : Standard_Complex_Vectors.Vector;
                        x : Complex_Number ) return Complex_Number is

    res : Complex_Number := double_float(p'last)*p(p'last);
    pow : double_float;

  begin
    for k in reverse p'first+1..p'last-1 loop
      pow := double_float(k);
      res := res*x + pow*p(k);
    end loop;
    return res;
  end Derivative;

  procedure Match ( f,g : in Poly; f1,g1 : in Complex_Number;
                    tol : in double_float; output : in boolean;
                    cff : out Complex_Number; exp : out integer32 ) is

    found : boolean := false;
    cf : constant Standard_Complex_Vectors.Vector := Coefficients(f);
    cg : constant Standard_Complex_Vectors.Vector := Coefficients(g);
    det : Complex_Number;
    i : integer32 := cf'first;
    j : integer32 := cg'first;

  begin
    cff := Create(0.0); exp := 0;
    while i <= cf'last and j <= cg'last loop
      while AbsVal(cf(i)) < tol loop
        i := i + 1;
        exit when i > cf'last;
      end loop;
      exit when i > cf'last;
      while AbsVal(cg(j)) < tol loop
        j := j + 1;
        exit when j > cg'last;
      end loop;
      exit when j > cg'last;
      while i < j and i <= cf'last loop
        i := i + 1;
      end loop;
      exit when i > cf'last;
      while i > j and j <= cg'last loop
        j := j + 1;
      end loop;
      exit when j > cg'last;
      if AbsVal(cf(i)) > tol and AbsVal(cg(j)) > tol then
        det := cf(i)*g1 - cg(j)*f1;
        if AbsVal(det) < tol
         then cff := (-cf(i))/f1; exp := i; found := true;
        end if;
        if output then
          put("Checking "); put(i,1); put(" and "); put(j,1);
          put(" det :"); put(det,3);
          if AbsVal(det) > tol then
            put_line("  no match");
          else
            put_line("  match");
            put("c1 (through f) : "); put(-cf(i)/f1); new_line;
            put("c2 (through g) : "); put(-cg(j)/g1); new_line;
            cff := (-cf(i))/f1; exp := i; found := true;
          end if;
        end if;
      end if;
      exit when found;
      i := i + 1; j := j + 1;
    end loop;
    if output then
      if found
       then put_line("found a matching term");
       else put_line("found no matching term");
      end if;
    end if;
  end Match;

  procedure Second_Term ( f,g : in Poly; t : in Term; tol : in double_float;
                          output : in boolean;
                          c : out Complex_Number; w : out integer32 ) is

    k : constant integer32 := Standard_Power_Transformations.Pivot(t.dg.all);
    m : constant Matrix(1..2,1..2)
      := Standard_Power_Transformations.Eliminate(t.dg.all,k);
    tf : Poly := Standard_Initial_Forms.Transform(f,m);
    tg : Poly := Standard_Initial_Forms.Transform(g,m);
    sf : Poly := Shift(tf,k);
    sg : Poly := Shift(tg,k);
    sf1 : constant Standard_Complex_Vectors.Vector
        := Initial_Coefficients(sf,k);
    sg1 : constant Standard_Complex_Vectors.Vector
        := Initial_Coefficients(sg,k);
    df1 : constant Complex_Number := Derivative(sf1,t.cf);
    dg1 : constant Complex_Number := Derivative(sg1,t.cf);
    ind : constant integer32 := Other_Index(k);
    pf : Poly := Eval(sf,t.cf,ind);
    pg : Poly := Eval(sg,t.cf,ind);

  begin
    Initialize_Symbol_Table("x","y");
    if output then
      put("The transformation with pretropism ");
      put(t.dg.all); put_line(" :"); 
    end if;
   -- put_line(tf); put_line(tg);
   -- put_line("After the shift : "); put_line(sf); put_line(sg);
    if k = 1
     then Initialize_Symbol_Table("y");
     else Initialize_Symbol_Table("x");
    end if;
   -- put_line("After projection f :"); put_line(pf);
   -- put_line("After projection g :"); put_line(pg);
   -- put_line("Initial coefficients of f :"); put_line(sf1);
   -- put_line("Initial coefficients of g :"); put_line(sg1);
    Match(pf,pg,df1,dg1,tol,output,c,w);
    Clear(tf); Clear(tg); Clear(sf); Clear(sg); Clear(pf); Clear(pg);
  end Second_Term;

  procedure Second_Terms
              ( f,g : in Poly; t : in List_of_Terms;
                tol : in double_float; output : in boolean;
                s : out List_of_Germs; fail : out boolean ) is

    tmp : List_of_Terms := t;
    s_last : List_of_Germs;
    trm : Term;
    c : Complex_Number;
    w : integer32;

  begin
    fail := true;
    while not Is_Null(tmp) loop
      trm := Head_Of(tmp);
      Second_Term(f,g,trm,tol,output,c,w);
      if w /= 0 then
        if output
         then put("the coefficient : "); put(c); new_line;
        end if;
        declare
          gt : Germ;
        begin
          Copy(trm,gt.t); gt.c := c; gt.w := w;
          Append(s,s_last,gt);
        end;
        fail := false;
      else
        if output
         then put_line("no second term found");
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Second_Terms;

-- II. evaluate 

  function Evaluate ( t : Term; g : Germ ) return Poly is

    res : Poly := Null_Poly;
    first,second : Term;
    h : constant Term := g.t;
    k : constant integer32 := Standard_Power_Transformations.Pivot(h.dg.all);
    m : constant Matrix(1..2,1..2)
      := Standard_Power_Transformations.Eliminate(h.dg.all,k);
    a : constant integer32 := m(2,1)*t.dg(1);
    b : constant integer32 := m(2,2)*t.dg(2);

  begin
    first := Standard_Binomial_Factors.Evaluate(t,h);
    res := Create(first);
    second.dg := new Standard_Integer_Vectors.Vector(1..1);
    second.dg(1) := first.dg(1) + g.w;
    second.cf := Create(0.0);
    if k = 2 then
      if t.dg(1) /= 0 then
        second.cf := Create(double_float(t.dg(1)));
        if t.dg(1) /= 1
         then Mul(second.cf,h.cf**integer(t.dg(1)-1));
        end if;
      end if;
    else
      if a = 0 then
        if b = 1 then
          second.cf := Create(1.0);
        elsif b /= 0 then
          second.cf := h.cf**integer(b-1);
          Mul(second.cf,double_float(b));
        end if;
      elsif b = 0 then
        if a = 1 then
          second.cf := Create(1.0);
        else
          second.cf := h.cf**integer(a-1);
          Mul(second.cf,double_float(a));
        end if;
      else
        if a+b-1 = 0
         then second.cf := Create(1.0);
         else second.cf := h.cf**integer(a+b-1);
        end if;
        Mul(second.cf,double_float(a+b));
      end if;
    end if;
    if second.cf /= Create(0.0) then
      Mul(second.cf,t.cf);
      Mul(second.cf,g.c);
      Add(res,second);
    end if;
    Clear(first); Clear(second);
    return res;
  end Evaluate;

  function Evaluate ( p : Poly; g : Germ ) return Poly is
   
    res : Poly := Null_Poly;

    procedure Eval_Term ( t : Term; continue : out boolean ) is
 
      v : Poly := Evaluate(t,g);

    begin
      Add(res,v);
      Clear(v);
      continue := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Evaluate;

  function Order ( p : Poly; tol : double_float ) return integer32 is

    res : integer32 := integer32'first;
    first_time : boolean := true;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      if AbsVal(t.cf) > tol then
        if first_time then
          res := t.dg(t.dg'first);
          first_time := false;
        elsif t.dg(t.dg'first) < res then
          res := t.dg(t.dg'first);
        end if;
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Order;

end Standard_Puiseux_Certificates;
