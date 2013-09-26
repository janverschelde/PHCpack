with unchecked_deallocation;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;

--with text_io,integer_io; use text_io,integer_io;

package body Generic_Laur_Poly_Functions is

-- NOTE : Evaluation is not guaranteed to work for negative exponents...

-- DATA STRUCTURES : 

  MAX_INT : constant integer32 := 100000;

--  type kind is (coefficient,polynomial);
--  type Poly_Rec ( k : kind := coefficient ) is record
--    case k is
--      when coefficient => c : number;
--      when polynomial  => p : Eval_Poly;
--    end case;
--  end record;
--  type Coeff_Poly_Rec ( k : kind := coefficient ) is record
--    case k is
--      when coefficient => c : integer;
--      when polynomial  => p : Eval_Coeff_Poly;
--    end case;
--  end record;

--  type Eval_Poly_Rep is array(integer range <>) of Poly_Rec;
--  type Eval_Coeff_Poly_Rep is array(integer range <>) of Coeff_Poly_Rec;

  procedure free is new unchecked_deallocation(Eval_Poly_Rep,Eval_Poly);
  procedure free is
    new unchecked_deallocation(Eval_Coeff_Poly_Rep,Eval_Coeff_Poly);

-- AUXILIARY OPERATIONS :

  function Maximum ( a,b : integer32 ) return integer32 is
  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

  function Minimum ( a,b : integer32 ) return integer32 is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Minimum;

  function Convert ( c : number; n : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the corresponding value for c, when it lies in 1..n,
  --   otherwise 0 is returned.

  begin
    for i in 1..n loop
      if c = Create(integer(i))
       then return i;
      end if;
    end loop;
    return 0;
  end Convert;

  procedure Initialize ( evpr : in out Eval_Poly_Rep ) is
  begin
    for i in evpr'range loop
      declare
        nullpr : Poly_Rec(polynomial);
      begin
        nullpr.p := null;
        evpr(i) := nullpr;
      end;
    end loop;
  end Initialize;

  procedure Initialize ( evpr : in out Eval_Coeff_Poly_Rep ) is
  begin
    for i in evpr'range loop
      declare
        nullpr : Coeff_Poly_Rec(polynomial);
      begin
        nullpr.p := null;
        evpr(i) := nullpr;
      end;
    end loop;
  end Initialize;

  procedure Initialize ( p : in Poly; cff : out Vector; dgm : out Matrix ) is

  -- DESCRIPTION :
  --   Returns a vector/matrix representation of the polynomial p.
  --   Starts filling in backwards, since the highest degrees are in front
  --   of the list.  This could reduce the sorting time later.

    ind : integer32 := cff'last+1;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind-1;
      Copy(t.cf,cff(ind));
      for j in t.dg'range loop
        dgm(ind,j) := t.dg(j);
      end loop;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
  end Initialize;

  procedure Swap ( cff : in out Vector; dgm : in out Matrix;
                   i,j,k : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the ith row with the jth row, starting from the kth column.

    tmpcff : number;
    tmpdgm : integer32;

  begin
    Copy(cff(i),tmpcff);
    Copy(cff(j),cff(i));
    Copy(tmpcff,cff(j));
    Clear(tmpcff);
    for kk in k..dgm'last(2) loop
      tmpdgm := dgm(i,kk);
      dgm(i,kk) := dgm(j,kk);
      dgm(j,kk) := tmpdgm;
    end loop;
  end Swap;

  procedure Sort ( cff : in out Vector; dgm : in out Matrix;
                   i1,i2,k : in integer32 ) is

  -- DESCRIPTION :
  --   Sorts the elements in the kth column of the degree matrix, in
  --   the range i1..i2.  The coefficient vector gets sorted along.

    ind,min : integer32;

  begin
    for i in i1..i2 loop                  -- sort by swapping minimal element
      min := dgm(i,k);
      ind := i;
      for j in i+1..i2 loop               -- search for minimal element
        if dgm(j,k) < min then
          min := dgm(j,k);
          ind := j;
        end if;
      end loop;
      if ind /= i                         -- swap cff and deg
       then Swap(cff,dgm,i,ind,k);
      end if;
    end loop;
  end Sort;

  procedure Create ( cff : in out Vector; dgm : in out Matrix;
                     i1,i2,k : in integer32; ep : out Eval_Poly ) is

  -- DESCRIPTION :
  --   Returns in ep a nested Horner scheme to evaluate a polynomial given
  --   in vector/matrix representation with coefficients in cff and degrees
  --   in deg.  The range being considered is i1..i2, from the kth column.

  -- REQUIRED :
  --   The entries in the kth column of deg in i1..i2 are sorted
  --   in increasing order.

    min : constant integer32 := Minimum(0,dgm(i1,k));
    max : constant integer32 := Maximum(0,dgm(i2,k));
    evpr : Eval_Poly_Rep(min..max);
    ind,j1,j2 : integer32;

  begin
    Initialize(evpr);
    if k = dgm'last(2) then                   -- polynomial in one unknown
      for i in i1..i2 loop
        declare
          pr : Poly_Rec(Coefficient);
        begin
          Copy(cff(i),pr.c);
          evpr(dgm(i,k)) := pr;
        end;
      end loop;
    else
      ind := i1;                              -- recursive call
      while ind <= i2 loop
        j1 := ind; j2 := ind;
        while j2 < i2 and then dgm(j1,k) = dgm(j2+1,k) loop
          j2 := j2 + 1;
        end loop;
        declare
          pr : Poly_Rec(Polynomial);
        begin
          Sort(cff,dgm,j1,j2,k+1);
          Create(cff,dgm,j1,j2,k+1,pr.p);
          evpr(dgm(ind,k)) := pr;
        end;
        ind := j2+1;
      end loop;
    end if;
    ep := new Eval_Poly_Rep'(evpr);
  end Create;

-- CONSTRUCTORS :

  function Create ( p : Poly; n : natural; max,min : integer32 )
                  return Eval_Poly is

  -- DESCRIPTION :
  --   An evaluable polynomial is returned for p,
  --   with max = Maximum(0,Maximal_Degree(p,x1)),
  --   with min = Minimum(0,Minimal_Degree(p,x1)),
  --   and n = Number_of_Unknowns(p) >= 1.

    res : Eval_Poly;
    evpr : Eval_Poly_Rep(min..max);
    terms : array(min..max) of Poly := (min..max => Null_Poly);

    procedure Add_Term1 ( t : in Term; cont : out boolean ) is

      pr : Poly_Rec(coefficient);

    begin
      copy(t.cf,pr.c);
      evpr(t.dg(t.dg'first)) := pr;
      cont := true;
    end Add_Term1;
    procedure Add_Terms1 is new Visiting_Iterator(Add_Term1);

    procedure Add_Term ( t : in Term; cont : out boolean ) is

      nt : Term;

    begin
      nt.cf := t.cf;
      nt.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-1);
      for i in nt.dg'range loop
        nt.dg(i) := t.dg(i+1);
      end loop;
      Clear(nt);
      Add(terms(t.dg(t.dg'first)),nt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    Initialize(evpr);
    if n = 1 then
      Add_Terms1(p);
    else
      Add_Terms(p);
      for i in terms'range loop
        declare
          pr : Poly_Rec(polynomial);
        begin
          pr.p := Create(terms(i));
          evpr(i) := pr;
          Clear(terms(i));
        end;
      end loop;
    end if;
    res := new Eval_Poly_Rep'(evpr);
    return res;
  end Create;

  function Create ( p : Poly; n,nb : integer32; max,min : integer32 )
                  return Eval_Coeff_Poly is

  -- DESCRIPTION :
  --   An evaluable polynomial is returned for p, with
  --   max = Maximum(0,Maximal_Degree(p,x1)), in case deg(p,x1) < 0,
  --   min = Minimum(0,Minimal_Degree(p,x1)), to have 0 in range,
  --   n = Number_of_Unknowns(p) >= 1 and nb = Number_of_Terms(p).
  --   The coefficients of p are converted natural numbers.

    res : Eval_Coeff_Poly;
    evpr : Eval_Coeff_Poly_Rep(min..max);
    used : array(evpr'range) of boolean := (evpr'range => false);
    terms : array(min..max) of Poly := (min..max => Null_Poly);

    procedure Add_Term1 ( t : in Term; cont : out boolean ) is

      pr : Coeff_Poly_Rec(coefficient);

    begin
      pr.c := Convert(t.cf,nb);
      evpr(t.dg(t.dg'first)) := pr;
      used(t.dg(t.dg'first)) := true;
      cont := true;
    end Add_Term1;
    procedure Add_Terms1 is new Visiting_Iterator(Add_Term1);

    procedure Add_Term ( t : in Term; cont : out boolean ) is

      nt : Term;

    begin
      nt.cf := t.cf;
      nt.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-1);
      for i in nt.dg'range loop
        nt.dg(i) := t.dg(i+1);
      end loop;
      Add(terms(t.dg(t.dg'first)),nt);
      Clear(nt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    Initialize(evpr);
    if n = 1 then
      for i in evpr'range loop   -- initialization
        declare
          nullpr : Coeff_Poly_Rec(polynomial);
        begin
          nullpr.p := null;
          evpr(i) := nullpr;
        end;
      end loop;
      Add_Terms1(p);
    else
      Add_Terms(p);
      for i in terms'range loop
        declare
          pr : Coeff_Poly_Rec(polynomial);
          ind,max,min : integer32;
        begin
          if terms(i) = Null_Poly then
            pr.p := null;
          else
            ind := Head(terms(i)).dg'first;
            max := Maximum(0,Maximal_Degree(terms(i),ind));
            min := Minimum(0,Minimal_Degree(terms(i),ind));
            pr.p := Create(terms(i),n-1,nb,max,min);
          end if;
          evpr(i) := pr;
          used(i) := true;
          Clear(terms(i));
        end;
      end loop;
    end if;
    for i in used'range loop
      if not used(i) then
        declare
          pr : Coeff_Poly_Rec(polynomial);
        begin
          pr.p := null;
          evpr(i) := pr;
        end;
      end if;
    end loop;
    res := new Eval_Coeff_Poly_Rep'(evpr);
    return res;
  end Create;

  function Create ( p : Poly ) return Eval_Poly is

    res : Eval_Poly := null;
    nbvar : constant integer32 := integer32(Number_of_Unknowns(p));

  begin
    if (p /= Null_Poly) and then (nbvar /= 0) then
      declare
        nbtms : constant integer32 := integer32(Number_of_Terms(p));
        cff : Vector(1..nbtms);
        dgm : Matrix(1..nbtms,1..nbvar);
      begin
        Initialize(p,cff,dgm);
        Sort(cff,dgm,1,nbtms,1);
        Create(cff,dgm,1,nbtms,1,res);
        Clear(cff);
      end;
    end if;
    return res;
  end Create;

  function Create ( p : Poly ) return Eval_Coeff_Poly is

    res : Eval_Coeff_Poly;
    lp : Poly := Null_Poly;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    nb : constant integer32 := integer32(Number_of_Terms(p));
    max,min,ind : integer32;
    cnt : integer32 := 0;

    procedure Label_Term ( t : in Term; cont : out boolean ) is

      lt : Term;

    begin
      cnt := cnt + 1;
      lt.cf := Create(integer(cnt));
      lt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      Add(lp,lt);
      Clear(lt);
      cont := true;
    end Label_Term;
    procedure Label_Terms is new Visiting_Iterator(Label_Term);

  begin
    if (p = Null_Poly) or else (nb = 0) then
      return null;
    else
      Label_Terms(p);
      ind := Head(p).dg'first;
      max := Maximum(0,Maximal_Degree(p,ind));
      min := Minimum(0,Minimal_Degree(p,ind));
      res := Create(lp,n,nb,max,min);
      Clear(lp);
    end if;
    return res;
  end Create;

  procedure Diff ( p : in Poly; i : in integer32;
                   cp : out Eval_Coeff_Poly; m : out Vector ) is

    nb : constant integer32 := integer32(Number_of_Terms(p));
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    ind,cnt,max,min : integer32;
    dp : Poly := Null_Poly;

    procedure Diff_Term ( t : in Term; cont : out boolean ) is

      dt : Term;

    begin
      cnt := cnt + 1;
      if t.dg(i) /= 0 then
        dt.cf := Create(integer(cnt));
        dt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
        if t.dg(i) < 0  -- this is a patch ...
         then m(cnt) := -Create(integer(-t.dg(i)));
         else m(cnt) := Create(integer(t.dg(i)));
        end if;
        dt.dg(i) := dt.dg(i) - 1;
        Add(dp,dt);
        Clear(dt);
      else
        m(cnt) := Create(0);
      end if;
      cont := true;
   -- exception
   --   when others =>
   --     put("exception occurred in Diff_Term at cnt = "); put(cnt,1);
   --     new_line;
   --     put("  m'first = "); put(m'first,1); 
   --     put("  m'last = "); put(m'last,1);
   --     put("  t.dg(i) = "); put(t.dg(i),1); new_line;
   --     raise;
    end Diff_Term;
    procedure Diff_Terms is new Visiting_Iterator(Diff_Term);

  begin
    cnt := m'first-1;  -- changed this from cnt := 0;
    Diff_Terms(p);
    if dp /= Null_Poly then
      ind := Head(dp).dg'first;
      max := Maximum(0,Maximal_Degree(dp,ind));
      min := Minimum(0,Minimal_Degree(dp,ind));
      cp := Create(dp,n,nb,max,min);
    end if;
  end Diff;

  function Coeff ( p : Poly ) return Vector is

    res : Vector(1..integer32(Number_of_Terms(p)));
    cnt : integer32 := 0;

    procedure Collect_Term ( t : in Term; cont : out boolean ) is
    begin
      cnt := cnt + 1;
      copy(t.cf,res(cnt));
      cont := true;
    end Collect_Term;
    procedure Collect_Terms is new Visiting_Iterator(Collect_Term);

  begin
    Collect_Terms(p);
    return res;
  end Coeff;

-- EVALUATORS :

  function Eval ( p : Poly; x : number; i : integer32 ) return Poly is

    res : Poly := Null_Poly;

    procedure Eval_Term ( t : in Term; cont : out boolean ) is

      nt : Term;

    begin
      copy(t.cf,nt.cf);
      nt.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-1);
      for j in t.dg'range loop
        if j < i then
          nt.dg(j) := t.dg(j);
        elsif j > i then
          nt.dg(j-1) := t.dg(j);
        elsif t.dg(i) > 0 then
          for k in 1..t.dg(i) loop
            Mul(nt.cf,x);
          end loop;
        elsif t.dg(i) < 0 then
          for k in 1..(-t.dg(i)) loop
            Div(nt.cf,x);
          end loop;
        end if;
      end loop;
      Add(res,nt);
      Clear(nt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( d : Degrees; c : number; x : Vector ) return number is

    res : number;

  begin
    copy(c,res);
    for i in d'range loop
      for j in 1..(-d(i)) loop
        Div(res,x(i));
      end loop;
      for j in 1..d(i) loop
        Mul(res,x(i));
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( t : Term; x : Vector ) return number is
  begin
    return Eval(t.dg,t.cf,x);
  end Eval;

  function Eval ( t : Term; c : number; x : Vector ) return number is
  begin
    return Eval(t.dg,c,x);
  end Eval;

  function Eval ( p : Poly; x : Vector ) return number is

    res : number;

    procedure Eval_Term ( t : in Term; cont : out boolean ) is

      tmp : number := Eval(t,x);

    begin
      Add(res,tmp);
      Clear(tmp);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Copy(zero,res);
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : Poly; c,x : Vector ) return number is

    res : number;
    cnt : integer32 := 1;

    procedure Eval_Term ( t : in Term; cont : out boolean ) is

      tmp : number := Eval(t,c(cnt),x);

    begin
      Add(res,tmp);
      Clear(tmp);
      cnt := cnt + 1;
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Copy(zero,res);
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( vp : Eval_Poly_Rep; x : Vector;
                  i : integer32 ) return number;

  function Eval ( vprec : Poly_Rec; x : Vector;
                  i : integer32 ) return number is

    res : number;

  begin
    if vprec.k = coefficient then
      copy(vprec.c,res);
    elsif vprec.p = null then
      copy(zero,res);
    else
      res := Eval(vprec.p.all,x,i);
    end if;
    return res;
  end Eval;

  function Eval ( vp : Eval_Poly_Rep; x : Vector; i : integer32 )
                return number is

    res,val,tmp : number;

  begin
    if vp'first = 0 and vp'last = 0 then
      return Eval(vp(0),x,i+1);
    else
      if vp'last >= 0 then
        res := Eval(vp(vp'last),x,i+1);
        for j in reverse 0..(vp'last-1) loop
          Mul(res,x(i));
          if (vp(j).k = coefficient) or else (vp(j).p /= null) then
            tmp := Eval(vp(j),x,i+1);
            Add(res,tmp); Clear(tmp);
          end if;
        end loop;
      else
        copy(zero,res);
      end if;
      if vp'first < 0 then
        val := Eval(vp(vp'first),x,i+1);
        for j in (vp'first+1)..(-1) loop
          Div(val,x(i));
          if (vp(j).k = coefficient) or else (vp(j).p /= null) then
            tmp := Eval(vp(j),x,i+1);
            Add(val,tmp); Clear(tmp);
          end if;
        end loop;
        Div(val,x(i));
        Add(res,val); Clear(val);
      end if;
      return res;
    end if;
  end Eval;

  function Eval ( p : Eval_Poly; x : Vector ) return number is
  begin
    if p = null then
      declare
        res : number;
      begin
        Copy(zero,res);
        return res;
      end;
    else
      return Eval(p.all,x,x'first);
    end if;
  end Eval;

  function Eval ( vp : Eval_Coeff_Poly_Rep; c,x : Vector; i : integer32 )
                return number;

  function Eval ( vprec : Coeff_Poly_Rec; c,x : Vector; i : integer32 )
                return number is

    res : number;

  begin
    if vprec.k = coefficient then
      copy(c(vprec.c),res);
    elsif vprec.p = null then
      copy(zero,res);
    else
      res := Eval(vprec.p.all,c,x,i);
    end if;
    return res;
  end Eval;

  function Eval ( vp : Eval_Coeff_Poly_Rep; c,x : Vector; i : integer32 )
                return number is

    res,val,tmp : number;

  begin
    if vp'first = 0 and vp'last = 0 then
      return Eval(vp(0),c,x,i+1);
    else
      if vp'last >= 0 then
        res := Eval(vp(vp'last),c,x,i+1);
        for j in reverse 0..(vp'last-1) loop
          Mul(res,x(i));
          if (vp(j).k = coefficient) or else (vp(j).p /= null) then
            tmp := Eval(vp(j),c,x,i+1);
            Add(res,tmp); Clear(tmp);
          end if;
        end loop;
      else
        copy(zero,res);
      end if;
      if vp'first < 0 then
        val := Eval(vp(vp'first),c,x,i+1);
        for j in (vp'first+1)..(-1) loop
          Div(val,x(i));
          if (vp(j).k = coefficient) or else (vp(j).p /= null) then
            tmp := Eval(vp(j),c,x,i+1);
            Add(val,tmp); Clear(tmp);
          end if;
        end loop;
        Div(val,x(i));
        Add(res,val); Clear(val);
      end if;
      return res;
    end if;
  end Eval;

  function Eval ( p : Eval_Coeff_Poly; c,x : Vector ) return number is
  begin
    if p = null then
      declare
        res : number;
      begin
        Copy(zero,res);
        return res;
      end;
    else
      return Eval(p.all,c,x,x'first);
    end if;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( p : in out Eval_Poly ) is
  begin
    if p /= null then
      declare
        vp : Eval_Poly_Rep renames p.all;
      begin
        for i in vp'range loop
          if vp(i).k = coefficient
           then Clear(vp(i).c);
           else Clear(vp(i).p);
          end if;
        end loop;
      end;
      free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Eval_Coeff_Poly ) is
  begin
    if p /= null then
      declare
        vp : Eval_Coeff_Poly_Rep renames p.all;
      begin
        for i in vp'range loop
          if vp(i).k /= coefficient
           then Clear(vp(i).p);
          end if;
        end loop;
      end;
      free(p);
    end if;
  end Clear;

end Generic_Laur_Poly_Functions;
