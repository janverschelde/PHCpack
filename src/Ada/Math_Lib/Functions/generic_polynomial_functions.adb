with unchecked_deallocation;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;

--with text_io,integer_io;                 use text_io,integer_io;
--with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;

package body Generic_Polynomial_Functions is

-- DATA STRUCTURES : 

--  type kind is (coefficient,polynomial);
--  type Poly_Rec is record
--    k : kind;
--    c : number;
--    p : Eval_Poly;
--  end record;
--  type Coeff_Poly_Rec is record
--    k : kind;
--    c : integer;
--    p : Eval_Coeff_Poly;
--  end record;

-- Note: the original variant records (with k as case) made "valgrind"
--   to report an error in the initialization to coefficient,
--   as c and p occupy the same memory location.

--  type Eval_Poly_Rep is array(integer range <>) of Poly_Rec;
--  type Eval_Coeff_Poly_Rep is array(integer range <>) of Coeff_Poly_Rec;

  procedure free is new unchecked_deallocation(Eval_Poly_Rep,Eval_Poly);
  procedure free is
    new unchecked_deallocation(Eval_Coeff_Poly_Rep,Eval_Coeff_Poly);

-- AUXILIARY OPERATIONS :

  function Convert ( c : number; n : natural32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the corresponding value for c, when it lies in 1..n,
  --   otherwise 0 is returned.

    diff : number;
    eqzero : boolean;

  begin
    for i in 1..n loop
     -- if c = Create(integer(i)) -- not okay for multiprecision
      diff := c - Create(integer(i));
      eqzero := Equal(diff,zero);
      if eqzero
       then return integer32(i);
      end if;
    end loop;
    return integer32(0);
  end Convert;

  procedure Initialize ( p : in Poly; cff : out Vector; deg : out Matrix ) is

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
        deg(ind,j) := t.dg(j);
      end loop;
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
   -- put_line("The degree matrix : "); put(deg);
  end Initialize;

  procedure Swap ( cff : in out Vector; deg : in out Matrix;
                   i,j,k : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the ith row with the jth row, starting from the kth column.

    tmpcff : number;
    tmpdeg : natural32;

  begin
    Copy(cff(i),tmpcff); -- tmpcff := cff(i);
    Copy(cff(j),cff(i)); -- cff(i) := cff(j);
    Copy(tmpcff,cff(j)); -- cff(j) := tmpcff;
    Clear(tmpcff);
    for kk in k..deg'last(2) loop
      tmpdeg := deg(i,kk);
      deg(i,kk) := deg(j,kk);
      deg(j,kk) := tmpdeg;
    end loop;
  end Swap;

  procedure Sort ( cff : in out Vector; deg : in out Matrix;
                   i1,i2,k : in integer32 ) is

  -- DESCRIPTION :
  --   Sorts the elements in the kth column of the degree matrix, in
  --   the range i1..i2.  The coefficient vector gets sorted along.

    min : natural32;
    ind : integer32;

  begin
    for i in i1..i2 loop                  -- sort by swapping minimal element
      min := deg(i,k);
      ind := i;
      for j in i+1..i2 loop               -- search for minimal element
        if deg(j,k) < min
         then min := deg(j,k);
              ind := j;
        end if;
      end loop;
      if ind /= i                         -- swap cff and deg
       then Swap(cff,deg,i,ind,k);
      end if;
    end loop;
   -- put_line("The sorted degree matrix : "); put(deg);
  end Sort;

 -- procedure Write_Structure ( p : in Eval_Poly; ind : in natural ) is
 -- begin
 --   for i in p'range loop
 --     for j in 1..ind loop
 --       put("  ");
 --     end loop;
 --     put("ep("); put(i,1); put(") is ");
 --     case p(i).k is
 --       when coefficient => put_line("coefficient");
 --       when polynomial  => if p(i).p = null
 --                            then put_line("empty");
 --                            else put_line("polynomial");
 --                                 Write_Structure(p(i).p,ind+1);
 --                           end if;
 --     end case;
 --   end loop;
 -- end Write_Structure;

  function Terminal_Degree
             ( deg : Matrix; i,k : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if deg(i,l) = 0 for all l in range k..deg'last(2).

  begin
    for j in k..deg'last(2) loop
      if deg(i,j) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Terminal_Degree;

  procedure Create ( cff : in out Vector; deg : in out Matrix;
                     i1,i2,k : in integer32; ep : out Eval_Poly ) is

  -- DESCRIPTION :
  --   Returns in ep a nested Horner scheme to evaluate a polynomial given
  --   in vector/matrix representation with coefficients in cff and degrees
  --   in deg.  The range being considered is i1..i2, from the kth column.

  -- REQUIRED :
  --   The entries in the kth column of deg in i1..i2 are sorted
  --   in increasing order.

    evpr : Eval_Poly_Rep(0..integer32(deg(i2,k)));
    ind,j1,j2 : integer32;

  begin
    for i in evpr'range loop
      evpr(i).k := Polynomial;
      evpr(i).p := null;
    end loop;
    if ((i2 = i1) and then Terminal_Degree(deg,i2,k)
                  and then Terminal_Degree(deg,i1,k)) then  -- constant term
      evpr(integer32(deg(i2,k))).k := Coefficient;
      Copy(cff(i2),evpr(integer32(deg(i2,k))).c);
    elsif k = deg'last(2) then                  -- polynomial in one unknown
      for i in i1..i2 loop
        evpr(integer32(deg(i,k))).k := Coefficient;
        Copy(cff(i),evpr(integer32(deg(i,k))).c);
      end loop;
    else
      ind := i1;                                           -- recursive call
      while ind <= i2 loop
        j1 := ind; j2 := ind;
        while j2 < i2 and then deg(j1,k) = deg(j2+1,k) loop
          j2 := j2 + 1;
        end loop;
        evpr(integer32(deg(ind,k))).k := Polynomial;
        Sort(cff,deg,j1,j2,k+1);
        Create(cff,deg,j1,j2,k+1,evpr(integer32(deg(ind,k))).p);
        ind := j2+1;
      end loop;
    end if;
    ep := new Eval_Poly_Rep'(evpr);
  --exception
  --  when others => put_line("exception in auxiliary creator"); raise;
  end Create;

-- CONSTRUCTORS :

  function Create ( p : Poly; n,nb : natural32; d : integer32 )
                  return Eval_Coeff_Poly is

  -- DESCRIPTION :
  --   An evaluable polynomial is returned for p, with d = Degree(p,x1),
  --   n = Number_of_Unknowns(p) >= 1 and nb = Number_of_Terms(p).
  --   The coefficients of p are converted natural numbers.

    res : Eval_Coeff_Poly;
    evpr : Eval_Coeff_Poly_Rep(0..d);
    terms : array(0..d) of Poly := (0..d => Null_Poly);

    procedure Add_Term1 ( t : in Term; cont : out boolean ) is

      pr : Coeff_Poly_Rec;

    begin
      pr.k := Coefficient;
      pr.c := Convert(t.cf,nb);
      evpr(integer32(t.dg(t.dg'first))) := pr;
      cont := true;
    end Add_Term1;
    procedure Add_Terms1 is new Visiting_Iterator(Add_Term1);

    procedure Add_Term ( t : in Term; cont : out boolean ) is

      nt : Term;

    begin
      Copy(t.cf,nt.cf);
      nt.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
      for i in nt.dg'range loop
        nt.dg(i) := t.dg(i+1);
      end loop;
      Add(terms(integer32(t.dg(t.dg'first))),nt);
      Clear(nt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    for i in evpr'range loop
      evpr(i).k := Polynomial;
      evpr(i).p := null;
    end loop;
    if n = 1 then
      Add_Terms1(p);
    else
      Add_Terms(p);
      for i in terms'range loop
        declare
          ind : integer32;
        begin
          evpr(i).k := Polynomial;
          if terms(i) = Null_Poly then
            evpr(i).p := null;
          else
            ind := Head(terms(i)).dg'first;
            evpr(i).p := Create(terms(i),n-1,nb,Degree(terms(i),ind));
            Clear(terms(i));
          end if;
        end;
      end loop;
    end if;
    res := new Eval_Coeff_Poly_Rep'(evpr);
    return res;
  --exception 
  --  when others => put_line("exception in second create"); raise;
  end Create;

  function Create ( p : Poly ) return Eval_Poly is

    res : Eval_Poly := null;
    nbvar : constant integer32 := integer32(Number_of_Unknowns(p));
    nbtms : constant integer32 := integer32(Number_of_Terms(p));
    cff : Vector(1..nbtms);
    deg : Matrix(1..nbtms,1..nbvar);

  begin
    if (p /= Null_Poly) and then (nbvar /= 0) then
      Initialize(p,cff,deg);
      Sort(cff,deg,1,nbtms,1);
      Create(cff,deg,1,nbtms,1,res);
      Clear(cff);
     -- put_line("The structure of the eval poly : ");
     -- Write_Structure(res,1);
    end if;
    return res;
  --exception
  --  when others => put_line("exception in main create"); raise;
  end Create;

  function Create ( p : Poly ) return Eval_Coeff_Poly is

    res : Eval_Coeff_Poly;
    lp : Poly := Null_Poly;
    n : constant natural32 := Number_of_Unknowns(p);
    nb : constant natural32 := Number_of_Terms(p);
    cnt : natural32 := 0;
    ind : integer32;

    procedure Label_Term ( t : in Term; cont : out boolean ) is

      lt : Term;

    begin
      cnt := cnt + 1;
      lt.cf := Create(integer(cnt));
      lt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
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
      res := Create(lp,n,nb,Degree(p,ind));
      Clear(lp);
    end if;
    return res;
  end Create;

  procedure Diff ( p : in Poly; i : in integer32;
                   cp : out Eval_Coeff_Poly; m : out Vector ) is

    nb : constant natural32 := Number_of_Terms(p);
    n : constant natural32 := Number_of_Unknowns(p);
    ind,cnt : integer32;
    dp : Poly := Null_Poly;

    procedure Diff_Term ( t : in Term; cont : out boolean ) is

      dt : Term;

    begin
      cnt := cnt + 1;
      if t.dg(i) > 0 then
        dt.cf := Create(integer(cnt));
        dt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
        m(cnt) := Create(integer(t.dg(i)));
        dt.dg(i) := dt.dg(i) - 1;
        Add(dp,dt);
        Clear(dt);
      else
        m(cnt) := Create(0);
      end if;
      cont := true;
    end Diff_Term;
    procedure Diff_Terms is new Visiting_Iterator(Diff_Term);

  begin
    cnt := 0;
    Diff_Terms(p);
    if dp /= Null_Poly then
      ind := Head(dp).dg'first;
      cp := Create(dp,n,nb,Degree(dp,ind));
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
      Copy(t.cf,nt.cf);
      nt.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
      for j in t.dg'range loop
        if j < i then
          nt.dg(j) := t.dg(j);
        elsif j > i then
          nt.dg(j-1) := t.dg(j);
        else
          for k in 1..t.dg(i) loop
            Mul(nt.cf,x);
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
    Copy(c,res);
    for i in d'range loop
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
    cnt : integer32 := 0;

    procedure Eval_Term ( t : in Term; cont : out boolean ) is

      tmp : number := Eval(t,c(cnt),x);

    begin
      cnt := cnt + 1;
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

  function Eval ( vp : Eval_Poly_Rep; x : Vector;
                  i : integer32 ) return number;

  function Eval ( vprec : Poly_Rec; x : Vector;
                  i : integer32 ) return number is

    res : number;

  begin
    if vprec.k = coefficient then
      Copy(vprec.c,res);
    elsif vprec.p = null then
      Copy(zero,res);
    else
      res := Eval(vprec.p.all,x,i);
    end if;
    return res;
  end Eval;

  function Eval ( vp : Eval_Poly_Rep; x : Vector; i : integer32 )
                return number is

    deg : constant integer32 := vp'length-1;
    res : number;

  begin
    if deg = 0 then
      res := Eval(vp(0),x,i+1);
    else
      res := Eval(vp(deg),x,i+1);
      for j in reverse 0..(deg-1) loop
        Mul(res,x(i));
        if (vp(j).k = coefficient) or else (vp(j).p /= null) then
          declare
            temp : number := Eval(vp(j),x,i+1);
          begin
            Add(res,temp);
            Clear(temp);
          end;
        end if;
      end loop;
    end if;
    return res;
  --exception
  --  when others => put("exception in eval poly rep"); raise;
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
      Copy(c(vprec.c),res);
    elsif vprec.p = null then
      Copy(zero,res);
    else
      res := Eval(vprec.p.all,c,x,i);
    end if;
    return res;
  end Eval;

  function Eval ( vp : Eval_Coeff_Poly_Rep; c,x : Vector;
                  i : integer32 ) return number is

    deg : constant integer32 := vp'length-1;
    res : number;

  begin
    if deg = 0 then
      res := Eval(vp(0),c,x,i+1);
    else
      res := Eval(vp(deg),c,x,i+1);
      for j in reverse 0..(deg-1) loop
        Mul(res,x(i));
        if (vp(j).k = coefficient) or else (vp(j).p /= null) then
          declare
            temp : number := Eval(vp(j),c,x,i+1);
          begin
            Add(res,temp);
            Clear(temp);
          end;
        end if;
      end loop;
    end if;
    return res;
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
      for i in p'range loop
        if p(i).k = coefficient
         then Clear(p(i).c);
         else Clear(p(i).p);
        end if;
      end loop;
      free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Eval_Coeff_Poly ) is
  begin
    if p /= null then
      for i in p'range loop
        if p(i).k /= coefficient
         then Clear(p(i).p);
        end if;
      end loop;
    free(p);
    end if;
  end Clear;

end Generic_Polynomial_Functions;
