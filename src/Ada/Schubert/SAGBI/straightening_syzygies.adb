with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;

package body Straightening_Syzygies is

-- AUXILIARIES :

  function Create ( v1,v2 : Vector ) return Bracket_Term is

    b1 : Bracket(v1'range);
    b2 : Bracket(v2'range);
    sig1,sig2 : integer32;
    bm : Bracket_Monomial;
    bt : Bracket_Term;

  begin
    Create(v1,b1,sig1);
    if Is_Zero(b1) then
      bt.coeff := Create(0.0);
    else
      Create(v2,b2,sig2);
      if Is_Zero(b2) then
        bt.coeff := Create(0.0);
      else
        bm := Create(b1);
        Multiply(bm,b2);
        if sig1*sig2 > 0
         then bt.coeff := Create(1.0);
         else bt.coeff := -Create(1.0);
        end if;
      end if;
    end if;
    bt.monom := bm;
    return bt;
  end Create;

  function Sort ( v : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the sorted vector v, in increasing order.

    res : Vector(v'range) := v;
    ind : integer32;
    min : natural32;

  begin
    for i in v'first..v'last-1 loop
      min := res(i);
      ind := i;
      for j in i+1..res'last loop
        if res(j) < min
         then min := res(j); ind := j;
        end if;
      end loop;
      if ind /= i
       then res(ind) := res(i); res(i) := min;
      end if;
    end loop;
    return res;
  end Sort;

  function Sign ( v1,v2 : Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns the sign of the permutation (v1,v2).

    d : constant integer32 := v1'length+v2'length;
    b : Bracket(1..d);
    v : Vector(1..d);
    sig : integer32;

  begin
    v(v1'range) := v1;
    v(v1'last+1..v'last) := v2;
    Create(v,b,sig);
    return sig;
  end Sign;

  function Complement ( n : integer32; v : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the complement of the vector w.r.t. the set 1..n.

    res : Vector(1..n-v'length);
    cnt : integer32 := 0;
    found : boolean;

  begin
    for i in 1..n loop
      found := false;
      for j in v'range loop
        if integer32(v(j)) = i
         then found := true; exit;
        end if;
      end loop;
      if not found 
       then cnt := cnt + 1; res(cnt) := natural32(i);
      end if;
    end loop;
    return res;
  end Complement;

  procedure Enumerate
              ( start,i,n : in integer32; accu : in out Vector;
                s : in integer32; b1,b2 : in Bracket;
                frame : in Vector; bp : in out Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Enumerates all subsets of 1..n, of size accu'length, starting to
  --   fill up accu(i) with entries in start..n.

    v1,v2 : Vector(1..b1'last);

  begin
    if i > accu'last then
      declare
        comp : constant Vector := Complement(n,accu);
      begin
        v1(b1'first..s-1)
          := Standard_Natural_Vectors.Vector(b1(b1'first..s-1));
        for j in accu'range loop
          v1(s+j-1) := frame(integer32(accu(j)));
        end loop;
        v2(s+1..b2'last) := Standard_Natural_Vectors.Vector(b2(s+1..b2'last));
        for j in comp'range loop
          v2(j) := frame(integer32(comp(j)));
        end loop;
        declare
          bt : Bracket_Term := Create(v1,v2);
        begin
          if bt.coeff /= Create(0.0) then
            if Sign(accu,comp) > 0
             then Frontal_Add(bp,bt);
             else Frontal_Min(bp,bt);
            end if;
          end if;
          Clear(bt);
        end;
      end;
    else
      for ell in start..n loop
        accu(i) := natural32(ell);
        Enumerate(ell+1,i+1,n,accu,s,b1,b2,frame,bp);
      end loop;
    end if;
  end Enumerate;

  function Create ( b1,b2 : Bracket ) return Bracket_Term is

    bm : Bracket_Monomial := Create(b1);
    bt : Bracket_Term;

  begin
    Multiply(bm,b2);
    bt.coeff := Create(1.0);
    bt.monom := bm;
    return bt;   
  end Create;

  procedure Enumerate2
              ( n,d,k,start : in integer32; b : in Bracket;
                accu : in out Bracket;
                cntstd,cntnonstd : in out integer32;
                nonstd : in out Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets, with b < accu and with
  --   a test whether the pair b*accu forms a Standard monomial.

  -- ON ENTRY :
  --   n             total number of indices to choose from;
  --   d             number of indices in the brackets;
  --   k             current entry in accu that is to be filled;
  --   start         indices will be taken in between start..n;
  --   b             previously enumerated bracket;
  --   accu          accumulating parameter, filled in upto (k-1)th entry;
  --   cntstd        current number of quadratic standard monomials;
  --   cntnonstd     current number of quadratic nonstandard monomials;
  --   nonstd        current polynomial of quadratic nonstandard monomials.

  -- ON RETURN :
  --   accu          updated accumulating parameter, accu(k) is filled in;
  --   cnstd         updated number of quadratic standard monomials;
  --   cntnonstd     updated number of quadratic nonstandard monomials;
  --   nonstd        updated polynomial of quadratic nonstandard monomials.

    s : natural32;

  begin
    if k > d then
      if b < accu then
        s := Brackets.Is_Standard(b,accu);
        if s = 0
         then cntstd := cntstd + 1;
         else cntnonstd := cntnonstd + 1; Add(nonstd,Create(b,accu));
        end if;
      end if;
    else
      for ell in start..n loop
        accu(k) := natural32(ell);
        Enumerate2(n,d,k+1,ell+1,b,accu,cntstd,cntnonstd,nonstd);
      end loop;
    end if;
  end Enumerate2;

  procedure Enumerate1
              ( n,d,k,start : integer32; acc1,acc2 : in out Bracket;
                cntstd,cntnonstd : in out integer32;
                nonstd : in out Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Lexicographic enumeration of all brackets, with acc1 < acc2 and with
  --   a test whether the pair acc1*acc2 forms a Standard monomial.
  --   Counts the standard and nonstandard monomials and adds every
  --   nonStandard monomial to the polynomial nonstd.

  begin
    if k > d then
      Enumerate2(n,d,1,integer32(acc1(1)),acc1,acc2,cntstd,cntnonstd,nonstd);
    else
      for ell in start..n loop
        acc1(k) := natural32(ell);
        Enumerate1(n,d,k+1,ell+1,acc1,acc2,cntstd,cntnonstd,nonstd);
      end loop;
    end if;
  end Enumerate1;

  procedure Enumerate3
              ( n,d,k,start : in integer32; accu : in out Vector;
                bp : in out Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Lexicographic enumerations of all vectors with d-entries, out
  --   of a set of n natural numbers.

  -- ON ENTRY :
  --   n             total number of indices to choose from;
  --   d             number of indices in the brackets;
  --   k             current entry in accu that is to be filled;
  --   start         indices will be taken in between start..n;
  --   accu          accumulating parameter, filled in upto (k-1)th entry;

  -- ON RETURN :
  --   accu          filled in upto the kth entry.

  begin
    if k > d then
      declare
        comp : constant Vector := Complement(n,accu);
        acc0 : Vector(1..d+1);
        t : Bracket_Term;
      begin
       -- put(" accu : "); put(accu); put(" complement : "); put(comp);
        acc0(1) := 0;
        acc0(2..d+1) := accu(1..d);
       --  put(" acc0 : "); put(acc0); new_line;
        t := Create(acc0,comp);
        if Sign(accu,comp) > 0
         then Frontal_Add(bp,t);
         else Frontal_Min(bp,t);
        end if;
        Clear(t);
      end;
    else
      for ell in start..n-d+k loop
        accu(k) := natural32(ell);
        Enumerate3(n,d,k+1,ell+1,accu,bp);
      end loop;
    end if;
  end Enumerate3;

-- TARGET OPERATIONS :

  function Laplace_Expansion ( n,d : natural32 ) return Bracket_Polynomial is

    res : Bracket_Polynomial;
    acc : Vector(1..integer32(d));

  begin
    Enumerate3(integer32(n),integer32(d),1,1,acc,res);
    return res;
  end Laplace_Expansion;

  function Straightening_Syzygy ( b1,b2 : Bracket )
                                return Bracket_Polynomial is

    s : constant integer32 := integer32(Is_Standard(b1,b2));
    bm : Bracket_Monomial;
    bp : Bracket_Polynomial;

  begin
    if s = 0 then
      bm := Create(b1);
      Multiply(bm,b2);
      bp := Create(bm);
    else
      declare
        d1 : constant integer32 := b1'last+1;
        frame : Vector(1..d1);
        accu : Vector(1..d1-s);
      begin
        for i in s..b1'last loop
          frame(i-s+1) := b1(i);
        end loop;
        for i in b2'first..s loop
          frame(d1-s+i) := b2(i);
        end loop;
        frame := Sort(frame);
        Enumerate(1,1,d1,accu,s,b1,b2,frame,bp);
      end;
    end if;
    return bp;
  end Straightening_Syzygy;

  function Straightening_Syzygy ( b : Bracket_Monomial )
                                return Bracket_Polynomial is

    b1,b2 : Link_to_Bracket;
    bp : Bracket_Polynomial;

    procedure Get_Bracket ( bb : in Bracket; continue : out boolean ) is
    begin
      if b1 = null
       then b1 := new Bracket'(bb);
       else b2 := new Bracket'(bb);
      end if;
      continue := true;
    end Get_Bracket;
    procedure Get_Brackets is new Enumerate_Brackets(Get_Bracket);

  begin
    Get_Brackets(b);
    bp := Straightening_Syzygy(b1.all,b2.all);
    Clear(b1); Clear(b2);
    return bp;
  end Straightening_Syzygy;

  function nonStandard_Monomials
             ( n,d : natural32 ) return Bracket_Polynomial is

    nonstd : Bracket_Polynomial;
    b1,b2 : Bracket(1..integer32(d));
    cntstd,cntnonstd : integer32 := 0;

  begin
    Enumerate1(integer32(n),integer32(d),1,1,b1,b2,cntstd,cntnonstd,nonstd);
    return nonstd;
  end nonStandard_Monomials;

  procedure Enumerate_Syzygies ( p : in Bracket_Polynomial ) is

    procedure Process_Syzygy ( t : in Bracket_Term; continue : out boolean ) is
    begin
      Process(Straightening_Syzygy(t.monom),continue);
    end Process_Syzygy;
    procedure Process_Syzygies is new Enumerate_Terms(Process_Syzygy);

  begin
    Process_Syzygies(p);
  end Enumerate_Syzygies;

  function Straighten ( b1,b2 : Bracket ) return Bracket_Polynomial is

    bp : constant Bracket_Polynomial := Straightening_Syzygy(b1,b2);

  begin
    return bp;
  end Straighten;

  function Straighten ( b : Bracket_Monomial ) return Bracket_Polynomial is

    bp : Bracket_Polynomial;

  begin
    return bp;
  end Straighten;

  function Straighten ( b : Bracket_Term ) return Bracket_Polynomial is

    bp : Bracket_Polynomial;

  begin
    return bp;
  end Straighten;

  function Straighten ( b : Bracket_Polynomial ) return Bracket_Polynomial is

    bp : Bracket_Polynomial;

  begin
    return bp;
  end Straighten;

end Straightening_Syzygies;
