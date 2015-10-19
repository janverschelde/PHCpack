package body Generic_Bracket_Polynomials is

-- CONSTRUCTORS :

  function Create ( m : Bracket_Monomial ) return Bracket_Polynomial is

    t : Bracket_Term;

  begin
    t.coeff := create(integer(1));
    t.monom := m;
    return Create(t);
  end Create;

  function Create ( t : Bracket_Term ) return Bracket_Polynomial is

    p : Bracket_Polynomial;

  begin
    Construct(t,p);
    return p;
  end Create;

  procedure Copy_Multiply ( t1 : in Bracket_Term;
                            t2 : in out Bracket_Term ) is
  begin
    t2.coeff := t1.coeff;
    Copy_Multiply(t1.monom,t2.monom);
  end Copy_Multiply;

  procedure Copy_Append ( t1 : in Bracket_Term; t2 : in out Bracket_Term ) is
  begin
    t2.coeff := t1.coeff;
    Copy_Append(t1.monom,t2.monom);
  end Copy_Append;

  procedure Copy ( p : in Bracket_Polynomial;
                   q : in out Bracket_Polynomial ) is

    tmp : constant Bracket_Polynomial := p;

  begin
    Clear(q);
    while not Is_Null(tmp) loop
      Add(q,Head_Of(tmp));
    end loop;
  end Copy;

-- SELECTORS :

  function Dimension ( t : Bracket_Term ) return natural32 is
  begin
    return Dimension(t.monom);
  end Dimension;

  function Dimension ( p : Bracket_Polynomial ) return natural32 is
  begin
    if Is_Null(p)
     then return 0;
     else return Dimension(Head_Of(p));
    end if;
  end Dimension;

-- COMPARISON OPERATIONS :

  function Is_Equal ( t1,t2 : Bracket_Term ) return boolean is
  begin
    return (t1.coeff = t2.coeff and then Is_Equal(t1.monom,t2.monom));
  end Is_Equal;

  function Is_Equal ( p,q : Bracket_Polynomial ) return boolean is

    tmp1 : Bracket_Polynomial := p;
    tmp2 : Bracket_Polynomial := q;

  begin
    while not Is_Null(tmp1) and not Is_Null(tmp2) loop
      if not Is_Equal(Head_Of(tmp1),Head_Of(tmp2))
       then return false;
       else tmp1 := Tail_Of(tmp1); tmp2 := Tail_Of(tmp2);
      end if;
    end loop;
    if Is_Null(tmp1) and Is_Null(tmp2)
     then return true;
     else return false;
    end if;
  end Is_Equal;

  function "<" ( t1,t2 : Bracket_Term ) return boolean is
  begin
    return t1.monom < t2.monom;
  end "<";

  function ">" ( t1,t2 : Bracket_Term ) return boolean is
  begin
    return t1.monom > t2.monom;
  end ">";

-- ARITHMETIC OPERATIONS :

  function "+" ( t : Bracket_Term; p : Bracket_Polynomial )
               return Bracket_Polynomial is

    res : Bracket_Polynomial;

  begin
    Copy(p,res);
    Add(res,t);
    return res;
  end "+";

  function "+" ( p : Bracket_Polynomial; t : Bracket_Term )
               return Bracket_Polynomial is

    res : Bracket_Polynomial;

  begin
    Copy(p,res);
    Add(res,t);
    return res;
  end "+";

  procedure Add ( p : in out Bracket_Polynomial; t : in Bracket_Term ) is

    tt : Bracket_Term;

  begin
    Copy_Multiply(t,tt);
    if Is_Null(p) then
      p := Create(tt);
    else
      declare
        first,second : Bracket_Polynomial;
        t1 : Bracket_Term;
      begin
        first := p; second := Tail_Of(p);
        t1 := Head_Of(first);
        if t > t1 then
          Construct(tt,p);
        elsif Is_Equal(t.monom,t1.monom) then
          t1.coeff := t1.coeff + t.coeff;
          if t1.coeff = Create(integer(0)) then
            Clear(t1);
            p := Tail_Of(p);
          else
            Set_Head(p,t1);
          end if;
        else
          while not Is_Null(second) loop     -- merge term in list
            t1 := Head_Of(second);
            if t > t1 then
              Construct(tt,second);
              Swap_Tail(first,second); exit;
            elsif Is_Equal(t.monom,t1.monom) then
              t1.coeff := t1.coeff + t.coeff;
              if t1.coeff = Create(integer(0)) then
                Clear(t1);
                Swap_Tail(first,second);
              else
                Set_Head(second,t1);
              end if;
              exit;
            end if;
            first := Tail_Of(first);
            second := Tail_Of(second);
          end loop;
          if Is_Null(second)                -- then first points to last
           then Append(p,first,tt);         --   element of the list p
          end if;
        end if;
      end;
    end if;
  end Add;

  procedure Frontal_Add ( p : in out Bracket_Polynomial;
                          t : in Bracket_Term ) is

    tt : Bracket_Term;

  begin
    Copy_Multiply(t,tt);   
    Construct(tt,p);
  end Frontal_Add;

  procedure Frontal_Construct ( p : in out Bracket_Polynomial;
                                t : in Bracket_Term ) is

    tt : Bracket_Term;

  begin
    Copy_Append(t,tt);   
    Construct(tt,p);
  end Frontal_Construct;

  procedure Frontal_Min ( p : in out Bracket_Polynomial;
                          t : in Bracket_Term ) is

    mt : constant Bracket_Term := -t;

  begin
    Construct(mt,p);
  end Frontal_Min;

  function "+" ( p,q : Bracket_Polynomial ) return Bracket_Polynomial is

    res : Bracket_Polynomial;

  begin
    Copy(p,res);
    Add(res,q);
    return res;
  end "+";

  procedure Add ( p : in out Bracket_Polynomial;
                  q : in Bracket_Polynomial ) is

    tmp : Bracket_Polynomial := q;

  begin
    while not Is_Null(tmp) loop
      Add(p,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end Add;

  function "-" ( t : Bracket_Term ) return Bracket_Term is

    res : Bracket_Term;

  begin
    Copy_Multiply(t.monom,res.monom);
    res.coeff := -t.coeff;
    return res;
  end "-";

  procedure Min ( t : in out Bracket_Term ) is
  begin
    t.coeff := -t.coeff;
  end Min;

  function "-" ( p : Bracket_Polynomial ) return Bracket_Polynomial is

    res : Bracket_Polynomial;

  begin
    Copy(p,res);
    Min(res);
    return res;
  end "-";

  procedure Min ( p : in out Bracket_Polynomial ) is

    tmp : Bracket_Polynomial := p;

  begin
    while not Is_Null(tmp) loop
      declare
        bt : Bracket_Term := Head_Of(tmp);
      begin
        Min(bt);
        Set_Head(tmp,bt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Min;

  function "-" ( t : Bracket_Term; p : Bracket_Polynomial )
               return Bracket_Polynomial is

    mp : Bracket_Polynomial := -p;
    res : constant Bracket_Polynomial := t+mp;

  begin
    Clear(mp);
    return res;
  end "-";

  function "-" ( p : Bracket_Polynomial; t : Bracket_Term )
               return Bracket_Polynomial is

    mt : Bracket_Term := -t;
    res : constant Bracket_Polynomial := p+mt;

  begin
    Clear(mt);
    return res;
  end "-";

  procedure Min ( p : in out Bracket_Polynomial; t : in Bracket_Term ) is

    mt : constant Bracket_Term := -t;

  begin
    Add(p,mt);
  end Min;

  function "-" ( p,q : Bracket_Polynomial ) return Bracket_Polynomial is

    mq : Bracket_Polynomial := -q;
    res : constant Bracket_Polynomial := p+mq;

  begin
    Clear(mq);
    return res;
  end "-";

  procedure Min ( p : in out Bracket_Polynomial;
                  q : in Bracket_Polynomial ) is

    mq : constant Bracket_Polynomial := -q;

  begin
    Add(p,mq);
  end Min;

-- ITERATORS OVER MONOMIALS :

  function Number_of_Monomials ( p : Bracket_Polynomial ) return natural32 is
  begin
    return Length_Of(p);
  end Number_of_Monomials;

  procedure Enumerate_Terms ( p : in Bracket_Polynomial ) is

    tmp : Bracket_Polynomial := p;
    continue : boolean := true;

  begin
    while not Is_Null(tmp) loop
      Process(Head_Of(tmp),continue);
      exit when not continue;
      tmp := Tail_Of(tmp);
    end loop;
  end Enumerate_Terms;

-- DESTRUCTOR :

  procedure Clear ( t : in out Bracket_Term ) is
  begin
    Clear(t.monom);
  end Clear;

  procedure Clear ( p : in out Bracket_Polynomial ) is

    tmp : Bracket_Polynomial := p;

  begin
    while not Is_Null(tmp) loop
      declare
        t : Bracket_Term := Head_Of(tmp);
      begin
        Clear(t);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Bracket_Terms.Clear(Lists_of_Bracket_Terms.List(p));
  end Clear;

end Generic_Bracket_Polynomials;
