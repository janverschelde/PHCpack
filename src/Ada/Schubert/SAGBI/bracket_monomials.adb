package body Bracket_Monomials is

-- INTERNAL SORTING ROUTINE TO MAINTAIN ORDER :

  procedure Swap_Heads ( bm1,bm2 : in out Bracket_Monomial;
                         lb1,lb2 : in out Link_to_Bracket ) is

    b1 : constant Bracket(lb1'range) := lb1.all;
    b2 : constant Bracket(lb2'range) := lb2.all;

  begin
    Clear(lb2); lb2 := new Bracket'(b1);
    Clear(lb1); lb1 := new Bracket'(b2);
    Set_Head(bm1,lb1);
    Set_Head(bm2,lb2);
  end Swap_Heads;

  procedure Sort ( bm : in out Bracket_Monomial ) is

    tmp1 : Bracket_Monomial := bm;

  begin
    while not Is_Null(tmp1) loop
      declare
        lb1 : Link_to_Bracket := Head_Of(tmp1);
        min : Link_to_Bracket := lb1;
        mintmp : Bracket_Monomial := tmp1;
        tmp2 : Bracket_Monomial := Tail_Of(tmp1);
      begin
        while not Is_Null(tmp2) loop
          if Head_Of(tmp2).all < min.all then
            min := Head_Of(tmp2);
            mintmp := tmp2;
          end if;
          tmp2 := Tail_Of(tmp2);
        end loop;
        if not Is_Equal(lb1.all,min.all)
         then Swap_Heads(tmp1,mintmp,lb1,min);
        end if;
      end;
      tmp1 := Tail_Of(tmp1);
    end loop;
  end Sort;

-- CONSTRUCTORS :

  function Create ( b : Bracket ) return Bracket_Monomial is

    bm : Bracket_Monomial;
    lb : constant Link_to_Bracket := new Bracket'(b);

  begin
    Construct(lb,bm);
    return bm;
  end Create;

  procedure Multiply ( bm : in out Bracket_Monomial; b : in Bracket ) is
  begin
    if Is_Null(bm) then
      bm := Create(b);
    else
      declare
        lb : constant Link_to_Bracket := new Bracket'(b);
      begin
        Construct(lb,bm);
        Sort(bm);
      end;
    end if;
  end Multiply;

  procedure Append ( bm : in out Bracket_Monomial; b : in Bracket ) is
  begin
    if Is_Null(bm) then
      bm := Create(b);
    else
      declare
        lb : constant Link_to_Bracket := new Bracket'(b);
        last : Bracket_Monomial := bm;
        tmp : Bracket_Monomial := Tail_Of(last);
      begin
        while not Is_Null(tmp) loop
          last := tmp;
          tmp := Tail_Of(tmp);
        end loop;
        Append(bm,last,lb);
      end;
    end if;
  end Append;

  procedure Copy_Multiply ( bm1 : in Bracket_Monomial;
                            bm2 : in out Bracket_Monomial ) is

    tmp : Bracket_Monomial := bm1;

  begin
    Clear(bm2);
    while not Is_Null(tmp) loop
      declare
        b : constant Bracket := Head_Of(tmp).all;
      begin
        Multiply(bm2,b);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Copy_Multiply;

  procedure Copy_Append ( bm1 : in Bracket_Monomial;
                          bm2 : in out Bracket_Monomial ) is

    tmp : Bracket_Monomial := bm1;
    last : Bracket_Monomial;

  begin
    Clear(bm2); last := bm2;
    while not Is_Null(tmp) loop
      declare
        b : constant Bracket := Head_Of(tmp).all;
        lb : constant Link_to_Bracket := new Bracket'(b);
      begin
        Append(bm2,last,lb);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Copy_Append;

  function Create ( bm : Bracket_Monomial ) return Array_of_Brackets is

    res : Array_of_Brackets(1..integer32(Length_Of(bm)));
    tmp : Bracket_Monomial := bm;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

-- OPERATIONS :

  function "*" ( b1,b2 : Bracket ) return Bracket_Monomial is

    bm : Bracket_Monomial := Create(b1);

  begin
    Multiply(bm,b2);
    return bm;
  end "*";

  function "*" ( bm : Bracket_Monomial; b : Bracket )
               return Bracket_Monomial is

    res : Bracket_Monomial;

  begin
    Copy(bm,res);
    Multiply(res,b);
    return res;
  end "*";

  function "*" ( b : Bracket; bm : Bracket_Monomial )
               return Bracket_Monomial is

    res : Bracket_Monomial;

  begin
    Copy(bm,res);
    Multiply(res,b);
    return res;
  end "*";

  function "*" ( bm1,bm2 : Bracket_Monomial ) return Bracket_Monomial is

    res : Bracket_Monomial;

  begin
    Copy(bm1,res);
    Multiply(res,bm2);
    return res;
  end "*";

  procedure Multiply ( bm1 : in out Bracket_Monomial;
                       bm2 : in Bracket_Monomial ) is

    tmp : Bracket_Monomial := bm2;

  begin
    while not Is_Null(tmp) loop
      declare
        b : constant Bracket := Head_Of(tmp).all;
      begin
        Multiply(bm1,b);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Multiply;

  function Is_Null ( bm : Bracket_Monomial ) return boolean is
  begin
    return Lists_of_Brackets.Is_Null(Lists_of_Brackets.List(bm));
  end Is_Null;

  function Is_Equal ( bm1,bm2 : Bracket_Monomial ) return boolean is

    tmp1 : Bracket_Monomial := bm1;
    tmp2 : Bracket_Monomial := bm2;

  begin
    if Length_Of(tmp1) /= Length_Of(tmp2) then
      return false;
    else
      while not Is_Null(tmp1) loop
        if not Is_Equal(Head_Of(tmp1).all,Head_Of(tmp2).all) then
          return false;
        else
          tmp1 := Tail_Of(tmp1);
          tmp2 := Tail_Of(tmp2);
        end if;
      end loop;
      return true;
    end if;
  end Is_Equal;

  function Dimension ( bm : Bracket_Monomial ) return natural32 is
  begin
    if Is_Null(bm)
     then return 0;
     else return natural32(Head_Of(bm)'last);
    end if;
  end Dimension;
   
  function "<" ( bm1,bm2 : Bracket_Monomial ) return boolean is

    tmp1 : Bracket_Monomial := bm1;
    tmp2 : Bracket_Monomial := bm2;
    lb1,lb2 : Link_to_Bracket;

  begin
    while not Is_Null(tmp1) and not Is_Null(tmp2) loop
      lb1 := Head_Of(tmp1); lb2 := Head_Of(tmp2);
      if lb1.all < lb2.all then
        return true;
      elsif lb1.all > lb2.all then
        return false;
      else 
        tmp1 := Tail_Of(tmp1); tmp2 := Tail_Of(tmp2);
      end if;
    end loop;
    if Is_Null(tmp1) and not Is_Null(tmp2)
     then return true;
     else return false;
    end if;
  end "<";

  function ">" ( bm1,bm2 : Bracket_Monomial ) return boolean is

    tmp1 : Bracket_Monomial := bm1;
    tmp2 : Bracket_Monomial := bm2;
    lb1,lb2 : Link_to_Bracket;

  begin
    while not Is_Null(tmp1) and not Is_Null(tmp2) loop
      lb1 := Head_Of(tmp1); lb2 := Head_Of(tmp2);
      if lb1.all > lb2.all then
        return true;
      elsif lb1.all < lb2.all then
        return false;
      else
        tmp1 := Tail_Of(tmp1); tmp2 := Tail_Of(tmp2);
      end if;
    end loop;
    if Is_Null(tmp2) and not Is_Null(tmp1)
     then return true;
     else return false;
    end if;
  end ">";

  function Divisible ( bm : Bracket_Monomial; b : Bracket ) return boolean is

    tmp : Bracket_Monomial := bm;

  begin
    while not Is_Null(tmp) loop
      if Is_Equal(Head_Of(tmp).all,b)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Divisible;

-- ITERATORS OVER THE BRACKETS :

  function Number_of_Brackets ( bm : Bracket_Monomial ) return natural32 is
  begin
    return Length_Of(bm);
  end Number_of_Brackets;

  procedure Enumerate_Brackets ( bm : in Bracket_Monomial ) is

    tmp : Bracket_Monomial := bm;
    continue : boolean := true;

  begin
    while not Is_Null(tmp) loop
      Process(Head_Of(tmp).all,continue);
      exit when not continue;
      tmp := Tail_Of(tmp);
    end loop;
  end Enumerate_Brackets;

-- DESTRUCTOR :

  procedure Clear ( bm : in out Bracket_Monomial ) is

    tmp : Bracket_Monomial := bm;
    lb : Link_to_Bracket;

  begin
    while not Is_Null(tmp) loop
      lb := Head_Of(tmp);
      Clear(lb);
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Brackets.Clear(Lists_of_Brackets.List(bm));
  end Clear;

end Bracket_Monomials;
