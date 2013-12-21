with Standard_Integer_Vectors;           use Standard_Integer_Vectors;

package body Symmetry_Group is

  use Lists_of_Permutations;

  procedure Add ( L : in out List_of_Permutations; p : in Permutation ) is

    lp : Link_To_Permutation;

  begin
    lp := new Vector'(Vector(p));
    Construct(lp,l);
  end Add;

  procedure Append ( first,last : in out List_of_Permutations;
                     p : in Permutation ) is

    lp : Link_To_Permutation;

  begin
    lp := new Vector'(Vector(p));
    Append(first,last,lp);
  end Append;

  function Union ( a,b : List_of_Permutations ) return List_of_Permutations is

    tmp,res : List_of_Permutations;

  begin
    res := a;
    tmp := b;
    while not Is_Null(tmp) loop
      Add(res,Permutation(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Union;

  function SymGrp ( n : integer32 ) return List_of_Permutations is

    sn : List_of_Permutations;
    p : Vector(1..n);

  begin
    for k in p'range loop
      p(k) := k;
    end loop;
    for k in reverse 1..n loop
      p(k) := 1; p(1) := k;
      declare
        lp : constant Link_to_Permutation := new Vector'(p);
      begin
        Construct(lp,sn);
      end;
      p(1) := 1; p(k) := k;
    end loop;
    return sn;
  end SymGrp;

  function Generate ( g : List_of_Permutations ) return List_of_Permutations is

    res : List_Of_Permutations;
    at_end : boolean;

  -- IMPORTANT :
  --   This routine assumes that permutations are added to the front
  --   of the list !!!

  begin
    if not Is_Null(g) then
      declare
        p1,p2,r : Permutation(Head_Of(g).all'range);
        temp1,temp2,nwres : List_Of_Permutations;
        nb,cnt : natural32;
      begin
        temp1 := g; -- INITIALIZE res :
        while not Is_Null(temp1) loop
          p1 := Permutation(Head_Of(temp1).all);
          Add(res,p1);
          temp1 := Tail_Of(temp1);
        end loop;
        at_end := false; -- CONSTRUCT res :
        nb := Length_Of(res);
        while not at_end loop
          temp1 := g; --res;
          while not Is_Null(temp1) loop
            p1 := Permutation(Head_Of(temp1).all);
            cnt := 0;
            temp2 := res;
            while cnt < nb loop
              p2 := Permutation(Head_Of(temp2).all);
              r := p1*p2;
              if not Is_In(res,r) and not Is_In(nwres,r)
               then Add(nwres,r);
              end if;
              cnt := cnt + 1;
              temp2 := Tail_Of(temp2);
            end loop;
            temp1 := Tail_Of(temp1);
          end loop;
          nb := Length_Of(nwres);
          at_end := (nb = 0);
          if not at_end then
            temp2 := nwres;
            while not Is_Null(temp2) loop
              Add(res,Permutation(Head_Of(temp2).all));
              temp2 := Tail_Of(temp2);
            end loop;
            Clear(nwres);
          end if;
        end loop;
      end;
    end if;
    return res;
  end Generate;

-- SELECTORS :

  function Number ( L : List_Of_Permutations ) return natural32 is
  begin
    return Length_Of(L);
  end Number;

  function Is_In ( L : List_of_Permutations; p : Permutation )
                 return boolean is

    tmp : List_Of_Permutations := L;

  begin
    while not Is_Null(tmp) loop
      if Equal(Permutation(Head_Of(tmp).all),p)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  procedure Iterator ( L : in List_of_Permutations ) is

    tmp : List_Of_Permutations := L;
    cont : boolean;

  begin
    cont := false;
    while not Is_Null(tmp) loop
      Process(Permutation(Head_Of(tmp).all),cont);
      exit when not cont;
      tmp := Tail_Of(tmp);
    end loop;
  end Iterator;

-- DESTRUCTOR :

  procedure Clear ( L : in out List_of_Permutations ) is
  begin
    Lists_of_Permutations.Clear(Lists_of_Permutations.List(L));
  end Clear;

end Symmetry_Group;
