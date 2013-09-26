with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Random_Numbers;
with Permutations;                       use Permutations;
with Symmetry_Group,Symmetry_Group_io;   use Symmetry_Group,Symmetry_Group_io;

procedure ts_group is

-- DESCRITPION :
--   Test on permutation groups.

  procedure Search ( g : in List_of_Permutations; i,j : in integer32;
                     fail : out boolean; s : out Permutation ) is

  -- DESCRIPTION :
  --   Searches the list g for a permutation s for which s(i) = j. 
  --   If such an s exists, then fail is false and s is on return,
  --   otherwise fail is true and s has no meaning.

    tmp : List_of_Permutations := g;
    lkp : Link_to_Permutation;

  begin
    fail := true;
    while not Is_Null(tmp) loop
      lkp := Head_Of(tmp);
      if lkp(i) = j
       then fail := false; s := Permutation(lkp.all); exit;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Search;

  procedure Sort ( g : in List_of_Permutations; p : in Permutation ) is

  -- DESCRIPTION :
  --   Attempts to sort the permutation p, i.e.: reduce to the
  --   identity permutation by application of the permutations in p.

    q : Permutation(p'range) := p;
    s : Permutation(p'range);
    fail : boolean;
    ind : integer32 := q'first;

  begin
    put("Attempting to sort "); put(q); put_line(" ...");
    while ind <= q'last loop
      if q(ind) = ind then
        ind := ind + 1;
      else
        Search(g,ind,p(ind),fail,s);
        if not fail then
          put("applying "); put(s); put(" to "); put(q);
          q := s*q;
          put(" : "); put(q); new_line;
         -- ind := q'first;
        else
          ind := q'last + 1;
        end if;
      end if;
    end loop;
  end Sort;

  function Random_Permutation ( n : integer32 ) return Permutation is

  -- DESCRIPTION :
  --   Returns a random permutation on n variables,
  --   obtained after n random swaps.

    res : Permutation(1..n);
    ind,tmp : integer32;

  begin
    for i in 1..n loop
      res(i) := i;
    end loop;
    for i in 1..n loop
      ind := Standard_Random_Numbers.Random(1,n);
      tmp := res(i);
      res(i) := res(ind);
      res(ind) := tmp; 
    end loop;
    return res;
  end Random_Permutation;

  function Sorting_Swaps ( p : Permutation ) return List_of_Permutations is

  -- DESCRIPTION :
  --   Returns the list of swaps used to sort the elements in p
  --   in ascending order.

    res,res_last : List_of_Permutations;
    swap : Permutation(p'range);
    work : Permutation(p'range) := p;

  begin
    for i in swap'range loop
      swap(i) := i;
    end loop;
    for i in work'range loop
      if work(i) /= i then
        for j in i+1..work'last loop
          if work(j) = i then
            swap(i) := j; swap(j) := i;
            Append(res,res_last,swap);
            swap(i) := i; swap(j) := j;
            put(work); put(" ->");
            work(j) := work(i); work(i) := i;
            put(work); new_line;
          end if;
        end loop;
      end if;
    end loop;
    put("the sorted permutation : "); put(work); new_line;
    return res;
  end Sorting_Swaps;

  function Unsort ( n : integer32; s : List_of_Permutations )
                  return Permutation is

  -- DESCRIPTION :
  --   Applies the swaps of s to the identity permutation of dimension n.
  --   The swaps are applied in reserse order.
  --   The resulting permutation is returned.
  --   This function is a check on Sorting_Swaps, as
  --   we must have that Unsort(Sorting_Swaps(p)) = p.

    res : Permutation(1..n);
 
    procedure Apply ( t : in List_of_Permutations ) is

    -- DESCRIPTION :
    --   Applies the swaps in t to res.

    begin
      if Length_Of(t) > 1 then
        Apply(Tail_Of(t));
        put(res); put(" ->");
        res := res*Permutation(Head_Of(t).all);
        put(res); new_line;
      else
        put(res); put(" ->");
        res := res*Permutation(Head_Of(t).all);
        put(res); new_line;
      end if;
    end Apply;
  
  begin
    for i in 1..n loop
      res(i) := i;
    end loop;
    if Length_Of(s) > 0
     then Apply(s);
    end if;
    return res;
  end Unsort;

  procedure Main is

    n : integer32 := 0;
    g : List_of_Permutations;

  begin
    put("give the number of variables : "); get(n);
    g := SymGrp(n);
    put("the transpositions for n = "); put(n,1); put_line(" :");
    put(g);
    declare
      p,q : Permutation(1..n);
      ans : character;
      s : List_of_Permutations;
    begin
      put("Random permutation ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then p := Random_Permutation(n);
       else put("Give a permutation : "); get(p);
      end if;
      put("-> our permutation : "); put(p);
      if Is_In(g,p) then
        put_line(" is a generating transposition");
      else
        put_line(" is not a generating transposition");
       -- Sort(g,p);
        s := Sorting_Swaps(p);
        put_line("The list of sorting swaps : "); put(s);
        put_line("Unsorting the identity ...");
        q := Unsort(n,s);
        put("p = "); put(p); new_line;
        put("q = "); put(q); new_line;
        if Equal(p,q)
         then put_line("Sorting swaps are okay.");
         else put_line("Sorting swaps are NOT okay, bug?!");
        end if;
      end if;
    end;
  end Main;

begin
  Main;
end ts_group;
