with unchecked_deallocation;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Permute_Operations;                 use Permute_Operations;

package body Orbits_of_Solutions is

-- AUXILIAIRIES :

  function Equal_upon_Sign
             ( x,y : Standard_Complex_Vectors.Vector; tol : double_float )
             return boolean is

  -- DESCRIPTION :
  --   Returns true if for all i : | x(i) +- y(i) | <= tol.

  begin
    for i in x'range loop
      if (AbsVal(x(i)-y(i)) > tol) and then (AbsVal(x(i)+y(i)) > tol)
       then return false;
      end if;
    end loop;
    return true;
  end Equal_upon_Sign;

  procedure Is_In ( sols : in out Solution_List; s : in Solution;
                    sign : in boolean; tol : in double_float;
                    isin : out boolean ) is

  -- DESCRIPTION :
  --   The output variable isin isin becomes true if one of the solutions
  --   in the list sols equals s.  When this is the case, the multiplicity
  --   of the corresponding solution is augmented.
  --   Otherwise, isin equals false on return.

    s1 : Solution(s.n);
    temp : Solution_List;

  begin
    temp := sols;
    while not Is_Null(temp) loop
      s1 := Head_Of(temp).all;
      if (sign and then Equal_upon_Sign(s.v,s1.v,tol))
        or else (not sign and then Equal(s.v,s1.v,tol))
       then Head_Of(temp).m := Head_Of(temp).m + 1;
	    isin := true;
	    return;
      end if;
      temp := Tail_Of(temp);
    end loop;
    isin := false;
  end Is_In;

-- CONSTRUCTORS :

  procedure Analyze ( l : in List_of_Permutations; sign : in boolean;
                      tol : in double_float; sols : in out Solution_List ) is
  begin
    if not Is_Null(sols) then
      declare
        res,res_last,tmpsols : Solution_List;
        n : constant integer32 := Head_Of(sols).n;
        s1,s2 : Solution(n);
        lp : Link_to_Permutation;
        tmpl : List_of_Permutations;
        found : boolean;
      begin
        tmpsols := sols;
        while not Is_Null(tmpsols) loop
          s1 := Head_Of(tmpsols).all;
          tmpl := l;
          while not Is_Null(tmpl) loop
            lp := Head_Of(tmpl);
            s2.t := s1.t;
            s2.m := 1;
            s2.v := Permutation(lp.all)*s1.v;
            Is_In(res,s2,sign,tol,found);
            exit when found;
            tmpl := Tail_Of(tmpl);
          end loop;
          if not found 
           then Append(res,res_last,s1);
          end if;
          tmpsols := Tail_Of(tmpsols);
        end loop;
        Clear(sols); sols := res;
      end;
    end if;
  end Analyze;

  function Permutable ( s : Solution; sign : boolean; tol : double_float;
                        sols : Solution_List ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the solutions s can be permuted into one of the
  --   solutions of the given list.

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      if sign and then Sign_Permutable(s.v,Head_Of(tmp).v,tol)
       then return true;
       elsif Permutable(s.v,Head_Of(tmp).v,tol)
           then return true;
           else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Permutable;

  procedure Permutable
                ( s : in Solution; sign : in boolean; tol : in double_float;
                  sols : in out Solution_List; isin : out boolean ) is

  -- DESCRIPTION :
  --   The output variable isin becomes true when the given solution can be
  --   permuted into one of the solutions of the given list.
  --   In this case, also the multiplicity of the corresponding solutions is
  --   augmented.  Otherwise, isin is false on return.

    tmp : Solution_List := sols;
    ls1 : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls1 := Head_Of(tmp);
      if (sign and then Sign_Permutable(s.v,ls1.v,tol))
        or else (not sign and then Permutable(s.v,ls1.v,tol))
       then isin := true;
            ls1.m := ls1.m + 1; Set_Head(tmp,ls1);
            return;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    isin := false;
  end Permutable;

  function Generating ( sols : Solution_List; sign : boolean;
                        tol : double_float ) return Solution_List is

    tmp,gensols,gensols_last : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        isin : boolean;
      begin
        Permutable(ls.all,sign,tol,gensols,isin);
        if not isin
         then Append(gensols,gensols_last,ls.all);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return gensols;
  end Generating;

  procedure Orbit_Structure ( s : in Solution; tol : in double_float;
			      orbit : in out Permutation;
			      nbdiff : out natural32 ) is

    nb : natural32 := 1;
    found : boolean;

  begin
    orbit(1) := integer32(nb);
    for i in 2..s.n loop
      found := false;
      for j in 1..(i-1) loop
	if AbsVal(s.v(i) - s.v(j)) < tol then
	  orbit(i) := orbit(j);
          found := true;
	  exit;
        end if;
      end loop;
      if not found then
        nb := nb + 1;
        orbit(i) := integer32(nb);
      end if;
    end loop;
    nbdiff := nb;
  end Orbit_Structure;

  function Orbits ( sols : Solution_List; tol : double_float )
                  return permutation is

    ind : integer32;
    tmp : Solution_List := sols;
    n : constant integer32 := Head_Of(sols).n;
    o,orb : Permutation(1..n);
    sol : Solution(n);

  begin
    orb := (1..n => 0);
    o := orb;
    while not Is_Null(tmp) loop
      sol := Head_Of(tmp).all;
      Orbit_Structure(sol,tol,o,natural32(ind));
      orb(ind) := orb(ind) + sol.m;
      tmp := Tail_Of(tmp);
    end loop;
    return orb;
  end Orbits;

  procedure Orbits ( grp : in List_of_Permutations; tol : in double_float;
		     sols : in out Solution_List;
		     lorb : in out List_of_Orbits ) is

    tmp : Solution_List;
    n : constant integer32 := Head_Of(sols).n;

  begin
    Analyze(grp,false,tol,sols);
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
	sol : constant Solution(n) := Head_Of(tmp).all;
	orbi : Orbit(n);
      begin
	Orbit_Structure(sol,tol,orbi.orb,orbi.nbdiff);
        orbi.nbgen := 1;
	orbi.nbsols := natural32(sol.m);
	declare
	  tmp2 : List_of_Orbits := lorb;
	  found : boolean := false;
        begin
	  while not Is_Null(tmp2) loop
	    declare
	      orbi2 : constant Orbit(n) := Head_Of(tmp2).all;
	      lorbi2 : Link_to_Orbit := Head_Of(tmp2);
	    begin
	      if orbi2.nbdiff = orbi.nbdiff 
	        and then Same_Orbit(orbi.orb,orbi2.orb)
               then lorbi2.nbgen := orbi2.nbgen + 1;
                    lorbi2.nbsols := orbi2.nbsols + natural32(sol.m);
                    found := true;
		    exit;
               else tmp2 := Tail_Of(tmp2);
              end if;
	    end;
	  end loop;
          if not found then
	    declare
              liorbi : constant Link_to_Orbit := new Orbit'(orbi);
            begin
              Construct(liorbi,lorb);
            end;
          end if;
        end;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Orbits;

-- SELECTOR :

  function Same_Orbit ( orb1,orb2 : Permutation ) return boolean is

    f1,f2 : Permutation(1..orb1'last);
    found : boolean;

  begin
   -- construct frequencies :
    f1 := (f1'range => 0);
    for i in orb1'range loop
      f1(orb1(i)) := f1(orb1(i)) + 1;
    end loop;
    f2 := (f2'range => 0);
    for i in orb2'range loop
      f2(orb2(i)) := f2(orb2(i)) + 1;
    end loop;
   -- compare the frequencies :
    for i in f1'range loop
      found := false;
      for j in f2'range loop
	if f1(i) = f2(j)
	 then found := true;
	      exit;
        end if;
      end loop;
      if not found
       then return false;
      end if;
    end loop;
    for i in f2'range loop
      found := false;
      for j in f1'range loop
	if f2(i) = f1(j)
	 then found := true;
	      exit;
        end if;
      end loop;
      if not found
       then return false;
      end if;
    end loop;
    return true;
  end Same_Orbit;

-- DESTRUCTOR :

  procedure Clear ( lorb : in out List_of_Orbits ) is

    procedure free is new unchecked_deallocation(Orbit,Link_to_Orbit);

    tmp : List_of_Orbits := lorb;
  begin
    while not Is_Null(tmp) loop
      declare
	lior : Link_to_Orbit := Head_Of(tmp);
      begin
	free(lior);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Orbits.Clear(Lists_of_Orbits.List(lorb));
  end Clear;

end Orbits_of_Solutions;
