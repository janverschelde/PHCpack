with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;

package body Standard_Multiple_Solutions is

  function Equal ( n : natural32; tol : double_float;
                   s1,s2 : Vector ) return boolean is
  begin
    for i in 1..integer32(n) loop
      if AbsVal(s1(i)-s2(i)) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  procedure Set_Multiplicity
              ( sols : in out Solution_List; s : in Solution;
                tol : in double_float; n,m : in natural32 ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if Equal(n,tol,ls.v,s.v)
       then ls.m := integer32(m); Set_Head(tmp,ls);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Set_Multiplicity;

  function Number_of_Occurrences 
              ( sols : Solution_List; s : Solution;
                tol : in double_float; n : in natural32 ) return natural32 is

    res : natural32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if Equal(n,tol,ls.v,s.v)
       then res := res+1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Number_of_Occurrences;

  function Is_In ( sols : Solution_List; v : Vector;
                   tol : double_float; n : natural32 ) return boolean is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if Equal(n,tol,ls.v,v)
       then return true;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return false;
  end Is_In;

  procedure Remove_Duplicates
              ( sols : in out Solution_List;
                tol : in double_float; n : in natural32 ) is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if not Is_In(res,ls.v,tol,n)
       then Append(res,res_last,ls.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(sols);
    sols := res;
  end Remove_Duplicates;

  procedure Merge_Multiple_Solutions
              ( sols : in out Solution_List; tol : in double_float ) is

    res,res_last,pt1,pt2 : Solution_List;
    ls1,ls2 : Link_to_Solution;
    found : boolean;

  begin
    pt1 := sols;
    while not Is_Null(pt1) loop
      ls1 := Head_Of(pt1);
      found := false;
      pt2 := res;
      while not Is_Null(pt2) loop
        ls2 := Head_Of(pt2);
        found := Equal(natural32(ls1.n),tol,ls1.v,ls2.v);
        exit when found;
        pt2 := Tail_Of(pt2);
      end loop;
      if not found then
        Append(res,res_last,ls1.all);
      elsif ls1.m > ls2.m then
        ls2.m := ls1.m;
        Set_Head(pt2,ls2);
      end if;
      pt1 := Tail_of(pt1);
    end loop;
    Clear(sols);
    sols := res;
  end Merge_Multiple_Solutions;

  procedure Compute_Multiplicities
              ( sols : in out Solution_List;
                tol : in double_float; n : in natural32 ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    m : natural32;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      m := Number_of_Occurrences(sols,ls.all,tol,n);
      Set_Multiplicity(sols,ls.all,tol,n,m);
      tmp := Tail_Of(tmp);
    end loop;
  end Compute_Multiplicities;

  procedure Compute_Multiplicities
              ( nd : in out Standard_Deflation_Trees.Node;
                tol : in double_float; n : in natural32 ) is

    use Standard_Deflation_Trees;

  begin
    if not Is_Null(nd.sols)
     then Compute_Multiplicities(nd.sols,tol,n);
          Remove_Duplicates(nd.sols,tol,n);
    end if;
    for i in nd.children'range loop
      if nd.children(i) /= null
       then Compute_Multiplicities(nd.children(i).all,tol,n);
      end if;
    end loop;
  end Compute_Multiplicities;

end Standard_Multiple_Solutions;
