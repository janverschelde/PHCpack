with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Norms;             use Standard_Integer_Norms;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Standard_Integer64_Linear_Solvers;  use Standard_Integer64_Linear_Solvers;

package body Lists_of_Vectors64_Utilities is

  procedure Compute_Normal ( v : in VecVec; n : out Link_to_Vector;
                             deg : out natural64 ) is

    d : Link_to_Vector renames v(v'last);
    im : matrix(d'range,d'range);
    res : Link_to_Vector;
    cnt : integer64;

  begin
    res := new Vector(d'range);
    cnt := integer64(im'first(1));
    for i in v'first..(v'last-1) loop
      for j in im'range(2) loop
        im(integer32(cnt),j) := v(i)(j) - d(j);
      end loop;
      cnt := cnt + 1;
    end loop;
    for i in integer32(cnt)..im'last(1) loop
      for j in im'range(2) loop
        im(i,j) := 0;
      end loop;
    end loop;
    Upper_Triangulate(im);
    cnt := 1;
    for k in im'first(1)..im'last(1)-1 loop
      cnt := cnt*im(k,k);
    end loop;
    if cnt < 0
     then deg := natural64(-cnt);
     else deg := natural64(cnt);
    end if;
    Scale(im);
    Solve0(im,res.all);
    Normalize(res.all);
    n := res;
  end Compute_Normal;

  function Compute_Normal ( v : VecVec ) return Link_to_Vector is

    deg : natural64;
    res : Link_to_Vector;

  begin
    Compute_Normal(v,res,deg);
    return res;
  end Compute_Normal;

  function Pointer_to_Last ( L : List ) return List is

    res : List := L;

  begin
    if not Is_Null(res) then
      while not Is_Null(Tail_Of(res)) loop
        res := Tail_Of(res);
      end loop;
    end if;
    return res;
  end Pointer_to_Last;

  procedure Move_to_Front ( L : in out List; v : in Vector ) is

    tmp : List := L;
    found : boolean := false;
    first,lv : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      if Equal(lv.all,v)
       then found := true;
       else tmp := Tail_Of(tmp);
      end if;
      exit when found;
    end loop;
    if found then
      first := Head_Of(l);
      if first /= lv then
        lv.all := first.all;  Set_Head(tmp,lv);
        first.all := v;       Set_Head(l,first);
      end if;
    end if;
  end Move_to_Front;

  function Difference ( L1,L2 : List ) return List is

    res,res_last : List;
    tmp : List := L1;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if not Is_In(L2,pt.all)
       then Append(res,res_last,pt.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Difference;

  function Different_Points ( L : List ) return List is

    tmp,res,res_last : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      Append_Diff(res,res_last,Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Different_Points;

  procedure Remove_Duplicates ( L : in out List ) is

    res : constant List := Different_Points(L);

  begin
    Deep_Clear(L);
    L := res;
  end Remove_Duplicates;

end Lists_of_Vectors64_Utilities;
