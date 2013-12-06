with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Permutations,Permute_Operations;    use Permutations,Permute_Operations;
with Standard_Linear_Product_System;

package body Linear_Symmetric_Reduction is

-- AUXILIARY DATA STRUCTURE AND OPERATIONS :

  type Lin_Sys is 
    array ( integer32 range <> ) of Standard_Complex_Vectors.Link_to_Vector;

-- ELEMENTARY OPERATIONS :

  function Linear_System ( pos : Vector ) return Lin_Sys is

  -- DESCRIPTION :
  --   Creates a linear system, by extracting the vectors that
  --   correspond to the entries in the given position.

    res : Lin_Sys(pos'range);
    use Standard_Linear_Product_System;

  begin
    for k in res'range loop
      res(k) := Get_Hyperplane(natural32(k),natural32(pos(k)));
    end loop;
    return res;
  end Linear_System;

  function Permute ( p : Permutation; ls : Lin_Sys ) return Lin_Sys is

  -- DESCRIPTION :
  --   Permutes the equations in the linear system.

    res : Lin_Sys(ls'range);
    use Standard_Complex_Vectors;

  begin
    for i in p'range loop
      if p(i) >= 0
       then res(i) := ls(p(i));
       else res(i) := -ls(-p(i));
      end if;
    end loop;
    return res;
  end Permute;

  function Permute ( ls : Lin_Sys; p : Permutation ) return Lin_Sys is

  -- DESCRIPTION :
  --   Permutes the unknowns in the linear system.

    res : Lin_Sys(ls'range);

  begin
    for k in res'range loop
      res(k) := new Standard_Complex_Vectors.Vector'(p*ls(k).all);
    end loop;
    return res;
  end Permute;

  function Permutable ( ls1,ls2 : Lin_Sys ) return boolean is

  -- DESCRIPTION :
  --   Returns true when there exists a permutation that permutes
  --   the first linear system into the second one.

    found : boolean := true;

  begin
    for i in ls1'range loop
      for j in ls2'range loop
        found := Permutable(ls1(i).all,ls2(j).all);
        exit when found;
      end loop;
      exit when not found;
    end loop;
    return found;
  end Permutable;

  function Sign_Permutable ( ls1,ls2 : Lin_Sys ) return boolean is

  -- DESCRIPTION :
  --   Returns true when there exists a permutation that permutes 
  --   the first linear system into the second one, also w.r.t. sign
  --   permutations.

    found : boolean := true;

  begin
    for i in ls1'range loop
      for j in ls2'range loop
        found := Sign_Permutable(ls1(i).all,ls2(j).all);
        exit when found;
      end loop;
      exit when not found;
    end loop;
    return found;
  end Sign_Permutable;

--  procedure Clear ( ls : in out Lin_Sys ) is
--
--  -- DESCRIPTION :
--  --   Deallocation of the occupied memory space.
--
--  begin
--    for k in ls'range loop
--      Standard_Complex_Vectors.Clear(ls(k));
--    end loop;
--  end Clear;

-- UTILITIES :

  procedure Search_Permutable 
                  ( sub : in Lin_Sys; pos : in Vector;
                    res,res_last : in out List ) is

  -- DESCRIPTION :
  --   In the list of positions, already in res, it will be searched
  --   whether there exists a linear system that is permutable with the
  --   given linear system.

    tmp : List := res;
    found : boolean := false;
    ls2 : Lin_Sys(sub'range);

  begin
    while not Is_Null(tmp) loop
      ls2 := Linear_System(Head_Of(tmp).all);
      found := Permutable(sub,ls2);
      exit when found;
      tmp := Tail_Of(tmp);
    end loop;
    if not found
     then Append(res,res_last,pos);
    end if;
  end Search_Permutable;

  procedure Search_Sign_Permutable 
                  ( sub : in Lin_Sys; pos : in Vector;
                    res,res_last : in out List ) is

  -- DESCRIPTION :
  --   In the list of positions, already in res, it will be searched
  --   whether there exists a linear system that is sign permutable
  --   with the given linear system.

    tmp : List := res;
    found : boolean := false;
    ls2 : Lin_Sys(sub'range);

  begin
    while not Is_Null(tmp) loop
      ls2 := Linear_System(Head_Of(tmp).all);
      found := Sign_Permutable(sub,ls2);
      exit when found;
      tmp := Tail_Of(tmp);
    end loop;
    if not found
     then Append(res,res_last,pos);
    end if;
  end Search_Sign_Permutable;

  function Search_Position ( sub : Lin_Sys ) return Vector is

  -- DESCRIPTION :
  --   Returns the position of the system in the product system.

    res : Vector(sub'range);
    lh : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for k in 1..Standard_Linear_Product_System.Dimension loop
      res(integer32(k)) := 0;
      for L in 1..Standard_Linear_Product_System.Number_of_Hyperplanes(k) loop
        lh := Standard_Linear_Product_System.Get_Hyperplane(k,L);
        if Standard_Complex_Vectors.Equal(sub(integer32(k)).all,lh.all)
         then res(integer32(k)) := integer32(L);
        end if;
        exit when res(integer32(k)) /= 0;
      end loop;
    end loop;
    return res;
  end Search_Position;

  procedure Permute_and_Search 
                ( v,w : List_of_Permutations; sub : in Lin_Sys;
                  pos : in Vector; res,res_last : in out List ) is

  -- DESCRIPTION :
  --   The permutations are applied to the subsystem.
  --   If none of the positions of the permuted systems already
  --   belongs to res, then its position pos will be added to res.

    lv,lw : List_of_Permutations;
    found : boolean := false;

  begin
    lv := v;  lw := w;
   -- put_line("The permuted positions : ");
    while not Is_Null(lv) loop
      declare
        vpersub : constant Lin_Sys(sub'range)
                := Permute(sub,Permutation(Head_Of(lv).all));
        wpersub : constant Lin_Sys(sub'range) 
                := Permute(Permutation(Head_Of(lw).all),vpersub);
        perpos : constant Vector(pos'range) := Search_Position(wpersub);
      begin
        if Is_In(res,perpos)
         then found := true;
        end if;
      end;
      exit when found;
      lv := Tail_Of(lv);
      lw := Tail_Of(lw);
    end loop;
    if not found
     then Append(res,res_last,pos);
    end if;
  end Permute_and_Search;

  function Generate_Positions return List is

    res,res_last : List;
    use Standard_Linear_Product_System;
    n : constant integer32 := integer32(Dimension);
    pos : Vector(1..n) := (1..n => 1);

    procedure Generate_Positions ( k : natural32 ) is
    begin
      if integer32(k) > n then
        Append(res,res_last,pos);
      else
        for ell in 1..Number_of_Hyperplanes(k) loop
          pos(integer32(k)) := integer32(ell);
          Generate_Positions(k+1);
        end loop;
      end if;
    end Generate_Positions;

  begin
    Generate_Positions(1);
    return res;
  end Generate_Positions;

-- TARGET ROUTINES :

  function Linear_Symmetric_Reduce ( sign : boolean ) return List is

    res : List;

  begin
    res := Generate_Positions;
    Linear_Symmetric_Reduce(res,sign);
    return res;
  end Linear_Symmetric_Reduce;

  function Linear_Symmetric_Reduce 
              ( v,w : List_of_Permutations ) return List is

    res : List;

  begin
    res := Generate_Positions;
    Linear_Symmetric_Reduce(v,w,res);
    return res;
  end Linear_Symmetric_Reduce;

  procedure Linear_Symmetric_Reduce ( lp : in out List; sign : in boolean ) is

    res,res_last : List;
    sub : Lin_Sys(1..integer32(Standard_Linear_Product_System.Dimension));
    pos : Vector(sub'range);
    tlp : List := lp;

  begin
    while not Is_Null(tlp) loop
      pos := Head_Of(tlp).all;
      sub := Linear_System(pos);
      if not sign
       then Search_Permutable(sub,pos,res,res_last);
       else Search_Sign_Permutable(sub,pos,res,res_last);
      end if;
      tlp := Tail_Of(tlp);
    end loop;
    Clear(lp);
    lp := res;
  end Linear_Symmetric_Reduce;

  procedure Linear_Symmetric_Reduce
              ( v,w : in List_of_Permutations; lp : in out List ) is

    res,res_last : List;
    sub : Lin_Sys(1..integer32(Standard_Linear_Product_System.Dimension));
    pos : Vector(sub'range);
    tlp : List := lp;

  begin
    while not Is_Null(tlp) loop
      pos := Head_Of(tlp).all;
      sub := Linear_System(pos);
      Permute_and_Search(v,w,sub,pos,res,res_last);
      tlp := Tail_Of(tlp);
    end loop;
    Clear(lp);
    lp := res;
  end Linear_Symmetric_Reduce;

end Linear_Symmetric_Reduction;
