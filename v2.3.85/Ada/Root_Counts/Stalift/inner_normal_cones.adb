with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Standard_Integer32_Transformations; 
 use Standard_Integer32_Transformations;
with Integer32_Vectors_Utilities;        use Integer32_Vectors_Utilities;
with Lists_of_Vectors32_Utilities;       use Lists_of_Vectors32_Utilities;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Standard_Integer32_Vertices;        use Standard_Integer32_Vertices;

package body Inner_Normal_Cones is

-- NOTE ON THE IMPLEMENTATION :
--   There is a safety mode in which after each computation of the generators,
--   it is automatically tested whether the generators satisfy the 
--   inequalities of the normal cone.  If not, then the bug is reported and
--   a program_error is raised.

-- AUXILIAIRIES :

  function Lower ( a1,b1,a2,b2 : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if a1/b1 < a2/b2, false otherwise. 

  -- REQUIRED : b1 /= 0 and b2 /= 0.

  begin
    if b1*b2 > 0
     then return (a1*b2 - a2*b1 < 0);
     else return (a1*b2 - a2*b1 > 0);
    end if;
  end Lower;

  function Higher ( a1,b1,a2,b2 : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if a1/b1 > a2/b2, false otherwise. 

  -- REQUIRED : b1 /= 0 and b2 /= 0.

  begin
    if b1*b2 > 0
     then return (a1*b2 - a2*b1 > 0);
     else return (a1*b2 - a2*b1 < 0);
    end if;
  end Higher;

  function Compute_Facets ( L : List; x : Vector ) return Faces is

  -- DESCRIPTION :
  --   Returns all facets of conv(l) that all contain x.
  --   Checks whether x belongs to the list or not.

    res : Faces;
    n : constant integer32 := x'length;

  begin
    if Is_In(L,x) then
      res := Create(n-1,n,L,x);
    else
      declare
        wrk : List := L;
        lx : Link_to_Vector;
      begin
        lx := new Vector'(x);
        Construct(lx,wrk);
        res := Create(n-1,n,wrk,x);
      end;
    end if;
    return res;
  end Compute_Facets;

  procedure Shifted_Points ( L : in List; x : in Vector;
                             res,res_last : in out List ) is

  -- DESCRIPTION :
  --   Appends the list of shifted points w.r.t. x to the list res.

    tmp : List := L;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      Append_Diff(res,res_last,Head_Of(tmp).all-x);
      tmp := Tail_Of(tmp);
    end loop;
  end Shifted_Points;

  function In_Hull ( L : List; x : Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the vector x belongs the convex hull spanned by
  --   the points in l minus the point x itself.

    res : boolean;
    lx : List;

  begin
    Copy(L,lx);
    Remove(lx,x);
    res := Is_In_Hull(x,lx);
    Clear(lx);
    return res;
  end In_Hull;

-- AUXILIAIRIES TO CONSTRUCT THE NORMALS :

  procedure Normal ( m : in out Matrix; v : out Vector ) is

  -- DESCRIPTION :
  --   Computes the normal to the vectors stored in the rows of the matrix m.

  -- REQUIRED : the matrix m is square!

    res : Vector(m'range(2));

  begin
    Upper_Triangulate(m); Scale(m);
    res := (res'range => 0);
    Solve0(m,res);
    v := res;
  end Normal;

  procedure Orientate_Normal ( L : in List; f : in Face; 
                               normal : in out Vector ) is

  -- DESCRIPTION :
  --   Orientates the normal of the face such that it becomes an inner normal
  --   w.r.t. the points in the list l.

    tmp : List := L;
    ip : constant integer32 := f(f'first).all*normal;
    pt : Vector(Head_Of(L)'range);
    done : boolean := false;
    ptip : integer32;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp).all;
      if not Is_In(f,pt) then
        ptip := pt*normal;
        if ptip /= ip then
          if ptip < ip
           then normal := -normal;
          end if;
          done := true;
        end if;
      end if;
      exit when done;
      tmp := Tail_Of(tmp);
    end loop;
  end Orientate_Normal;

  function Inner_Normal ( L : List; facet : Face ) return Vector is

  -- DESCRIPTION :
  --   Returns the inner normal to the facet.

    res,fst : Vector(facet(facet'first)'range);
    m : Matrix(facet'first+1..facet'last,facet(facet'first)'range);

  begin
    fst := facet(facet'first).all;
    for i in m'first(1)..m'last(1) loop
      for j in m'range(2) loop
        m(i,j) := facet(i)(j) - fst(j);
      end loop;
    end loop;
    Normal(m,res);
    Orientate_Normal(l,facet,res);
    return res;
  end Inner_Normal;

  procedure Inner_Normals ( L : in List; facets : in Faces;
                            iv,iv_last : in out List ) is

  -- DESCRIPTION :
  --   Returns the list of all inner normals to the facets of conv(l).
  --   If there is only one facet, then both `inner' and `outer' normal
  --   will be returned, because there is nothing to orientate with.
  --   This means that also the negation of the normal will be in the list iv.

    tmp : Faces := facets;

  begin
    while not Is_Null(tmp) loop
      declare
        v : constant Vector := Inner_Normal(L,Head_Of(tmp));
      begin
        Append(iv,iv_last,v);
        if Length_Of(facets) = 1
         then Append(iv,iv_last,-v);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Inner_Normals;

  function Number_of_Rows ( L : List ) return integer32 is

  -- DESCRIPTION :
  --   Returns Maximum(dimension,Length_Of(l)).

  -- REQUIRED : not Is_Null(l).

    n : constant integer32 := Head_Of(l)'length;
    m : constant integer32 := integer32(Length_Of(l));

  begin
    if n > m
     then return n;
     else return m;
    end if;
  end Number_of_Rows;

  function Normal ( L : List ) return Vector is

  -- DESCRIPTION :
  --   Return a normal to the vectors in the list.

  -- REQUIRED : Length_Of(l) <= n = Head_Of(l)'length
  --   or the points in l all lie on the same facet.

    fst : constant Link_to_Vector := Head_Of(l);
    m : Matrix(1..Number_of_Rows(l),fst'range);
    cnt : integer32 := 0;
    res : Vector(fst'range);
    tmp : List := Tail_Of(L);

  begin
    while not Is_Null(tmp) loop
      res := Head_Of(tmp).all;
      cnt := cnt+1;
      for i in res'range loop
        m(cnt,i) := res(i) - fst(i);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    for i in cnt+1..m'last(1) loop
      for j in m'range(2) loop
        m(i,j) := 0;
      end loop;
    end loop;
    Normal(m,res);
    return res;
  end Normal;

  procedure Normals ( L : in List; x : in Vector; iv,iv_last : in out List ) is

  -- DESCRIPTION :
  --   Computes a list of all normals to the vectors in the list.

  -- REQUIRED : Length_Of(l) <= n = x'length.
  --   or all points in l lie on the same facet.

    nor : constant Vector := Normal(l);
    piv : constant integer32 := Pivot(nor);
    pivnor : Vector(nor'range);   -- normal with pivnor(piv) > 0

  begin
    if piv <= nor'last then
      if x'length > 2 then
        if nor(piv) < 0
         then pivnor := -nor;
         else pivnor := nor;
        end if;
        declare
          t : Transfo := Build_Transfo(pivnor,piv);
          tx : constant Vector := Reduce(t*x,piv);
          tt : Transfo := Transpose(t);
          tl : List := Transform_and_Reduce(t,piv,l);
        begin
          if Length_Of(tl) <= nor'length-1 then
            Normals(tl,tx,iv,iv_last);
          else
            declare
              facets : constant Faces := Compute_Facets(tl,tx);
            begin
              if Length_Of(facets) > 1
               then Inner_Normals(tl,facets,iv,iv_last);
               else Normals(tl,tx,iv,iv_last);
              end if;
            end;
          end if;
          Insert_and_Transform(iv,piv,0,tt);
          iv_last := Pointer_to_Last(iv);
          Clear(t); Deep_Clear(tl); Clear(tt);
        end;
      end if;
      Append(iv,iv_last,nor); Append(iv,iv_last,-nor);
    end if;
  end Normals;

-- CONSTRUCTORS FOR PRIMAL REPRESENTATION :

  function Generators ( L : List; facets : Faces; x : Vector ) return List is

    res,res_last : List;

  begin
    if Length_Of(facets) > 1 then
      Inner_Normals(l,facets,res,res_last);     -- compute inner normals
    elsif In_Hull(L,x) then
      return res;                               -- empty normal cone
    else
      Normals(l,x,res,res_last);                -- compute normals
      if Length_Of(l) = 2
       then Shifted_Points(l,x,res,res_last);
      end if;
    end if;
   -- SAFETY MODE :
   -- if not Contained_in_Cone(l,x,res)
   --  then put_line("Bug raised for computing generators for list"); put(l);
   --       put(" and the point : "); put(x); new_line;
   --       raise PROGRAM_ERROR;
   -- end if;
    return res;
  end Generators;

  function Generators ( L : List; x : Vector ) return List is

    res,res_last : List;
    n : constant integer32 := x'length;
    facets : Faces;

  begin
    if In_Hull(L,x) then
      return res;                                -- empty normal cone
    else
      if Length_Of(L) > natural32(n)
       then facets := Compute_Facets(L,x);
      end if;
      if Length_Of(facets) > 1 then
        res := Generators(L,facets,x);
      else
        Normals(L,x,res,res_last);
        if Length_Of(L) = 2
         then Shifted_Points(L,x,res,res_last);
        end if;
      end if;
      -- SAFETY MODE :
      -- if not Contained_in_Cone(l,x,res)
      --  then put_line("Bug raised for computing generators for list");
      --       put(l); put(" and the point : "); put(x); new_line;
      --       raise PROGRAM_ERROR;
      -- end if;
      return res;
    end if;
  end Generators;

-- CONSTRUCTOR FOR DUAL REPRESENTATION :

  function Inner_Normal_Cone ( L : List; x : Vector ) return Matrix is

    len : constant integer32 := integer32(Length_Of(L));
    res : Matrix(x'range,1..len-1); -- CAUTION : Is_In(l,x) must hold !!
    y : Vector(x'range);
    cnt : integer32 := 0;
    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      y := Head_Of(tmp).all;
      if not Equal(x,y) then                  -- avoid trivial inequality
        cnt := cnt+1;
        for i in x'range loop
          res(i,cnt) := y(i) - x(i);
        end loop;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Inner_Normal_Cone;

  function Included_Normal_Cone ( gx : List; x : Vector ) return Matrix is

    len : constant integer32 := integer32(Length_Of(gx));
    res : Matrix(x'first-1..x'last,1..len);
    tmp : List := gx;
    v : Vector(x'range);

  begin
    for j in res'range(2) loop
      v := Head_Of(tmp).all;
      res(res'first(1),j) := x*v;
      for i in res'first(1)+1..res'last(1) loop
        res(i,j) := v(i);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Included_Normal_Cone;

-- PRIMITIVE SELECTORS :

  function Evaluate ( m : Matrix; i : integer32;
                      v : Vector ) return integer32 is

    sum : integer32 := 0;

  begin
    for j in v'range loop
      sum := sum + m(j,i)*v(j);
    end loop;
    return sum;
  end Evaluate;

  function Satisfies
             ( m : Matrix; i : integer32; v : Vector ) return boolean is
  begin
    if m'first(1) = v'first-1
     then return (Evaluate(m,i,v) >= m(0,i));
     else return (Evaluate(m,i,v) >= 0);
    end if;
  end Satisfies;

  function Strictly_Satisfies
             ( m : Matrix; i : integer32; v : Vector ) return boolean is
  begin
    if m'first(1) = v'first-1
     then return (Evaluate(m,i,v) > m(0,i));
     else return (Evaluate(m,i,v) > 0);
    end if;
  end Strictly_Satisfies;

  function Satisfies ( m : Matrix; v : Vector ) return boolean is
  begin
    for i in m'range(2) loop
      if not Satisfies(m,i,v)
       then return false;
      end if;
    end loop;
    return true;
  end Satisfies;

  function Strictly_Satisfies ( m : Matrix; v : Vector ) return boolean is
  begin
    for i in m'range(2) loop
      if not Strictly_Satisfies(m,i,v)
       then return false;
      end if;
    end loop;
    return true;
  end Strictly_Satisfies;

  function Satisfies ( m : Matrix; v : List ) return boolean is

    tmp : List := v;

  begin
    while not Is_Null(tmp) loop
      if not Satisfies(m,Head_Of(tmp).all)
       then return false;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return true;
  end Satisfies;

  function Strictly_Satisfies ( m : Matrix; v : List ) return boolean is

    tmp : List := v;

  begin
    while not Is_Null(tmp) loop
      if not Strictly_Satisfies(m,Head_Of(tmp).all)
       then return false;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return true;
  end Strictly_Satisfies;

-- SECONDARY SELECTORS :

  function Satisfy ( m : Matrix; L : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      if Satisfies(m,Head_Of(tmp).all)
       then Append(res,res_last,Head_Of(tmp).all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Satisfy;

  function Strictly_Satisfy ( m : Matrix; L : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := L;
    while not Is_Null(tmp) loop
      if Strictly_Satisfies(m,Head_Of(tmp).all)
       then Append(res,res_last,Head_Of(tmp).all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Strictly_Satisfy;

  function Contained_in_Cone ( L : List; x : Vector; v : List )
                             return boolean is

    ineq : constant Matrix(x'range,1..integer32(Length_Of(L))-1)
         := Inner_Normal_Cone(L,x);

  begin
    return Satisfies(ineq,v);
  end Contained_in_Cone;

  function Strictly_Contained_in_Cone ( L : List; x : Vector; v : List )
                                      return boolean is

    ineq : constant Matrix(x'range,1..integer32(Length_Of(L))-1) 
         := Inner_Normal_Cone(L,x);

  begin
    return Strictly_Satisfies(ineq,v);
  end Strictly_Contained_in_Cone;

  function Contained_in_Cone ( L,v : List ) return boolean is

    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      if Contained_in_Cone(L,Head_Of(tmp).all,v)
       then -- put("true for vector : "); put(Head_Of(tmp).all); new_line;
            return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Contained_in_Cone;

  function Strictly_Contained_in_Cone ( L,v : List ) return boolean is

    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      if Strictly_Contained_in_Cone(L,Head_Of(tmp).all,v)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Strictly_Contained_in_Cone;

-- CONCERNING THE UNION OF CONES :

  function In_Union ( v1,v2 : Vector; k1,k2 : Matrix ) return boolean is

  -- ALGORITHM : consider x(t) = (1-t)*v1 + t*v2, for t in [0,1]
  --   and compute t1: smallest t such that Satisfies(k1,x(t))
  --   and t2: largest t such that Satisfies(k2,x(t)).
  --   Then t1 >= t2 makes the test returning true, false otherwise.

    diff : constant Vector := v1-v2;
    t1a,t1b,t2a,t2b,a,b : integer32;        -- t1 = t1a/t1b and t2 = t2a/t2b

  begin
    t1a := 1; t1b := 1;                            -- initially t1 = 1 = 1/1
    for i in k1'range(2) loop
      b := Evaluate(k1,i,diff);
      if b /= 0 then
        a := Evaluate(k1,i,v1);
        if Lower(a,b,t1a,t1b)
         then t1a := a; t1b := b;
        end if;
      end if;
    end loop;
    t2a := 0; t2b := 1;                            -- initially t2 = 0 = 0/1
    for i in k2'range(2) loop
      b := Evaluate(k2,i,diff);
      if b /= 0 then
        a := Evaluate(k2,i,v1);
        if Higher(a,b,t2a,t2b)
         then t2a := a; t2b := b;
        end if;
      end if;
    end loop;
   -- put("In_Union for "); put(v1); put(" and "); put(v2);
   -- put(" : t1 = "); put(t1a,1); put("/"); put(t1b,1);
   -- put(" : t2 = "); put(t2a,1); put("/"); put(t2b,1); new_line;
   -- put_line("k1 : "); put(k1); put_line("k2 : "); put(k2);
    return not Lower(t1a,t1b,t2a,t2b);
  end In_Union;

  function In_Union ( v1,v2 : List; k1,k2 : Matrix ) return boolean is

    tmpv1,tmpv2 : List;

  begin
    tmpv1 := v1;
    while not Is_Null(tmpv1) loop
      tmpv2 := v2;
      while not Is_Null(tmpv2) loop
        if not In_Union(Head_Of(tmpv1).all,Head_Of(tmpv2).all,k1,k2)
         then return false;
         else tmpv2 := Tail_Of(tmpv2);
        end if;
      end loop;
      tmpv1 := Tail_Of(tmpv1);
    end loop;
    return true;
  end In_Union;

end Inner_Normal_Cones;
