with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;

package body Cayley_Embedding is

-- AUXILIARIES :

  function Is_Good_Point
             ( cnt,n : integer32; pt : Link_to_Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the point pt is a point of the type indicated
  --   by the parameter cnt, i.e. whether it belongs to the polytope
  --   placed on the vertex with number cnt.

    goodpoint : boolean;

  begin
    if cnt = 0 then
      goodpoint := true;
      for L in pt'first..pt'last-n-1 loop
        if pt(L) /= 0
         then goodpoint := false;
        end if;
        exit when not goodpoint;
      end loop;
    else
      goodpoint := (pt(cnt) = 1);
    end if;
    return goodpoint;
  end Is_Good_Point;

  procedure Project ( n : integer32; v : in out Link_to_Vector ) is

  -- DESCRIPTION :
  --   After application, v points to a vector of length n+1.

    newv : Link_to_Vector;

  begin
    newv := new Vector(1..n+1);
    newv(1..n+1) := v(v'last-n..v'last);
    Clear(v);
    v := newv;
  end Project;

-- TARGET ROUTINES :

  function Embedding_Before_Lifting 
               ( supports : Array_of_Lists ) return List is

    tmp,res,res_last : List;
    r1 : constant integer32 := supports'length-1;
    pt : Link_to_Vector;
    cnt : integer32 := 0;

  begin
    for k in supports'range loop
      tmp := supports(k);
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        declare
          npt : Vector(pt'first..pt'last+r1);
        begin
          npt(npt'last-pt'length+1..npt'last) := pt.all;
          npt(npt'first..npt'first+r1-1) := (npt'first..npt'first+r1-1 => 0);
          if cnt > 0
           then npt(cnt) := 1;
          end if;
          Append(res,res_last,npt);
        end;
        tmp := Tail_Of(tmp);
      end loop;
      cnt := cnt + 1;
    end loop;
    return res;
  end Embedding_Before_Lifting;

  function Extract ( vtp,n : integer32; pts : VecVec ) return List is

    res,res_last : List;

  begin
    for k in pts'range loop
      if Is_Good_Point(vtp,n,pts(k))
       then Append(res,res_last,pts(k).all);
      end if;
    end loop;
    return res;
  end Extract;

  function Extract ( vtp,n : integer32; pts : List ) return List is

   -- DESCRIPTION :
   --   Extracts the points out of the list that are of the type
   --   indicated by vtp.

    tmp,res,res_last : List;
    pt : Link_to_Vector;
 
  begin
    tmp := pts;
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if Is_Good_Point(vtp,n,pt)
       then Append(res,res_last,pt.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Extract;

  function Extract_Mixed_Cell 
             ( n : integer32; mix : Vector; s : Simplex ) return Mixed_Cell is

    res : Mixed_Cell;
    work : Array_of_Lists(mix'range);
    cnt : integer32 := 0;
    iscell : boolean;
    pts : constant VecVec := Vertices(s);

  begin
    for k in mix'range loop
      work(k) := Extract(cnt,n,pts);
      iscell := (integer32(Length_Of(work(k))) = mix(k)+1);
      exit when not iscell;
      cnt := cnt + 1;
    end loop;
    if iscell then
      res.pts := new Array_of_Lists'(work);
      res.nor := new vector'(Normal(s));
    else
      Deep_Clear(work);
    end if;
    return res;
  end Extract_Mixed_Cell;

  function Extract_Mixed_Cells
             ( n : integer32; mix : Vector; t : Triangulation ) 
             return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    s : Simplex;
    tmp : Triangulation;
  
  begin
    tmp := t;
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      declare
        mic : Mixed_Cell := Extract_Mixed_Cell(n,mix,s);
      begin
        if mic.nor /= null
         then Append(res,res_last,mic);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Extract_Mixed_Cells;

  function Extract_non_Flat_Mixed_Cells
             ( n : integer32; mix : Vector; t : Triangulation )
             return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    s : Simplex;
    tmp : Triangulation;
 
  begin
    tmp := t;
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      exit when Is_Flat(s);
      declare
        mic : Mixed_Cell := Extract_Mixed_Cell(n,mix,s);
      begin
        if mic.nor /= null
         then Append(res,res_last,mic);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Extract_non_Flat_Mixed_Cells;

  procedure Deflate ( n : integer32; L : in out List ) is

    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      declare
        pt : Link_to_Vector := Head_Of(tmp);
      begin
        Project(n,pt);
        Set_Head(tmp,pt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Deflate;

  procedure Deflate ( n : integer32; mic : in out Mixed_Cell ) is
  begin
    Project(n,mic.nor);
    for k in mic.pts'range loop
      Deflate(n,mic.pts(k));
    end loop;
  end Deflate;

  procedure Deflate ( n : integer32; mixsub : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        Deflate(n,mic);
        Set_Head(tmp,mic);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Deflate;

end Cayley_Embedding;
