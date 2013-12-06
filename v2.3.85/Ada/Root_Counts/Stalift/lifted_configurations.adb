with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;

package body Lifted_Configurations is

-- ROUTINES TO COLLECT LIFTING AND SUPPORTS :

  function Index ( L : List; pt : Standard_Floating_Vectors.Vector )
                 return natural is

  -- DESCRIPTION :
  --   Returns 0 when the points does not belong to the list.
  --   Otherwise, returns the position of the point in the list.

  -- NOTE :
  --   Only the range of pt is checked, the vectors in the list may be longer,
  --   but these additional coordinates are ignored in the comparison.

    tmp : List := L;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
    equ : boolean;

  begin
    for i in 1..integer32(Length_Of(L)) loop
      lpt := Head_Of(tmp);
      equ := true;
      for j in pt'range loop
        if lpt(j) /= pt(j)
         then equ := false;
        end if;
        exit when not equ;
      end loop;
      if equ
       then return integer(lpt(lpt'last-1));
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end Index;

  procedure Strip_Lifting
              ( pts : in out Array_of_Lists;
                nbp : in natural32;
                lif : out Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   The points on entry contain as their two last coordinates the lifting
  --   and their index.  The points are stripped from their lifting and the
  --   lifting vector is constructed.

  -- REQUIRED :
  --   The last component of the points in pts is the lifting value,
  --   the next to last component of each point is the index of the point.

    lv : Standard_Floating_Vectors.Vector(1..integer32(nbp));
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop 
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        declare
          ind : constant integer32 := integer32(lpt(lpt'last-1));
          lip : constant double_float := lpt(lpt'last);
          stp : Standard_Floating_Vectors.Vector(lpt'first..lpt'last-1);
        begin
          stp := lpt(stp'range);
          Standard_Floating_Vectors.Clear(lpt);
          lpt := new Standard_Floating_Vectors.Vector'(stp);
          Set_Head(tmp,lpt);
          lv(ind) := lip;
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    lif := new Standard_Floating_Vectors.Vector'(lv);
  end Strip_Lifting;

-- ROUTINES TO CONSTRUCT INCIDENCE MATRIX :

  function Is_In ( L : List; i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the list contains a point with last coordinate
  --   equal to i.

  -- REQUIRED :
  --   The meaning of the last coordinate of each point in the list
  --   must be equal to the index of that point.

    tmp : List := L;
    pt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if integer32(pt(pt'last)) = i
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( mic : Mixed_Cell; i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the cell contains a point with last coordinate
  --   equal to i.

  -- REQUIRED :
  --   The meaning of the last coordinate of each point in the cell
  --   must be equal to the index of that point.

  begin
    for k in mic.pts'range loop
      if Is_In(mic.pts(k),i)
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Last_Entry ( im : Incidence_Matrix; i,k : integer32 ) 
                      return integer32 is

  -- DESCRIPTION :
  --   Returns the last index j >= k for which im(i,j) is true.
  --   If the index on return is zero, then the ith point does 
  --   not occur in the remaining cells k..im'last(2).

    res : integer32 := 0;

  begin
    for j in k..im'last(2) loop
      if im(i,j)
       then res := j;
      end if;
    end loop;
    return res;
  end Last_Entry;

  function Is_In ( v : Standard_Integer_Vectors.Vector;
                   i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there is a j such that v(j) = i.

  begin
    for j in v'range loop
      if v(j) = i
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  procedure Minimal_Last_Entry
              ( im : in Incidence_Matrix; k : in integer32;
                done : in Standard_Integer_Vectors.Vector;
                row,col : out integer32 ) is

  -- DESCRIPTION :
  --   Returns the index to the row with the minimal last entry,
  --   as well as the colum where that last true value occured,
  --   under the restriction that row does not belong to done.

    minrow,mincol,val : integer32;

  begin
    put("Starting column : "); put(k,1);
    put("  Points done : "); put(done); new_line;
    mincol := im'last(2)+1;
    for i in im'range(1) loop
      if not Is_In(done,i) then
        val := Last_Entry(im,i,k);
        if val < mincol
         then mincol := val; minrow := i;
        end if;
      end if;
      exit when (mincol = 0);
    end loop;
    row := minrow; col := mincol;
    put("  minrow : "); put(minrow,1); put("  mincol : "); put(mincol,1);
    new_line;
  end Minimal_Last_Entry;

  function Number_of_Points ( im : Incidence_Matrix; i : integer32 )
                            return natural32 is

  -- DESCRIPTION :
  --   Returns the number of points in the ith cell.

    cnt : natural32 := 0;

  begin
    for j in im'range(1) loop
      if im(j,i)
       then cnt := cnt + 1;
      end if;
    end loop;
    return cnt;
  end Number_of_Points;

-- ROUTINES TO ASSIGN THE LIFTING :

  procedure Assign_Lifting
              ( pt : in out Standard_Floating_Vectors.Link_to_Vector;
                lifting : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns a lifting to the point pt.

  -- REQUIRED :
  --   The last coordinate of pt must correspond to its index.

  begin
    pt(pt'last) := lifting(integer32(pt(pt'last)));
  end Assign_Lifting;

  procedure Assign_Lifting
              ( L : in out List;
                lifting : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns a lifting to each point in the list.

  -- REQUIRED :
  --   The last coordinate of each point must correspond to its index.

    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      declare
        lpt : Standard_Floating_Vectors.Link_to_Vector := Head_Of(tmp);
      begin
        Assign_Lifting(lpt,lifting);
        Set_Head(tmp,lpt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Assign_Lifting;

  procedure Assign_Lifting
              ( pts : in out Array_of_Lists;
                lifting : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns the lifting to points in the lists pts.

  -- REQUIRED :
  --   The last coordinate of each point must correspond to its index.

  begin
    for i in pts'range loop
      Assign_Lifting(pts(i),lifting);
    end loop;
  end Assign_Lifting;

  function Compute_Inner_Normal
              ( n : in integer32; pts : in Array_of_Lists )
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..n+1);
    mat : Standard_Floating_Matrices.Matrix(1..n,1..n);
    rhs : Standard_Floating_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    row : integer32 := 0;
    shi,lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop
      shi := Head_Of(pts(i));
      declare
        tmp : List := Tail_Of(pts(i));
      begin
        while not Is_Null(tmp) loop
          lpt := Head_Of(tmp);
          row := row+1;
          for i in 1..n loop
            mat(row,i) := lpt(i) - shi(i);
          end loop;
          rhs(row) := shi(shi'last) - lpt(lpt'last);
          tmp := Tail_Of(tmp);
        end loop;
      end;
    end loop;
    lufac(mat,n,ipvt,info);
    lusolve(mat,n,ipvt,rhs);
    res(rhs'range) := rhs;
    res(res'last) := 1.0;
    return res;
  end Compute_Inner_Normal;

  procedure Assign_Lifting
              ( n : in integer32;
                mic : in out Floating_Mixed_Subdivisions.Mixed_Cell;
                lifting : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns the lifting to points in the cells of mixsub.

  -- REQUIRED :
  --   The last coordinate of each point must be correspond to its index.

  begin
    for i in mic.pts'range loop
      Assign_Lifting(mic.pts(i),lifting);
    end loop;
  end Assign_Lifting;

-- ROUTINES TO SET UP THE LINEAR-OPTIMIZATION MODEL :

  procedure Collect_Points_and_Lifting
               ( n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 mixsub : in out Mixed_Subdivision;
                 pts : out Array_of_Lists; nbp : out natural32;
                 lifting : out Standard_Floating_Vectors.Link_to_Vector ) is

    res,res_last : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    tmpsub : Mixed_Subdivision := mixsub;
    cnt : natural32 := 0;

    procedure Collect_Data_from_Cell ( mic : in out Mixed_Cell ) is

      tmp : Lists_of_Floating_Vectors.List;

    begin
      for i in mic.pts'range loop
        tmp := mic.pts(i);
        while not Is_Null(tmp) loop
          declare
            lpt : constant Standard_Floating_Vectors.Link_to_Vector
                := Head_Of(tmp);
            ind : constant natural := Index(res(i),lpt(lpt'first..lpt'last-1));
          begin
            if ind /= 0
             then
               lpt(lpt'last) := double_float(ind);
             else
               cnt := cnt+1;
               declare
                 fip : Standard_Floating_Vectors.Vector(lpt'first..lpt'last+1);
               begin
                 fip(lpt'first..lpt'last-1) := lpt(lpt'first..lpt'last-1);
                 fip(fip'last) := lpt(lpt'last);     -- copy lifting
                 lpt(lpt'last) := double_float(cnt); -- replace lifting
                 fip(lpt'last) := lpt(lpt'last);     -- set index
                 Append(res(i),res_last(i),fip);
               end;
            end if;
            Set_Head(tmp,lpt);
            tmp := Tail_Of(tmp);
          end;
        end loop;
      end loop;
    end Collect_Data_from_Cell;

  begin
    while not Is_Null(tmpsub) loop
      declare
        mic : Mixed_Cell := Head_Of(tmpsub);
      begin
        Collect_Data_from_Cell(mic);
        Set_Head(tmpsub,mic);
      end;
      tmpsub := Tail_Of(tmpsub);
    end loop;
    Strip_Lifting(res,cnt,lifting);
    pts := res;
    nbp := cnt;
  end Collect_Points_and_Lifting;

  procedure Get_Point ( pts : in Array_of_Lists; k : in integer32;
                        lpk : out Standard_Floating_Vectors.Link_to_Vector ) is

    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop 
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        if integer32(lpt(lpt'last)) = k
         then lpk := lpt; exit;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Get_Point;

  function Create ( np : integer32; mixsub : Mixed_Subdivision )
                  return Incidence_Matrix is

    res : Incidence_Matrix(1..np,1..integer32(Length_Of(mixsub)));
    tmp : Mixed_Subdivision;

  begin
    for i in 1..np loop
      tmp := mixsub;
      for j in 1..integer32(Length_Of(mixsub)) loop
        res(i,j) := Is_In(Head_Of(tmp),i);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Create;

  procedure Sort ( im : in Incidence_Matrix;
                   order : out Standard_Integer_Vectors.Vector;
                   fail : out boolean ) is

    pts : Standard_Integer_Vectors.Vector(im'range(1)) := (im'range(1) => 0);
    cnt : integer32 := 0;
   -- wrk : Incidence_Matrix(im'range(1),im'range(2)) := im;
    col : integer32 := 1;
    nextrow,nextcol : integer32;

  begin
    for i in pts'range loop
      Minimal_Last_Entry(im,col,pts(1..cnt),nextrow,nextcol);
      if nextcol = 0 then
        fail := true; exit;
      elsif nextcol < im'last(2) then
        cnt := cnt+1;
        pts(cnt) := nextrow;
        col := nextcol+1;
      else -- nextcol = im'last(2)
        put("adding ");
        if natural32(cnt) = natural32(pts'last)
                          - Number_of_Points(im,nextcol) then
          put("points ");
          for j in im'range(1) loop
            if im(j,nextcol) then
              cnt := cnt+1;
              put(" "); put(j,1);
              pts(cnt) := j;
            end if;
          end loop;
          fail := false; exit;
        else
          fail := true; exit;
        end if;
      end if;
      put("Points : "); put(pts(1..cnt)); new_line;
    end loop;
    order := pts;
  end Sort;

  procedure Assign_Lifting
              ( n : in integer32; mixsub : in out Mixed_Subdivision;
                lifting : in Standard_Floating_Vectors.Vector ) is

    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        Assign_Lifting(n,mic,lifting);
        mic.nor.all := Compute_Inner_Normal(n,mic.pts.all);
        Set_Head(tmp,mic);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Assign_Lifting;

end Lifted_Configurations;
