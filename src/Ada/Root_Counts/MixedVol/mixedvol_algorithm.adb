with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Supports_of_Polynomial_Systems;
with Prepare_for_MV,Mixed_Volume;       use Prepare_for_MV,Mixed_Volume;

-- for testing :
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;

package body MixedVol_Algorithm is

  function Flatten ( n,m : integer32; v : Standard_Integer_VecVecs.VecVec )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n*m);
    cnt : integer32 := 0;

  begin
    for i in v'range loop
      for j in v(i)'range loop
        cnt := cnt + 1;
        res(cnt) := v(i)(j);
      end loop;
    end loop;
    return res;
  end Flatten;

  procedure Flatten_Supports
             ( s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               ind : out Standard_Integer_Vectors.Vector;
               sup : out Standard_Integer_VecVecs.VecVec ) is

    use Lists_of_Integer_Vectors;

    ptr : List;
    cnt : integer32 := 0;

  begin
    for i in s'range loop
      ptr := s(i); 
      ind(i) := cnt + 1;
      while not Is_Null(ptr) loop
        cnt := cnt + 1;
        sup(cnt) := Head_Of(ptr);
        ptr := Tail_Of(ptr);
      end loop;
    end loop;
  end Flatten_Supports;

  procedure Extract_Supports
             ( n : in integer32; p : in Poly_Sys; m : out integer32;
               ind,cnt : out Standard_Integer_Vectors.Vector;
               sup : out Standard_Integer_Vectors.Link_to_Vector ) is

    s : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);

  begin
    Extract_Supports(n,s,m,ind,cnt,sup);
  end Extract_Supports;

  procedure Extract_Supports
             ( n : in integer32; p : in Laur_Sys; m : out integer32;
               ind,cnt : out Standard_Integer_Vectors.Vector;
               sup : out Standard_Integer_Vectors.Link_to_Vector ) is

    s : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);

  begin
    Extract_Supports(n,s,m,ind,cnt,sup);
  end Extract_Supports;

  procedure Extract_Supports
             ( n : in integer32; 
               s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
	       m : out integer32;
               ind,cnt : out Standard_Integer_Vectors.Vector;
               sup : out Standard_Integer_Vectors.Link_to_Vector ) is
  begin
   -- new_line;
   -- put_line("The supports of your system :"); put(s);
    m := 0;
    for i in s'range loop
      cnt(i) := integer32(Lists_of_Integer_Vectors.Length_Of(s(i)));
      m := m + cnt(i);
    end loop;
   -- put("Cardinality of supports :"); put(cnt); new_line;
   -- put(" total number of points : "); put(m,1); put_line(".");
    declare
      spv : Standard_Integer_VecVecs.VecVec(1..m);
    begin
     -- Flatten_Supports(n,m,s,ind,spv);
      Flatten_Supports(s,ind,spv);
     -- put_line("Flattened supports :"); put(spv);
     -- put("Index set :"); put(ind); new_line;
      sup := new Standard_Integer_Vectors.Vector'(Flatten(n,m,spv));
    end;
  end Extract_Supports;

  procedure Write_Supports 
                ( nSpt : in integer32;
                  Idx : in Standard_Integer_Vectors.Link_to_Vector;
                  Spt : in Standard_Integer_VecVecs.Link_to_VecVec ) is
  begin
    put("The index set : "); put(Idx.all); new_line;
    for i in 1..nSpt loop
      put("Support "); put(i,1); put_line(" :");
      for j in Idx(i-1)..Idx(i)-1 loop
        put(Spt(j)); new_line;
      end loop;
    end loop;
  end Write_Supports;

  function Is_In ( k : integer32; v : Standard_Integer_Vectors.Vector;
                   Idx : Standard_Integer_Vectors.Link_to_Vector;
                   Spt : Standard_Integer_VecVecs.Link_to_VecVec )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if the vector v belongs to the kth support.

  begin
    for i in Idx(k-1)..Idx(k)-1 loop
      if Standard_Integer_Vectors.Equal(Spt(i).all,v)
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Number_of_Occurrences
             ( nSpt : integer32;
               Idx : Standard_Integer_Vectors.Link_to_Vector;
               Spt : Standard_Integer_VecVecs.Link_to_VecVec;
               v : Standard_Integer_Vectors.Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of times the vector belongs to the supports.

    res : integer32 := 0;

  begin
    for i in 1..nSpt loop
      if Is_In(i,v,Idx,Spt)
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Occurrences;

  procedure Copy_Support
              ( k : in integer32;
                Idx : in Standard_Integer_Vectors.Link_to_Vector;
                Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                newIdx : in Standard_Integer_Vectors.Link_to_Vector;
                newSpt : in Standard_Integer_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Copies k-th support from Spt to newSpt.

    ind : integer32 := newIdx(k-1);

  begin
    for i in Idx(k-1)..Idx(k)-1 loop
      newSpt(ind) := new Standard_Integer_Vectors.Vector'(Spt(i).all);
      ind := ind + 1;
    end loop;
    newIdx(k) := ind;
  end Copy_Support;

  procedure Add_Point
              ( k : in integer32; v : in Standard_Integer_Vectors.Vector;
                Idx : in Standard_Integer_Vectors.Link_to_Vector;
                Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                newIdx : in Standard_Integer_Vectors.Link_to_Vector;
                newSpt : in Standard_Integer_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Copies k-th support from Spt to newSpt, adding the point v
  --   to the front of the new k-th support.

    ind : integer32 := newIdx(k-1);

  begin
    newSpt(ind) := new Standard_Integer_Vectors.Vector'(v);
    ind := ind + 1;
    for i in Idx(k-1)..Idx(k)-1 loop
      newSpt(ind) := new Standard_Integer_Vectors.Vector'(Spt(i).all);
      ind := ind + 1;
    end loop;
    newIdx(k) := ind;
  end Add_Point;
 
  procedure Add_Artificial_Origins
                ( nVar,nSpt : in integer32;
                  Idx : in out Standard_Integer_Vectors.Link_to_Vector;
                  Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  nbadd : out integer32;
                  added : out Standard_Integer_Vectors.Link_to_Vector ) is

    origin : constant Standard_Integer_Vectors.Vector(0..nVar-1)
           := (0..nVar-1 => 0);
    newIdx : Standard_Integer_Vectors.Link_to_Vector;
    newSpt : Standard_Integer_VecVecs.Link_to_VecVec;
    addind : integer32 := 0;

  begin
    nbadd := nSpt - Number_of_Occurrences(nSpt,Idx,Spt,origin);
    if nbadd > 0 then
      newIdx := new Standard_Integer_Vectors.Vector'(Idx.all);
      added := new Standard_Integer_Vectors.Vector(1..nbadd);
      newSpt := new Standard_Integer_VecVecs.VecVec(0..Spt'last+nbadd);
      for i in 1..nSpt loop
        if Is_In(i,origin,Idx,Spt) then
          Copy_Support(i,Idx,Spt,newIdx,newSpt);
        else
          Add_Point(i,origin,Idx,Spt,newIdx,newSpt);
          addind := addind + 1;
          added(addind) := newIdx(i-1);
        end if;
      end loop;
      Standard_Integer_Vectors.Clear(Idx); Idx := newIdx;
      Standard_Integer_VecVecs.Deep_Clear(Spt); Spt := newSpt;
    end if;
  end Add_Artificial_Origins;

--  function max ( x,y : integer ) return integer is
--
--  -- DESCRIPTION :
--  --   Returns the maximum of x and y.
--
--  begin
--    if x >= y
--     then return x;
--     else return y;
--    end if;
--  end max;

--  function min ( x,y : integer ) return integer is
--
--  -- DESCRIPTION :
--  --   Returns the minimum of x and y.
--
--  begin
--    if x <= y
--     then return x;
--     else return y;
--    end if;
--  end min;

  function quick_return
            -- ( nVar : natural; 
            --   SptIdx : Standard_Integer_Vectors.Link_to_Vector;
            --   Spt : Standard_Integer_VecVecs.Link_to_VecVec )
             return integer32 is

  -- DESCRITION :
  --   Returns 1 if the system is univariate, or linear, or if some
  --   supports have fewer than two terms in them, printing a message;
  --   otherwise 0 is returned.

  begin
    return 0;
  end quick_return;

  procedure mv ( nVar,nPts : in integer32;
                 ind,cnt,sup : in Standard_Integer_Vectors.Vector;
                 stlb : in double_float; nSpt : out integer32;
                 SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : out Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nbCells : out integer32; cells : out CellStack;
                 mixvol : out natural32;
                 multprec_hermite : in boolean := false ) is

    SptIdx : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..nVar);
    Spt : Standard_Integer_VecVecs.Link_to_VecVec
        := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
    NuIdx2OldIdx,added : Standard_Integer_Vectors.Link_to_Vector;
    idx,nbadd : integer32;
    MVol : natural32;

  begin
   -- put("The vector cnt : "); put(cnt); new_line;
    for i in 0..(nVar-1) loop
      SptIdx(i) := cnt(i+1);
    end loop;
    SptIdx(nVar) := nPts;
    for i in reverse 0..nVar-1 loop
      SptIdx(i) := SptIdx(i+1) - SptIdx(i);
    end loop;
    idx := 0;
    for i in Spt'range loop
      Spt(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
      for j in 0..(nVar-1) loop
        idx := idx + 1;
        Spt(i)(j) := sup(idx);
      end loop;
    end loop;
    if quick_return = 1 -- quick_return(nVar,SptIdx,Spt) = 1
     then mixvol := 0; nSpt := nPts; CellSize := 0; nbCells := 0; return;
    end if;
    SptType := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    VtxIdx := new Standard_Integer_Vectors.Vector(0..nVar);
    Vtx := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
    for i in 0..(nPts-1) loop
      Vtx(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    end loop;
    NuIdx2OldIdx := new Standard_Integer_Vectors.Vector(0..(nPts-1));
    nSpt := nVar;
   -- Write_Supports(nSpt,SptIdx,Spt);
    Pre4MV(nVar,nSpt,nSpt,SptType,Spt,SptIdx,Vtx,VtxIdx,NuIdx2OldIdx,perm);
    if stlb /= 0.0 then
     -- put_line("Supports before adding artificial origins :");
     -- Write_Supports(nSpt,VtxIdx,Vtx);
      Add_Artificial_Origins(nVar,nSpt,VtxIdx,Vtx,nbadd,added);
     -- put_line("Supports after adding artificial origins :");
     -- Write_Supports(nSpt,VtxIdx,Vtx);
    else
      nbadd := 0;
    end if;
   -- Write_Supports(nSpt,VtxIdx,Vtx);
    lft := new Standard_Floating_Vectors.Vector(0..(VtxIdx(nSpt)-1));
    if nbadd = 0 then
      for i in 0..(VtxIdx(nSpt)-1) loop
        lft(i) := 2.0*(1.5 + Standard_Random_Numbers.Random);
      end loop;
    else
      for i in 0..(VtxIdx(nSpt)-1) loop
        lft(i) := Standard_Random_Numbers.Random;
      end loop;
      for i in added'range loop
        lft(added(i)) := stlb;
      end loop;
    end if;
   -- put_line("The lifting values : "); put_line(lft);
    CellSize := cell_size(nSpt,SptType);
    Cs_Init(cells,CellSize);
    MixedVol(nVar,nSpt,CellSize,SptType,VtxIdx,Vtx,lft,nbCells,cells,MVol,
             multprec_hermite);
    Standard_Integer_Vectors.Clear(SptIdx);
    Standard_Integer_Vectors.Clear(NuIdx2OldIdx);
    Standard_Integer_VecVecs.Deep_Clear(Spt);
    mixvol := MVol;
  end mv;

  procedure mv_upto_pre4mv
             ( nVar,nPts : in integer32;
               ind,cnt,sup : in Standard_Integer_Vectors.Vector;
               nSpt : out integer32;
               SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
               VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
               Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
               SptIdx : out Standard_integer_Vectors.Link_to_Vector;
               Spt : out Standard_Integer_VecVecs.Link_to_VecVec; 
               NuIdx2OldIdx : out Standard_Integer_Vectors.Link_to_Vector ) is

    idx : integer32;

  begin
    SptIdx := new Standard_Integer_Vectors.Vector(0..nVar);
    Spt := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
   -- put("The vector cnt : "); put(cnt); new_line;
    for i in 0..(nVar-1) loop
      SptIdx(i) := cnt(i+1);
    end loop;
    SptIdx(nVar) := nPts;
    for i in reverse 0..nVar-1 loop
      SptIdx(i) := SptIdx(i+1) - SptIdx(i);
    end loop;
    idx := 0;
    for i in Spt'range loop
      Spt(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
      for j in 0..(nVar-1) loop
        idx := idx + 1;
        Spt(i)(j) := sup(idx);
      end loop;
    end loop;
   -- if quick_return = 1 -- quick_return(nVar,SptIdx,Spt) = 1
   --  then nSpt := 0; return; -- no quick return here
   -- end if;
    SptType := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    VtxIdx := new Standard_Integer_Vectors.Vector(0..nVar);
    Vtx := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
    for i in 0..(nPts-1) loop
      Vtx(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    end loop;
    NuIdx2OldIdx := new Standard_Integer_Vectors.Vector(0..(nPts-1));
    nSpt := nVar;
    Pre4MV(nVar,nSpt,nSpt,SptType,Spt,SptIdx,Vtx,VtxIdx,NuIdx2OldIdx,perm);
  end mv_upto_pre4mv;

  procedure mv_lift
              ( nVar : in integer32;
                stlb : in double_float; nSpt : in integer32;
                VtxIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                Vtx : in out Standard_Integer_VecVecs.Link_to_VecVec;
                lft : out Standard_Floating_Vectors.Link_to_Vector ) is

    nbadd : integer32 := 0;
    added : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if stlb /= 0.0 then
     -- put_line("Supports before adding artificial origins :");
     -- Write_Supports(nSpt,VtxIdx,Vtx);
      Add_Artificial_Origins(nVar,nSpt,VtxIdx,Vtx,nbadd,added);
     -- put_line("Supports after adding artificial origins :");
     -- Write_Supports(nSpt,VtxIdx,Vtx);
    end if;
   -- Write_Supports(nSpt,VtxIdx,Vtx);
    lft := new Standard_Floating_Vectors.Vector(0..(VtxIdx(nSpt)-1));
    if nbadd = 0 then
      for i in 0..(VtxIdx(nSpt)-1) loop
        lft(i) := 2.0*(1.5 + Standard_Random_Numbers.Random);
      end loop;
    else
      for i in 0..(VtxIdx(nSpt)-1) loop
        lft(i) := Standard_Random_Numbers.Random;
      end loop;
      for i in added'range loop
        lft(added(i)) := stlb;
      end loop;
    end if;
   -- put_line("The lifting values : "); put_line(lft);
   Standard_Integer_Vectors.Clear(added);
  end mv_lift;

  procedure mv_with_callback
               ( nVar,nPts : in integer32;
                 ind,cnt,sup : in Standard_Integer_Vectors.Vector;
                 stlb : in double_float; nSpt : out integer32;
                 SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : out Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nbCells : out integer32; cells : out CellStack;
                 mixvol : out natural32;
                 multprec_hermite : in boolean := false;
                 next_cell : access procedure
                   ( idx : Standard_Integer_Vectors.Link_to_Vector )
                   := null ) is

    SptIdx : Standard_Integer_Vectors.Link_to_Vector;
    Spt : Standard_Integer_VecVecs.Link_to_VecVec;
    NuIdx2OldIdx : Standard_Integer_Vectors.Link_to_Vector;
    MVol : natural32;

  begin
    mv_upto_pre4mv
      (nVar,nPts,ind,cnt,sup,nSpt,SptType,perm,VtxIdx,Vtx,
       SptIdx,Spt,NuIDX2OldIdx);
    mv_lift(nVar,stlb,nSpt,VtxIdx,Vtx,lft);
    CellSize := cell_size(nSpt,SptType);
    Cs_Init(cells,CellSize);
    MixedVol_with_Callback
      (nVar,nSpt,CellSize,SptType,VtxIdx,Vtx,lft,nbCells,cells,MVol,
       multprec_hermite,next_cell);
    Standard_Integer_Vectors.Clear(SptIdx);
    Standard_Integer_Vectors.Clear(NuIdx2OldIdx);
    Standard_Integer_VecVecs.Deep_Clear(Spt);
    mixvol := MVol;
  end mv_with_callback;

  procedure uliftmv 
               ( nVar,nPts : in integer32;
                 ind,cnt,sup : in Standard_Integer_Vectors.Vector;
                 stlb : in double_float; nSpt : out integer32;
                 SptType,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 VtxIdx : out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : out Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : out Standard_Floating_Vectors.Link_to_Vector;
                 CellSize,nbCells : out integer32; cells : out CellStack;
                 mixvol : out natural32 ) is

    SptIdx : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..nVar);
    Spt : Standard_Integer_VecVecs.Link_to_VecVec
        := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
    NuIdx2OldIdx,added : Standard_Integer_Vectors.Link_to_Vector;
    idx,nbadd,tmpind : integer32;
    MVol : natural32;

  begin
   -- put("The vector cnt : "); put(cnt); new_line;
    for i in 0..(nVar-1) loop
      SptIdx(i) := cnt(i+1);
    end loop;
    SptIdx(nVar) := nPts;
    for i in reverse 0..nVar-1 loop
      SptIdx(i) := SptIdx(i+1) - SptIdx(i);
    end loop;
    idx := 0;
    for i in Spt'range loop
      Spt(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
      for j in 0..(nVar-1) loop
        idx := idx + 1;
        Spt(i)(j) := sup(idx);
      end loop;
    end loop;
    if quick_return = 1 -- quick_return(nVar,SptIdx,Spt) = 1
     then mixvol := 0; nSpt := nPts; CellSize := 0; nbcells := 0; return;
    end if;
    SptType := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    VtxIdx := new Standard_Integer_Vectors.Vector(0..nVar);
    Vtx := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
    for i in 0..(nPts-1) loop
      Vtx(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    end loop;
    NuIdx2OldIdx := new Standard_Integer_Vectors.Vector(0..(nPts-1));
    nSpt := nVar;
   -- Write_Supports(nSpt,SptIdx,Spt);
    Pre4MV(nVar,nSpt,nSpt,SptType,Spt,SptIdx,Vtx,VtxIdx,NuIdx2OldIdx,perm);
    if stlb /= 0.0 then
     -- put_line("Supports before adding artificial origins :");
     -- Write_Supports(nSpt,VtxIdx,Vtx);
      Add_Artificial_Origins(nVar,nSpt,VtxIdx,Vtx,nbadd,added);
     -- put_line("Supports after adding artificial origins :");
     -- Write_Supports(nSpt,VtxIdx,Vtx);
    else
      nbadd := 0;
    end if;
   -- Write_Supports(nSpt,VtxIdx,Vtx);
    lft := new Standard_Floating_Vectors.Vector(0..(VtxIdx(nSpt)-1));
    if nbadd = 0 then
      for i in 0..(VtxIdx(nSpt)-1) loop
       -- lft(i) := 2.0*(1.5 + Standard_Random_Numbers.Random);
       -- put("Give lifting for "); put(i,1);
       -- put(" with coordinates "); put(vtx(i)); put(" : ");
       -- get(lft(i));
        put("lifting for "); put(vtx(i)); put(" : ");
        tmpind := vtx(i)'first;
        lft(i) := double_float(vtx(i)(tmpind));
        vtx(i)(tmpind) := Standard_Random_Numbers.Random(0,9);
        put(vtx(i)); put("|"); put(lft(i)); new_line;
      end loop;
    else
      for i in 0..(VtxIdx(nSpt)-1) loop
       -- lft(i) := Standard_Random_Numbers.Random;
        put("Give lifting value for "); put(i,1); put(" : ");
        get(lft(i));
      end loop;
      for i in added'range loop
        lft(added(i)) := stlb;
      end loop;
    end if;
   -- put_line("The lifting values : "); put_line(lft);
    CellSize := cell_size(nSpt,SptType);
    Cs_Init(cells,CellSize);
    MixedVol(nVar,nSpt,CellSize,SptType,VtxIdx,Vtx,lft,nbCells,cells,MVol);
    Standard_Integer_Vectors.Clear(SptIdx);
    Standard_Integer_Vectors.Clear(NuIdx2OldIdx);
    Standard_Integer_VecVecs.Deep_Clear(Spt);
    mixvol := MVol;
  end uliftmv;

  function Supports_of_Mixed_Cell
               ( nVar,nSpt : integer32;
                 SptType,labels :  Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : Standard_Floating_Vectors.Link_to_Vector )
               return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    use Lists_of_Floating_Vectors;

    res,last : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..nSpt);
    ind : integer32 := -1;
    point : Standard_Floating_Vectors.Vector(1..nVar+1);

  begin
    for i in 0..(nSpt-1) loop
      for j in 0..SptType(i) loop
        ind := ind + 1;
        for k in 0..(nVar-1) loop
          point(k+1) := double_float(Vtx(labels(ind))(k));
        end loop;
        point(nVar+1) := lft(labels(ind));
        Append(res(i+1),last(i+1),point);
      end loop;
    end loop;
    return res;
  end Supports_of_Mixed_Cell;

  function Supports_of_Mixed_Cell
               ( nVar,nSpt : integer32;
                 SptType,perm :  Standard_Integer_Vectors.Link_to_Vector;
                 labels :  Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : Standard_Floating_Vectors.Link_to_Vector )
               return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

  -- NOTE : If the system is semi-mixed, then we have to be careful
  --   in selecting the points of the supports, see the use of offset.

    use Lists_of_Floating_Vectors;

    res,last : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..nSpt);
    ind : integer32 := -1;
    point : Standard_Floating_Vectors.Vector(1..nVar+1);
    offset : integer32 := 0;
    first : boolean;

  begin
   -- put_line("in Supports_of_Mixed_Cell ...");
    for i in 0..(nSpt-1) loop
      first := true;
      for j in 0..SptType(i) loop
        ind := ind + 1;
        for k in 0..(nVar-1) loop
          point(k+1) := double_float(Vtx(labels(ind))(k));
        end loop;
        point(nVar+1) := lft(labels(ind));
        if first then
         -- put("storing points in list("); put(perm(i+offset)+1,1);
         -- put_line(")");
          first := false;
        end if;
        Append(res(perm(i+offset)+1),last(perm(i+offset)+1),point);
      end loop;
      offset := offset + SptType(i) - 1;
    end loop;
    return res;
  end Supports_of_Mixed_Cell;

  function Labels_to_Mixed_Cell
             ( nVar,nSpt : in integer32;
               SptType : in Standard_Integer_Vectors.Link_to_Vector;
               labels : in Standard_Integer_Vectors.Link_to_Vector;
               Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
               lft : in Standard_Floating_Vectors.Link_to_Vector )
             return Mixed_Cell is

    use Arrays_of_Floating_Vector_Lists;

    res : Mixed_Cell;
    pts : constant Array_of_Lists(1..nSpt)
        := Supports_of_Mixed_Cell(nVar,nSpt,SptType,labels,Vtx,lft);
    normal : constant Standard_Floating_Vectors.Vector(0..nVar-1)
           := Inner_Normal(nVar,nSpt,SptType,labels,Vtx,lft);

  begin
    res.nor := new Standard_Floating_Vectors.Vector(1..nVar+1);
    for i in normal'range loop
      res.nor(i+1) := normal(i);
    end loop;
    res.nor(nVar+1) := 1.0;
    res.pts := new Array_of_Lists'(pts);
    res.sub := null;
    return res;
  end Labels_to_Mixed_Cell;

  function Labels_to_Mixed_Cell
             ( nVar,nSpt : in integer32;
               SptType,perm : in Standard_Integer_Vectors.Link_to_Vector;
               labels : in Standard_Integer_Vectors.Link_to_Vector;
               Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
               lft : in Standard_Floating_Vectors.Link_to_Vector )
             return Mixed_Cell is

    use Arrays_of_Floating_Vector_Lists;

    res : Mixed_Cell;
    pts : constant Array_of_Lists(1..nSpt)
        := Supports_of_Mixed_Cell(nVar,nSpt,SptType,perm,labels,Vtx,lft);
    normal : constant Standard_Floating_Vectors.Vector(0..nVar-1)
           := Inner_Normal(nVar,nSpt,SptType,labels,Vtx,lft);

  begin
    res.nor := new Standard_Floating_Vectors.Vector(1..nVar+1);
    for i in normal'range loop
      res.nor(i+1) := normal(i);
    end loop;
    res.nor(nVar+1) := 1.0;
    res.pts := new Array_of_Lists'(pts);
    res.sub := null;
    return res;
  end Labels_to_Mixed_Cell;

  procedure Create_Mixed_Cell_Configuration
               ( nVar,nSpt,CellSize,nbCells : in integer32;
                 SptType : in Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 cells : in out CellStack; sub : out Mixed_Subdivision ) is

    labels : Standard_Integer_Vectors.Link_to_Vector;
    last : Mixed_Subdivision;

  begin
    for k in 1..nbCells loop
      labels := Cs_Cur(cells);
     -- put("cell "); put(k,1); put(" has labels "); put(labels); new_line;
      declare
        mic : constant Mixed_Cell
            := Labels_to_Mixed_Cell(nVar,nSpt,SptType,labels,Vtx,lft);
      begin
        Append(sub,last,mic);
      end;
      if k /= nbCells
       then Cs_Pop(cells);
      end if;
    end loop;
  end Create_Mixed_Cell_Configuration;

  procedure Create_Mixed_Cell_Configuration
               ( nVar,nSpt,CellSize,nbCells : in integer32;
                 SptType,perm : in Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                 lft : in Standard_Floating_Vectors.Link_to_Vector;
                 cells : in out CellStack; sub : out Mixed_Subdivision ) is

    labels : Standard_Integer_Vectors.Link_to_Vector;
    last : Mixed_Subdivision;

  begin
   -- put_line("inside MixedVol_Algorithm.create_mcc with perm ...");
    for k in 1..nbCells loop
      labels := Cs_Cur(cells);
     -- put("cell "); put(k,1); put(" has labels "); put(labels); new_line;
      declare
        mic : constant Mixed_Cell
            := Labels_to_Mixed_Cell(nVar,nSpt,SptType,perm,labels,Vtx,lft);
      begin
        Append(sub,last,mic);
      end;
      if k /= nbCells
       then Cs_Pop(cells);
      end if;
    end loop;
  end Create_Mixed_Cell_Configuration;

end MixedVol_Algorithm;
