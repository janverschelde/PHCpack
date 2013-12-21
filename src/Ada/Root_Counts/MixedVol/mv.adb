with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Cell_Stack;                        use Cell_Stack;
with Prepare_for_MV,Mixed_Volume;       use Prepare_for_MV,Mixed_Volume;

-- for testing :
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;

function mv ( fn : in string; nVar,nPts : in integer32;
              ind,cnt,sup : in Standard_Integer_Vectors.Vector )
            return natural32 is

  procedure Write_Supports 
                ( nSpt : in integer32;
                  Idx : in Standard_Integer_Vectors.Link_to_Vector;
                  Spt : in Standard_Integer_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Writes the supports to screen.

  -- ON ENTRY :
  --   nSpt       number of supports;
  --   Idx        index set to the supports;
  --   Spt        coordinates of the supports.

  begin
    put("The index set : "); put(Idx.all); new_line;
    for i in 1..nSpt loop
      put("Support "); put(i,1); put_line(" :");
      for j in Idx(i-1)..Idx(i)-1 loop
        put(Spt(j)); new_line;
      end loop;
    end loop;
  end Write_Supports;

  function max ( x,y : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the maximum of x and y.

  begin
    if x >= y
     then return x;
     else return y;
    end if;
  end max;

  function min ( x,y : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the minimum of x and y.

  begin
    if x <= y
     then return x;
     else return y;
    end if;
  end min;

  function quick_return
             ( nVar : integer32; 
               SptIdx : Standard_Integer_Vectors.Link_to_Vector;
               Spt : Standard_Integer_VecVecs.Link_to_VecVec )
             return integer32 is

  -- DESCRITION :
  --   Returns 1 if the system is univariate, or linear, or if some
  --   supports have fewer than two terms in them, printing a message;
  --   otherwise 0 is returned.

  begin
    return 0;
  end quick_return;

  function Main return natural32 is

    SptIdx : Standard_Integer_Vectors.Link_to_Vector
           := new Standard_Integer_Vectors.Vector(0..nVar);
    Spt : Standard_Integer_VecVecs.Link_to_VecVec
        := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
    SptType,VtxIdx,NuIdx2OldIdx : Standard_Integer_Vectors.Link_to_Vector;
    perm : Standard_Integer_Vectors.Link_to_Vector;
    Vtx : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;
    MCells : CellStack;
    ind,nSpt,CellSize,nbCells : integer32;
    MVol : natural32;

  begin
    put("The vector cnt : "); put(cnt); new_line;
    for i in 0..(nVar-1) loop
      SptIdx(i) := cnt(i+1);
    end loop;
    SptIdx(nVar) := nPts;
    for i in reverse 0..nVar-1 loop
      SptIdx(i) := SptIdx(i+1) - SptIdx(i);
    end loop;
    ind := 0;
    for i in Spt'range loop
      Spt(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
      for j in 0..(nVar-1) loop
        ind := ind + 1;
        Spt(i)(j) := sup(ind);
      end loop;
    end loop;
    if quick_return(nVar,SptIdx,Spt) = 1 then return 0; end if;
    SptType := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    VtxIdx := new Standard_Integer_Vectors.Vector(0..nVar);
    Vtx := new Standard_Integer_VecVecs.VecVec(0..(nPts-1));
    for i in 0..(nPts-1) loop
      Vtx(i) := new Standard_Integer_Vectors.Vector(0..(nVar-1));
    end loop;
    NuIdx2OldIdx := new Standard_Integer_Vectors.Vector(0..(nPts-1));
    nSpt := nVar;
    put_line("Supports before Pre4MV :"); Write_Supports(nSpt,SptIdx,Spt);
    Pre4MV(nVar,nSpt,nSpt,SptType,Spt,SptIdx,Vtx,VtxIdx,NuIdx2OldIdx,perm);
    put_line("Supports after Pre4MV :"); Write_Supports(nSpt,SptIdx,Spt);
    put_line("Vertices of supports :"); Write_Supports(nSpt,VtxIdx,Vtx);
    lft := new Standard_Floating_Vectors.Vector(0..(VtxIdx(nSpt)-1));
    for i in 0..(VtxIdx(nSpt)-1) loop
      lft(i) := 2.0*(1.5 + Standard_Random_Numbers.Random);
     -- put("Give lifting value for "); put(i,1); put(" : ");
     -- get(lft(i));
    end loop;
    put_line("The lifting values : "); put_line(lft);
    CellSize := cell_size(nSpt,SptType);
    Cs_Init(MCells,CellSize);
    MixedVol(nVar,nSpt,CellSize,SptType,VtxIdx,Vtx,lft,nbCells,MCells,MVol);
    write_cells(fn,nVar,nSpt,SptType,Vtx,lft,CellSize,nbCells,MCells);
    put("The mixed volume of this support is "); put(MVol,1);
    put_line(".");
    put("See the file " & fn);
    put_line(" for a regular mixed-cell configuration.");
    while not Cs_IsEmpty(MCells) loop               -- clear the memory
      Cs_Pop(MCells);
    end loop;
    Standard_Integer_Vectors.Clear(SptIdx);
    Standard_Integer_Vectors.Clear(VtxIdx);
    Standard_Integer_Vectors.Clear(SptType);
    Standard_Integer_Vectors.Clear(NuIdx2OldIdx);
    Standard_Integer_Vectors.Clear(perm);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(Spt);
    Standard_Integer_VecVecs.Deep_Clear(Vtx);
    return 0;
  end Main;

begin
  return Main;
end mv;
