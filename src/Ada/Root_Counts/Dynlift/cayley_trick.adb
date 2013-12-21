with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Dynamic32_Triangulations;  use Standard_Dynamic32_Triangulations;
with Cayley_Embedding;                   use Cayley_Embedding;
with Flatten_Mixed_Subdivisions;         use Flatten_Mixed_Subdivisions;

package body Cayley_Trick is

-- UTILITIES :

  function Extract ( n : integer32; mix : Vector; lifted : in List )
                   return Array_of_Lists is

  -- DESCRIPTION :
  --   Extracts from the list of lifted points to compute the Cayley
  --   triangulation, the tuple of lifted points.

    res : Array_of_Lists(mix'range);

  begin
    for k in res'range loop
      res(k) := Extract(k-1,n,lifted);
      Deflate(n,res(k));
    end loop;
    return res;
  end Extract;
 
  procedure Extract ( n : in integer32; mix : in Vector;
                      t : in Triangulation; liftedt : in List;
                      mixsub : out Mixed_Subdivision;
                      lifted : out Array_of_Lists ) is

  -- DESCRIPTION :
  --   Extracts the useful information from the Cayley polytope.

    res : Mixed_Subdivision;

  begin
    lifted := Extract(n,mix,liftedt);
    res := Extract_Mixed_Cells(n,mix,t);
    Deflate(n,res);
    mixsub := res;
  end Extract;

  procedure Extract_and_Clear
                ( n : in integer32; mix : in Vector;
                  t : in out Triangulation; liftedt : in out List;
                  lent : out natural32; mixsub : out Mixed_Subdivision;
                  lifted : out Array_of_Lists ) is

  -- DESCRIPTION :
  --   Extracts the useful information from the Cayley polytope.
  --   All intermediate data structures will be cleared.

  begin
    lent := Length_Of(t);
    Extract(n,mix,t,liftedt,mixsub,lifted);
    Clear(t); Clear(liftedt);
  end Extract_and_Clear;

-- BASIC VERSION :

  procedure Dynamic_Cayley
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 ) is

    tmpsub,lastcells : Mixed_Subdivision;
    L,liftedl,liftedl_last : list;
    t : Triangulation;

    procedure Col_Flat ( nt : in Triangulation; L : List ) is

    -- DESCRIPTION :
    --   Updates the subdivision mixsub with the flattened cells.
    --   The triangulation on entry contains the whole triangulation,
    --   not just the new cells.

      cells : Mixed_Subdivision;

    begin
      if Is_Null(tmpsub) then
        cells := Extract_Mixed_Cells(n,mix,nt);
        Deflate(n,cells);
      else
        cells := Extract_non_Flat_Mixed_Cells(n,mix,nt);
        Deflate(n,cells);
        Construct(Head_Of(tmpsub),cells);
      end if;
      Flatten(cells);
      tmpsub := cells;
    end Col_Flat;
    procedure C_Dynamic_Lifting is new Dynamic_Lifting_with_Flat(Col_Flat);

  begin
    L := Embedding_before_Lifting(supports);
    C_Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    if Is_Null(tmpsub) then
      Extract_and_Clear(n,mix,t,liftedl,numtri,mixsub,lifted);
    else
      lastcells := Extract_non_Flat_Mixed_Cells(n,mix,t);
      Deflate(n,lastcells);
      Construct(Head_Of(tmpsub),lastcells);
      mixsub := lastcells;
      lifted := Extract(n,mix,liftedl);
      numtri := Length_Of(t);
    end if;
  end Dynamic_Cayley;

  procedure Dynamic_Cayley
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation ) is

    L,liftedl,liftedl_last : list;

  begin
    L := Embedding_before_Lifting(supports);
    Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    lifted := Extract(n,mix,liftedl); Clear(liftedl);
  end Dynamic_Cayley;

-- EXTENDED VERSIONS :

  procedure Dynamic_Cayley_with_Flat
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 ) is

    L,liftedl,liftedl_last : list;
    t : Triangulation;
    tmpsub,lastcells : Mixed_Subdivision;

    procedure Bef_Flat ( tt : in Triangulation; lft : in List ) is

      cells,cells1 : Mixed_Subdivision;
      lftpts : Array_of_Lists(mix'range);
      
    begin
      Extract(n,mix,tt,lft,cells,lftpts);
      Before_Flattening(cells,lftpts);
      if Is_Null(tmpsub) then
        cells := Extract_Mixed_Cells(n,mix,tt);
        Deflate(n,cells);
      else
        cells := Extract_non_Flat_Mixed_Cells(n,mix,tt);
        Deflate(n,cells);
        Construct(Head_Of(tmpsub),cells);
      end if;
      Flatten(cells);
      tmpsub := cells;
    end Bef_Flat;
    procedure C_Dynamic_Lifting is new Dynamic_Lifting_with_Flat (Bef_Flat);

  begin
    L := Embedding_before_Lifting(supports);
    C_Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    if Is_Null(tmpsub) then
      Extract_and_Clear(n,mix,t,liftedl,numtri,mixsub,lifted);
    else
      lastcells := Extract_non_Flat_Mixed_Cells(n,mix,t);
      Deflate(n,lastcells);
      Construct(Head_Of(tmpsub),lastcells);
      mixsub := lastcells;
      lifted := Extract(n,mix,liftedl);
      numtri := Length_Of(t);
    end if;
  end Dynamic_Cayley_with_Flat;

  procedure Dynamic_Cayley_with_Flatt
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation ) is

    L,liftedl,liftedl_last : list;

    procedure Bef_Flat ( tt : in Triangulation; lft : in List ) is

      cells : Mixed_Subdivision;
      lftpts : Array_of_Lists(supports'range);

    begin
      Extract(n,mix,tt,lft,cells,lftpts);
      Before_Flattening(cells,lftpts);
    end Bef_Flat;
    procedure C_Dynamic_Lifting is new Dynamic_Lifting_with_Flat (Bef_Flat);

  begin
    L := Embedding_before_Lifting(supports);
    C_Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    lifted := Extract(n,mix,liftedl); Clear(liftedl);
  end Dynamic_Cayley_with_Flatt;

  procedure Dynamic_Cayley_with_New
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 ) is

    L,liftedl,liftedl_last : list;
    t : Triangulation;
    tmpsub,lastcells : Mixed_Subdivision;

    procedure Col_Flat ( nt : in Triangulation; L : List ) is

    -- DESCRIPTION :
    --   Updates the subdivision mixsub with the flattened cells.
    --   The triangulation on entry contains the whole triangulation,
    --   not just the new cells.

      cells : Mixed_Subdivision;

    begin
      if Is_Null(tmpsub) then
        cells := Extract_Mixed_Cells(n,mix,nt);
        Deflate(n,cells);
      else
        cells := Extract_non_Flat_Mixed_Cells(n,mix,nt);
        Deflate(n,cells);
        Construct(Head_Of(tmpsub),cells);
      end if;
      Flatten(cells);
      tmpsub := cells;
    end Col_Flat;

    procedure New_Cell ( tt : in Triangulation; pt : in vector ) is

      cells : Mixed_Subdivision := Extract_Mixed_Cells(n,mix,tt);
      index : integer32 := 1;

    begin
      Deflate(n,cells);
      for i in 1..mix'last-1 loop
        if pt(i+n) /= 0
         then index := i+1;
        end if;
        exit when index > 1;
      end loop;
      Process_New_Cells(cells,index,pt);
    end New_Cell;
    procedure C_Dynamic_Lifting is new Dynamic_Lifting_with_Flat_and_New
      (Before_Flattening => Col_Flat, Process_New_Simplices => New_Cell);

  begin
    L := Embedding_before_Lifting(supports);
    C_Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    if Is_Null(tmpsub) then
      Extract_and_Clear(n,mix,t,liftedl,numtri,mixsub,lifted);
    else
      lastcells := Extract_non_Flat_Mixed_Cells(n,mix,t);
      Deflate(n,lastcells);
      Construct(Head_Of(tmpsub),lastcells);
      mixsub := lastcells;
      lifted := Extract(n,mix,liftedl);
      numtri := Length_Of(t);
    end if;
  end Dynamic_Cayley_with_New;

  procedure Dynamic_Cayley_with_Newt
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation ) is

    L,liftedl,liftedl_last : list;

    procedure New_Cell ( tt : in Triangulation; pt : in vector ) is

      cells : Mixed_Subdivision := Extract_Mixed_Cells(n,mix,tt);
      index : integer32 := 1;

    begin
      Deflate(n,cells);
      for i in 1..mix'last-1 loop
        if pt(i+n) /= 0
         then index := i+1;
        end if;
        exit when index > 1;
      end loop;
      Process_New_Cells(cells,index,pt);
    end New_Cell;
    procedure C_Dynamic_Lifting is new Dynamic_Lifting_with_New(New_Cell);

  begin
    L := Embedding_before_Lifting(supports);
    C_Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    lifted := Extract(n,mix,liftedl); Clear(liftedl);
  end Dynamic_Cayley_with_Newt;

  procedure Dynamic_Cayley_with_Flat_and_New
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  mixsub : out Mixed_Subdivision; numtri : out natural32 ) is

    L,liftedl,liftedl_last : list;
    t : Triangulation;
    tmpsub,lastcells : Mixed_Subdivision;

    procedure Bef_Flat ( tt : in Triangulation; lft : in List ) is

      cells : Mixed_Subdivision;
      lftpts : Array_of_Lists(mix'range);

    begin
      Extract(n,mix,tt,lft,cells,lftpts);
      Before_Flattening(cells,lftpts);
      if Is_Null(tmpsub) then
        cells := Extract_Mixed_Cells(n,mix,tt);
        Deflate(n,cells);
      else
        cells := Extract_non_Flat_Mixed_Cells(n,mix,tt);
        Deflate(n,cells);
        Construct(Head_Of(tmpsub),cells);
      end if;
      Flatten(cells);
      tmpsub := cells;
    end Bef_Flat;

    procedure New_Cell ( tt : in Triangulation; pt : in vector ) is

      cells : Mixed_Subdivision := Extract_Mixed_Cells(n,mix,tt);
      index : integer32 := 1;

    begin
      Deflate(n,cells);
      for i in 1..mix'last-1 loop
        if pt(i+n) /= 0
         then index := i+1;
        end if;
        exit when index > 1;
      end loop;
      Process_New_Cells(cells,index,pt);
    end New_Cell;

    procedure C_Dynamic_Lifting is new Dynamic_Lifting_with_Flat_and_New
      (Before_Flattening => Bef_Flat, Process_New_Simplices => New_Cell);

  begin
    L := Embedding_before_Lifting(supports);
    C_Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    if Is_Null(tmpsub) then
      Extract_and_Clear(n,mix,t,liftedl,numtri,mixsub,lifted);
    else
      lastcells := Extract_non_Flat_Mixed_Cells(n,mix,t);
      Deflate(n,lastcells);
      Construct(Head_Of(tmpsub),lastcells);
      mixsub := lastcells;
      numtri := Length_Of(t);
    end if;
  end Dynamic_Cayley_with_Flat_and_New;

  procedure Dynamic_Cayley_with_Flat_and_Newt
                ( n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : out Array_of_Lists;
                  t : in out Triangulation ) is

    L,liftedl,liftedl_last : list;

    procedure Bef_Flat ( tt : in Triangulation; lft : in List ) is

      cells : Mixed_Subdivision;
      lftpts : Array_of_Lists(supports'range);

    begin
      Extract(n,mix,tt,lft,cells,lftpts);
      Before_Flattening(cells,lftpts);
    end Bef_Flat;

    procedure New_Cell ( tt : in Triangulation; pt : in vector ) is

      cells : Mixed_Subdivision := Extract_Mixed_Cells(n,mix,tt);
      index : integer32 := 1;

    begin
      Deflate(n,cells);
      for i in 1..mix'last-1 loop
        if pt(i+n) /= 0
         then index := i+1;
        end if;
        exit when index > 1;
      end loop;
      Process_New_Cells(cells,index,pt);
    end New_Cell;

    procedure C_Dynamic_Lifting is new Dynamic_Lifting_with_Flat_and_New
        (Before_Flattening => Bef_Flat, Process_New_Simplices => New_Cell);

  begin
    L := Embedding_before_Lifting(supports);
    C_Dynamic_Lifting(L,order,inter,maxli,liftedl,liftedl_last,t);
    lifted := Extract(n,mix,liftedl); Clear(liftedl);
  end Dynamic_Cayley_with_Flat_and_Newt;

end Cayley_Trick;
