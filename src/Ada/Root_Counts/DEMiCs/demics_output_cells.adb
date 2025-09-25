with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;   use Arrays_of_Floating_Vector_Lists;

package body DEMiCs_Output_Cells is

-- DATA STRUCTURES :

  lifting : Standard_Floating_VecVecs.Link_to_VecVec := null;
  first,last,cellptr : Lists_of_Integer_Vectors.List;
  cells,cells_last,allcellptr : Mixed_Subdivision;
  dimension : integer32 := -1;
  mixture : Standard_Integer_Vectors.Link_to_Vector;
  setcellptr,setallcellptr : boolean;
  labels_size : integer32 := -1;
  cellcnt : integer32 := 0;

-- CONSTRUCTORS :

  procedure Initialize_Lifting
              ( crdsup : in Standard_Integer_Vectors.Vector ) is

    use Standard_Floating_VecVecs;

  begin
    if lifting /= null
     then Clear_Lifting;
    end if;
    lifting := new Standard_Floating_VecVecs.VecVec(crdsup'range);
    for i in crdsup'range loop
      declare
        zeros : Standard_Floating_Vectors.Vector(1..crdsup(i));
      begin
        zeros := (zeros'range => 0.0);
        lifting(i) := new Standard_Floating_Vectors.Vector'(zeros);
      end;
    end loop;
  end Initialize_Lifting;

  procedure Assign_Lifting 
              ( idxsup,idxpnt : in integer32; val : in double_float ) is
  begin
    lifting(idxsup)(idxpnt) := val;
  end Assign_Lifting;

  procedure Store_Dimension_and_Mixture
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector ) is
  begin
    dimension := dim;
    mixture := mix;
    labels_size := 0;
    for i in mix'range loop
      labels_size := labels_size + 1 + mix(i);
    end loop;
    cellcnt := 0;
  end Store_Dimension_and_Mixture;

  procedure Allocate_Mixed_Cell is

  -- FOR MULTITASKING :
  --   The allocation is done by the producer task in a 2-stage pipeline.
  --   In this pipeline, only the producer allocates, while all other
  --   consumer tasks have allocated their memory space already.

    mic : Mixed_Cell;
    zeros : constant Standard_Floating_Vectors.Vector(1..dimension+1)
          := (1..dimension+1 => 0.0);
    last : Lists_of_Floating_Vectors.List;

  begin
    mic.nor := new Standard_Floating_Vectors.Vector'(1..dimension+1 => 0.0);
    mic.pts := new Array_of_Lists(mixture'range);
    for i in mic.pts'range loop
      last := mic.pts(i);
      for j in 1..mixture(i)+1 loop
        Lists_of_Floating_Vectors.Append(mic.pts(i),last,zeros);
      end loop; 
    end loop;
    mic.sub := null;
    Append(cells,cells_last,mic);
  end Allocate_Mixed_Cell;

  procedure Add_Cell_Indices
              ( idx : in Standard_Integer_Vectors.Link_to_Vector ) is

    idxcnt : integer32;

  begin
    if allocate 
     then Allocate_Mixed_Cell;
    end if;
    Lists_of_Integer_Vectors.Append(first,last,idx);
    if monitor then
      cellcnt := cellcnt + 1;
      put("# "); put(cellcnt,1); put(" : ");
      idxcnt := 0;
      for i in mixture'range loop
        put(i,1); put(" : (");
        for j in 1..mixture(i)+1 loop
          idxcnt := idxcnt + 1;
          put(" "); put(idx(idxcnt),1);
        end loop;
        put(" ) ");
      end loop;
      new_line;
    end if;
  end Add_Cell_Indices;

  procedure Initialize_Cell_Indices_Pointer is
  begin
    cellptr := first;
  end Initialize_Cell_Indices_Pointer;

  procedure Initialize_Allocated_Cell_Pointer is
  begin
    allcellptr := cells;
  end Initialize_Allocated_Cell_Pointer;

-- SELECTORS :

  function Get_Labels_Size return integer32 is
  begin
    return labels_size;
  end Get_Labels_Size;

  function Get_Mixture return Standard_Integer_Vectors.Link_to_Vector is
  begin
    return mixture;
  end Get_Mixture;

  function Retrieve_Lifting
             ( idxsup,idxpnt : integer32) return double_float is
  begin
    return lifting(idxsup)(idxpnt);
  end Retrieve_Lifting;

  function Lifting_Values return Standard_Floating_VecVecs.Link_to_VecVec is
  begin
    return lifting;
  end Lifting_Values;

  function Number_of_Cell_Indices return natural32 is
  begin
    return Lists_of_Integer_Vectors.Length_Of(first);
  end Number_of_Cell_Indices;

  function Get_Cell_Indices
             ( index : integer32 )
             return Standard_Integer_Vectors.Link_to_Vector is

    res : Standard_Integer_Vectors.Link_to_Vector := null;
    tmp : Lists_of_Integer_Vectors.List := first;
    cnt : integer32 := 0;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = index then
        res := Lists_of_Integer_Vectors.Head_Of(tmp);
        return res;
      end if;
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Get_Cell_Indices;

  function Retrieve_Cell_Indices return Lists_of_Integer_Vectors.List is
  begin
    return first;
  end Retrieve_Cell_Indices;

  function Get_Next_Cell_Indices
             return Standard_Integer_Vectors.Link_to_Vector is

  -- INVARIANT :
  --   Either cellptr is null,
  --   or cellptr points to the cell that was last returned.
  --   When multitasking, this function is executed in a critical section,
  --   i.e.: only one task can access the data manipulated by this function.

    res : Standard_Integer_Vectors.Link_to_Vector := null;
    tailofcellptr : Lists_of_Integer_Vectors.List;

  begin
    if Lists_of_Integer_Vectors.Is_Null(cellptr) then
      if not setcellptr then        -- cell pointer not set to first yet
        if not Lists_of_Integer_Vectors.Is_Null(first) then -- return 1st cell
          res := Lists_of_Integer_Vectors.Head_Of(first);
          cellptr := first;
          setcellptr := true;    -- we will have returned the first cell 
        end if;
      end if;
    else  -- if the tail of cellptr is not null, then we return its cell
      tailofcellptr := Lists_of_Integer_Vectors.Tail_Of(cellptr);
      if not Lists_of_Integer_Vectors.Is_Null(tailofcellptr) then
        res := Lists_of_Integer_Vectors.Head_Of(tailofcellptr);
        cellptr := tailofcellptr; -- invariant: cellptr -> last returned
      end if;
    end if;
    return res;
  end Get_Next_Cell_Indices;

  function Get_Allocated_Cells return Mixed_Subdivision is
  begin
    return cells;
  end Get_Allocated_Cells;
  
  function Get_Next_Allocated_Cell return Mixed_Subdivision is

  -- INVARIANT :
  --   Either allcellptr is null,
  --   or allcellptr points to the cell that was last returned.
  --   When multitasking, this function is executed in a critical section,
  --   i.e.: only one task can access the data manipulated by this function.

    res : Mixed_Subdivision;
    tailofallcellptr : Mixed_Subdivision;

  begin
    if Is_Null(allcellptr) then
      if not setallcellptr then       -- if allcellptr not set to cells yet
        if not Is_Null(cells) then    -- then we will return the first cell 
          res := cells;
          allcellptr := cells;
          setallcellptr := true;              -- the first cell is returned
        end if;
      end if;
    else  -- if the tail of allcellptr is not null, then we return its cell
      tailofallcellptr := Tail_Of(allcellptr);
      if not Is_Null(tailofallcellptr) then
        res := tailofallcellptr;
        allcellptr := tailofallcellptr; -- invariant : allcellptr points to 
      end if;                           -- the last returned cell
    end if;
    return res;
  end Get_Next_Allocated_Cell;

-- DESTRUCTORS :

  procedure Clear_Lifting is
  begin
    Standard_Floating_VecVecs.Deep_Clear(lifting);
  end Clear_Lifting;

  procedure Clear_Cell_Indices is

    tmp : Lists_of_Integer_Vectors.List := first;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      lpt := Lists_of_Integer_Vectors.Head_Of(tmp);
      Standard_Integer_Vectors.Clear(lpt);
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    Lists_of_Integer_Vectors.Clear(first);
    last := first;
  end Clear_Cell_Indices;

  procedure Clear_Allocated_Cells is
  begin
    Clear(cells);
  end Clear_Allocated_Cells;

begin
  setcellptr := false;
  setallcellptr := false;
end DEMiCs_Output_Cells;
