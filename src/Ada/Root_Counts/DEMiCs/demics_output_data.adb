with text_io;
with Standard_Floating_Vectors;

package body DEMiCs_Output_Data is

-- DATA STRUCTURES :

  lifting : Standard_Floating_VecVecs.Link_to_VecVec := null;
  first,last,cellptr : Lists_of_Strings.List;

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

  procedure Add_Cell_Indices ( strcell : in string ) is

    link2strcell : constant String_Splitters.Link_to_String
                 := new string'(strcell);

  begin
    Lists_of_Strings.Append(first,last,link2strcell);
    if monitor
     then text_io.put_line(strcell);
    end if;
  end Add_Cell_Indices;

  procedure Initialize_Cell_Pointer is
  begin
    cellptr := first;
  end Initialize_Cell_Pointer;

-- SELECTORS :

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
    return Lists_of_Strings.Length_Of(first);
  end Number_of_Cell_Indices;

  function Get_Cell_Indices
             ( index : integer32 ) return String_Splitters.Link_to_String is

    res : String_Splitters.Link_to_String := null;
    tmp : Lists_of_Strings.List := first;
    cnt : integer32 := 0;

  begin
    while not Lists_of_Strings.Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = index then
        res := Lists_of_Strings.Head_Of(tmp);
        return res;
      end if;
      tmp := Lists_of_Strings.Tail_Of(tmp);
    end loop;
    return res;
  end Get_Cell_Indices;

  function Retrieve_Cell_Indices return Lists_of_Strings.List is
  begin
    return first;
  end Retrieve_Cell_Indices;

  function Get_Next_Cell return String_Splitters.Link_to_String is

    res : String_Splitters.Link_to_String := null;

  begin
    if not Lists_of_Strings.Is_Null(cellptr) then
      res := Lists_of_Strings.Head_Of(cellptr);
      cellptr := Lists_of_Strings.Tail_Of(cellptr);
    end if;
    return res;
  end Get_Next_Cell;

-- DESTRUCTORS :

  procedure Clear_Lifting is
  begin
    Standard_Floating_VecVecs.Deep_Clear(lifting);
  end Clear_Lifting;

  procedure Clear_Cell_Indices is

    tmp : Lists_of_Strings.List := first;
    ls : String_Splitters.Link_to_String;

  begin
    while not Lists_of_Strings.Is_Null(tmp) loop
      ls := Lists_of_Strings.Head_Of(tmp);
      String_Splitters.Clear(ls);
      tmp := Lists_of_Strings.Tail_Of(tmp);
    end loop;
    Lists_of_Strings.Clear(first);
    last := first;
  end Clear_Cell_Indices;

end DEMiCs_Output_Data;
