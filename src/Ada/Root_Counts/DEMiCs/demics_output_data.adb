with String_Splitters;
with Standard_Floating_Vectors;

package body DEMiCs_Output_Data is

  lifting : Standard_Floating_VecVecs.Link_to_VecVec;
  first,last : Lists_of_Strings.List;

  procedure Initialize_Lifting
              ( crdsup : in Standard_Integer_Vectors.Vector ) is
  begin
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

  function Retrieve_Lifting
             ( idxsup,idxpnt : integer32) return double_float is
  begin
    return lifting(idxsup)(idxpnt);
  end Retrieve_Lifting;

  function Lifting_Values return Standard_Floating_VecVecs.Link_to_VecVec is
  begin
    return lifting;
  end Lifting_Values;

  procedure Add_Cell_Indices ( strcell : in string ) is

    link2strcell : constant String_Splitters.Link_to_String
                 := new string'(strcell);

  begin
    Lists_of_Strings.Append(first,last,link2strcell);
  end Add_Cell_Indices;

  function Retrieve_Cell_Indices return Lists_of_Strings.List is
  begin
    return first;
  end Retrieve_Cell_Indices;

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
