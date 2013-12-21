with unchecked_deallocation;

package body Cell_Stack is

  procedure Cell_Init1 ( c : in out Link_to_Cell ) is
  begin
    c.next := null;
  end Cell_Init1;

  procedure Cell_Init2 
               ( c : in out Link_to_Cell; n : in integer32;
                 J : in Standard_Integer_Vectors.Link_to_Vector;
                 ptr : in Link_to_Cell ) is
  begin
    c.idx := new Standard_Integer_Vectors.Vector(0..n-1);
    for i in 0..n-1 loop
      c.idx(i) := J(i);
    end loop;
    c.next := ptr;
  end Cell_Init2;

  procedure Cs_Init ( cs : in out CellStack; n : in integer32 ) is
  begin
    cs.size := n;
    cs.count := 0;
    cs.cur := null;
    cs.top := null;
  end Cs_Init;
     
  procedure Cs_Push ( cs : in out CellStack;
                      J : in Standard_Integer_Vectors.Link_to_Vector ) is

    c : Link_to_Cell := new Cell;

  begin
    Cell_Init2(c,cs.size,J,cs.top);
    cs.cur := c;
    cs.top := c;
    cs.count := cs.count + 1;
  end Cs_Push;

  procedure Cs_Pop ( cs : in out CellStack ) is

    ptr : Link_to_Cell := cs.top;  

  begin
    cs.top := cs.top.next;
    cs.cur := cs.top;
    Standard_Integer_Vectors.Clear(ptr.idx);
    Clear(ptr);
    cs.count := cs.count - 1;
  end Cs_Pop;
       
  procedure Cs_Next ( cs : in out CellStack; okay : out integer32 ) is
  begin
    if cs.cur.next /= null
     then cs.cur := cs.cur.next;
          okay := 1;
     else okay := 0;
    end if;
  end Cs_Next;

  function Cs_Cur ( cs : CellStack )
                  return Standard_Integer_Vectors.Link_to_Vector is
  begin
    return cs.cur.idx;
  end Cs_Cur;

  procedure Cs_Top ( cs : in out CellStack ) is
  begin
    cs.cur := cs.top;
  end Cs_Top;
      
  function Cs_IsEmpty ( cs : CellStack ) return boolean is
  begin
    if cs.top = null
     then return true;
     else return false;
    end if;
  end Cs_IsEmpty;

  function Cs_Count ( cs : CellStack ) return integer32 is
  begin
    return cs.count;
  end Cs_Count;
      
  procedure Csi ( cs : in out CellStack; i : in integer32;
                  J : out Standard_Integer_Vectors.Link_to_Vector ) is

    k : integer32;

  begin
    if (i < 0 or i >= cs.count) then
      J := null;
    else
      cs.cur := cs.top;
      k := 0;
      while k < i loop
        cs.cur := cs.cur.next;
        k := k + 1;
      end loop;
      J := cs.cur.idx;
    end if;
  end Csi;

  procedure Cs_Del ( cs : in out CellStack ) is
  begin
    while not Cs_IsEmpty(cs) loop
      Cs_Pop(cs);
    end loop;
    cs.count := 0;
  end Cs_Del;

  procedure Clear ( c : in out Link_to_Cell ) is

    procedure free is new unchecked_deallocation(Cell,Link_to_Cell);

  begin
    free(c);
  end Clear;

end Cell_Stack;
