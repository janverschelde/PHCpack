with Semaphore;

package body Static_Columns_Queue is

  sem : Semaphore.Lock;
  idx : integer32 := 0;
  col : Standard_Integer_VecVecs.Link_to_VecVec;

  procedure Initialize ( stc : in Standard_Integer_VecVecs.VecVec ) is
  begin
    col := new Standard_Integer_VecVecs.VecVec'(stc);
    idx := 0;
  end Initialize;

  function Next_Columns return Standard_Integer_Vectors.Link_to_Vector is

    res : Standard_Integer_Vectors.Link_to_Vector := null;

  begin
    Semaphore.Request(sem);
    idx := idx + 1;
    if idx <= col'last
     then res := col(idx);
    end if;
    Semaphore.Release(sem);
    return res;
  end Next_Columns;

  function Next_Counter return integer32 is
  begin
    return idx;
  end Next_Counter;

  procedure Clear is
  begin
    Standard_Integer_VecVecs.Shallow_Clear(col);
    idx := 0;
  end Clear;

end Static_Columns_Queue;
