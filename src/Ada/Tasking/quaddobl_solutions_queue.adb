with Semaphore;

package body QuadDobl_Solutions_Queue is

  ptr : Solution_List;
  cnt : integer32 := 0;
  sem : Semaphore.Lock;

  -- CORRESPONDENCE BETWEEN cnt AND ptr :
  --   cnt = 0                  <=> ptr not set to sols yet
  --   cnt in 1..Length_Of(ptr) <=> ptr points to current solution
  --   cnt = Length_Of(ptr) + 1 <=> Is_Null(ptr)

  procedure Initialize ( sols : in Solution_List ) is
  begin
    ptr := sols;
    cnt := 0;
  end Initialize;

  function Next return Solution_List is

    res : Solution_List := ptr;

  begin
    Semaphore.Request(sem);
    if cnt = 0 then
      cnt := 1;
    else
      cnt := cnt + 1;
      if not Is_Null(ptr)
       then ptr := Tail_Of(ptr);
      end if;
    end if;
    res := ptr;
    Semaphore.Release(sem);
    return res;
  end Next;

  function Next_Counter return integer32 is
  begin
    return cnt;
  end Next_Counter;

end QuadDobl_Solutions_Queue;
