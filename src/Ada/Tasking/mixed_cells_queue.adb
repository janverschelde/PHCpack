with Semaphore;

package body Mixed_Cells_Queue is

  ptr : Mixed_Subdivision;
  cnt : integer32 := 0;
  sem : Semaphore.Lock;

  -- CORRESPONDENCE BETWEEN cnt AND ptr :
  --   cnt = 0                  <=> ptr not set to sols yet
  --   cnt in 1..Length_Of(ptr) <=> ptr points to current cell
  --   cnt = Length_Of(ptr) + 1 <=> Is_Null(ptr)

  procedure Initialize ( cells : in Mixed_Subdivision ) is
  begin
    ptr := cells;
    cnt := 0;
  end Initialize;

  function Next return Mixed_Subdivision is

    res : Mixed_Subdivision := ptr;

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

end Mixed_Cells_Queue;
