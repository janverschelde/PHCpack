with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Solutions_Pool;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;

function use_solpool ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is -- initialize pool with n = a[0]

    v : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(v(v'first));

  begin
    Solutions_Pool.Initialize(n);
    return 0;
  end Job0;

  function Job1 return integer32 is -- returns the size of the pool in a[0]

    n : constant natural32 := Solutions_Pool.Size;

  begin
    Assign(integer32(n),a);
    return 0;
  end Job1;

  function Job2 return integer32 is -- returns length of solution list 

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    res : constant natural32 := Solutions_Pool.Length(k);

  begin
    Assign(integer32(res),b);
    return 0;
  end Job2;

  function Job3 return integer32 is -- returns dimension of solution list 

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    res : constant natural32 := Solutions_Pool.Dimension(k);

  begin
    Assign(integer32(res),b);
    return 0;
  end Job3;

  function Job4 return integer32 is -- appends a solution to a list
  
    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
   -- put("Thread "); put(k,1); put_line(" appends a solution ...");
    Solutions_Pool.Append(k,ls);
   -- delay(1.0);
   -- put("Number of solutions in pool "); put(k,1);
   -- put(" is "); put(Solutions_Pool.Length(k),1); new_line;
    return 0;
  end Job4;

  function Job5 return integer32 is -- retrieves a solution from a list

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant integer32 := integer32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;

  begin
    Solutions_Pool.Retrieve(k,i,ls,fail);
    if fail then
      return 325;
    else
      Assign_Solution(ls,b,c);
      return 0;
    end if;
  end Job5;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- initialize pool with a[0]
      when 1 => return Job1; -- returns size of pool in a[0]
      when 2 => return Job2; -- returns length of solution list
      when 3 => return Job3; -- returns dimension of solution list
      when 4 => return Job4; -- appends a solution to a list
      when 5 => return Job5; -- retrieves a solution from a list
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_solpool;
