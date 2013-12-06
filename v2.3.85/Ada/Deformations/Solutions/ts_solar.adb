with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Random_Vectors;          use Standard_Random_Vectors;
with Standard_Complex_Solutions;       use Standard_Complex_Solutions;
with Standard_Solution_Array_Lists;    use Standard_Solution_Array_Lists;

procedure ts_solar is

-- DESCRIPTION :
--   Test on lists of solution arrays.

  function Random_Solution ( n : natural32 ) return Solution is

  -- DESCRIPTION :
  --   Returns a random solution vector of size n.

    res : Solution(integer32(n));

  begin
    res.t := Create(0.0);
    res.m := 1;
    res.v := Random_Vector(1,integer32(n));
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Random_Solution;

  function Random_Solutions ( m,n : natural32 ) return Solution_List is

  -- DESCRIPTION :
  --   Generates a solution list of m vector of size n.

    res,res_last : Solution_List;

  begin
    for i in 1..m loop
      Append(res,res_last,Random_Solution(n));
    end loop;
    return res;
  end Random_Solutions;

  procedure Test_Creation ( m,n,size : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the creation of a list of solution arrays of the given size,
  --   from a list of m random solution vectors of dimension n.

    sols : constant Solution_List := Random_Solutions(m,n);
    sal : constant Solution_Array_List := Create(sols,size);
    tmp : Solution_Array_List := sal;

  begin
    put("Length of solution arrays :");
    while not Is_Null(tmp) loop
      declare
        lsa : constant Link_to_Solution_Array := Head_Of(tmp);
        len : constant natural32 := natural32(lsa'length);
      begin
        put(" "); put(len,1);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    new_line;
  end Test_Creation;

  procedure Main is

    m,n,size : natural32 := 0;

  begin
    new_line;
    put_line("Generating a list of random solutions...");
    new_line;
    put("Give the length of solution list : "); get(m);
    put("Give the dimension of vectors : "); get(n);
    put("Give size of arrays of solutions : "); get(size);
    Test_Creation(m,n,size);
  end Main;

begin
  Main;
end ts_solar;
