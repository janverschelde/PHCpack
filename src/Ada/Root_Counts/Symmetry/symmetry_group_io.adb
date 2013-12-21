with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;

package body Symmetry_Group_io is

  procedure get ( p : out Permutation ) is
  begin
    get(Standard_Input,p);
  end get;

  procedure get ( file : in file_type; p : out Permutation ) is
  begin
    for i in p'range loop
      p(i) := 0; get(file,p(i));
    end loop;
  end get;

  procedure get ( L : in out List_of_Permutations; n,nb : in integer32 ) is
  begin
    get(Standard_Input,l,n,nb);
  end get;

  procedure get ( file : in file_type;
		  L : in out List_of_Permutations; n,nb : in integer32 ) is

    p : Permutation(1..n);
    l2 : List_of_Permutations;

  begin
    for i in 1..nb loop
      get(file,p);
      if Is_Permutation(p)
       then Append(l,l2,p);
      end if;
    end loop;
  end get;
      
  procedure put ( p : in Permutation ) is
  begin
    put(Standard_Output,p);
  end put;

  procedure put ( file : in file_type; p : in Permutation ) is
  begin
    for i in p'range loop
      put(file,' '); put(file,p(i),1);
    end loop;
  end put;

  procedure put ( L : in List_of_Permutations ) is
  begin
    put(Standard_Output,l);
  end put;

  procedure put ( file : in file_type; L : in List_of_Permutations ) is

    temp : List_of_Permutations := L;

  begin
    while not Is_Null(temp) loop
      put(file,Permutation(Head_Of(temp).all));
      new_line(file);
      temp := Tail_Of(temp);
    end loop;
  end put;

end Symmetry_Group_io;
