package body Generic_Lists_of_Vectors_io is

  use Vectors,Vectors_io;

  procedure get ( n,m : in natural32; L : out List ) is
  begin
    get(Standard_Input,n,m,L);
  end get;

  procedure get ( file : in file_type;
                  n,m : in natural32; L : out List ) is

    res,res_last : List;

  begin
    for i in 1..m loop
      declare
	v : Link_to_Vector;
      begin
        get(file,n,v);
        Append(res,res_last,v);
      end;
    end loop;
    L := res;
  end get;

  procedure put ( L : in List ) is
  begin
    put(Standard_Output,l);
  end put;

  procedure put ( file : in file_type; L : in List ) is

    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end put;

end Generic_Lists_of_Vectors_io;
