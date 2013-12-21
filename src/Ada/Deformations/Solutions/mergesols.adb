with text_io,integer_io;                 use text_io,integer_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;

procedure mergesols is

-- DESCRIPTION :
--   Reads two solution lists from file.
--   Writes all different solutions to file.

  function Pointer_to_Last ( l : Solution_List ) return Solution_List is

  begin
    if Is_Null(l)
     then return l;
     elsif Is_Null(Tail_Of(l))
         then return l;
         else return Pointer_to_Last(Tail_Of(l));
    end if;
  end Pointer_to_Last;

  function Merge_Solutions ( l1,l2 : Solution_List; tol : double_float )
                           return Solution_List is

    res,res_last,tmp : Solution_List;
    sol : Link_to_Solution;

  begin
    Copy(l1,res);
    res_last := Pointer_to_Last(res);
    tmp := l2;
    while not Is_Null(tmp) loop
      sol := Head_Of(tmp);
      if not Is_In(l1,sol.all,tol)
       then Append(res,res_last,sol.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Merge_Solutions;

  procedure Main is

    infile,outfile : file_type;
    sols1,sols2 : Solution_List;
    tol : double_float;

  begin
    new_line;
    put_line("Merging two lists of solutions.");
    new_line;
    put_line("Reading the name of the file for the first list.");
    Read_Name_and_Open_File(infile);
    get(infile,sols1);
    Close(infile);
    new_line;
    put_line("Reading the name of the file for the second list.");
    Read_Name_and_Open_File(infile);
    get(infile,sols2);
    Close(infile);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Give the clustering tolerance : "); get(tol);
    declare
	  sols : Solution_List := Merge_Solutions(sols1,sols2,tol);
    begin
      put("There are "); put(Length_Of(sols),1);
      put_line(" distinct solutions.");
      put(outfile,Length_Of(sols),Head_Of(sols).n,sols);
    end;
    Close(outfile);
  end Main;

begin
  Main;
end mergesols;
