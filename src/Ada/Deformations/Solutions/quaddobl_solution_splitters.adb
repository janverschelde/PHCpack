with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Quad_Double_Numbers_io;              use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;            use Standard_Natural_Vectors;
with QuadDobl_Condition_Tables;           use QuadDobl_Condition_Tables;
with QuadDobl_Complex_Solutions_io;       use QuadDobl_Complex_Solutions_io;
with QuadDobl_Solution_Filters;           use QuadDobl_Solution_Filters;

package body QuadDobl_Solution_Splitters is

  procedure Filter_and_Split_Solutions
              ( file : in file_type; sols : in Solution_List;
                n,k : in integer32; tol : in double_float;
                sols0,sols1 : out Solution_List ) is

    tar_sols : constant Solution_List
             := On_Target_Filter(sols,Create(integer(1)),tol);
    van_sols : constant Solution_List := Vanishing_Filter(tar_sols,tol);

  begin
    new_line(file);
    put(file,"FILTERED "); put(file,Length_Of(sols),1);
    put(file," computed vectors and found "); 
    put(file,Length_Of(van_sols),1);
    put_line(file," vanishing solutions.");
    if k = 0 then
      sols1 := van_sols;
      if not Is_Null(van_sols)
       then put(file,Length_Of(sols1),natural32(Head_Of(sols1).n),sols1);
      end if;
    else
      sols0 := Zero_Component_Filter(van_sols,natural32(n+k),tol);
      if Is_Null(sols0) then
        put(file,"NO SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," = 0.");
      else
        put(file,"THE SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," = 0 :");
        put(file,Length_Of(sols0),natural32(Head_Of(sols0).n),sols0);
      end if;
      sols1 := Free_Component_Filter(van_sols,natural32(n+k),tol);
      if Is_Null(sols1) then
        put(file,"NO SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," /= 0.");
      else
        put(file,"THE SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," /= 0 :");
        put(file,Length_Of(sols1),natural32(Head_Of(sols1).n),sols1);
      end if;
    end if;
  end Filter_and_Split_Solutions;

  procedure Filter_and_Split_Solutions
                ( sols : in Solution_List;
                  n,k : in integer32; tol : in double_float;
                  sols0,sols1 : out Solution_List ) is

    tar_sols : constant Solution_List
             := On_Target_Filter(sols,Create(integer(1)),tol);
    van_sols : constant Solution_List := Vanishing_Filter(tar_sols,tol);

  begin
    if k = 0 then
      sols1 := van_sols;
    else
      sols0 := Zero_Component_Filter(van_sols,natural32(n+k),tol);
      sols1 := Free_Component_Filter(van_sols,natural32(n+k),tol);
    end if;
  end Filter_and_Split_Solutions;

  procedure S_Singular_Filter
              ( sols : in Solution_List; tol : in double_float;
                sinsols,regsols : out Solution_List ) is

    sin_last,reg_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      if ((ls.err <= tol) or (ls.res <= tol)) then
        if ls.rco > tol
         then Append(regsols,reg_last,ls.all);
         else Append(sinsols,sin_last,ls.all);
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end S_Singular_Filter;

  procedure R_Singular_Filter
              ( file : in file_type;
                sols : in Solution_List; tol : in double_float;
                sinsols,regsols : out Solution_List ) is

    sin_last,reg_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    table : Vector(0..15) := Create(15);

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put(file,"Solution "); put(file,i,1); put(file," : ");
      put(file,"  err ="); put(file,ls.err,3);
      put(file,"  rco ="); put(file,ls.rco,3);
      put(file,"  res ="); put(file,ls.res,3);
      if ((ls.err <= tol) or (ls.res <= tol)) then
        if ls.rco <= tol then
          put_line(file,"  singular");
          Append(sinsols,sin_last,ls.all);
        else
          put_line(file,"  regular");
          Append(regsols,reg_last,ls.all);
        end if;
      else
        put_line(file,"  no solution");
      end if;
      Update_Condition(table,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
    Write_Condition_Table(file,table);
  end R_Singular_Filter;

end QuadDobl_Solution_Splitters;
