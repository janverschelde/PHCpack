with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Natural_Vectors;            use Standard_Natural_Vectors;
with Standard_Condition_Tables;           use Standard_Condition_Tables;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;           use Standard_Solution_Filters;

package body Standard_Solution_Splitters is

  procedure Filter_and_Split_Solutions
              ( file : in file_type; sols : in Solution_List;
                n,k : in integer32; tol : in double_float;
                sols0,sols1 : out Solution_List ) is

    tar_sols : Solution_List
             := On_Target_Filter(sols,Create(1.0),tol);
    van_sols : Solution_List := Vanishing_Filter(tar_sols,tol);

  begin
    new_line(file);
    put(file,"FILTERED "); put(file,Length_Of(sols),1);
    put(file," computed vectors and found "); 
    put(file,Length_Of(van_sols),1);
    put_line(file," vanishing solutions.");
    if k = 0 then
      sols1 := van_sols;
      if not Is_Null(van_sols) then
        put_line(file,"FILTERED SOLUTIONS :");
        put(file,Length_Of(sols1),natural32(Head_Of(sols1).n),sols1);
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
      Clear(van_sols);
    end if;
    Clear(tar_sols);
  end Filter_and_Split_Solutions;

  procedure Zero_Singular_Split_Filter
              ( file : in file_type; sols : in Solution_List;
                n,k : in integer32; tolzero,tolsing : in double_float;
                sols0,sols1 : out Solution_List ) is

    tar_sols : Solution_List
             := On_Target_Filter(sols,Create(1.0),tolzero);
    van_sols : Solution_List := Vanishing_Filter(tar_sols,tolzero);
    regsols,sinsols : Solution_List;
    diflen : natural32;

  begin
    new_line(file);
    put(file,"FILTERED "); put(file,Length_Of(sols),1);
    put(file," computed vectors and found "); 
    put(file,Length_Of(van_sols),1);
    put_line(file," vanishing solutions.");
    if k = 0 then
      sols1 := van_sols;
      if not Is_Null(van_sols) then
        put_line(file,"FILTERED SOLUTIONS :");
        put(file,Length_Of(sols1),natural32(Head_Of(sols1).n),sols1);
      end if;
    else
      sols0 := Zero_Component_Filter(van_sols,natural32(n+k),tolzero);
      if Is_Null(sols0) then
        put(file,"NO SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," = 0.");
      else
        put(file,"THE SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," = 0 :");
        put(file,Length_Of(sols0),natural32(Head_Of(sols0).n),sols0);
      end if;
      sols1 := Free_Component_Filter(van_sols,natural32(n+k),tolzero);
      if Is_Null(sols1) then
        put(file,"NO SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," /= 0.");
      else
        Silent_Singular_Filter(sols1,tolsing,sinsols,regsols);
        diflen := Length_Of(sols1) - Length_Of(regsols);
        if diflen > 0 then
          put(file,"Removed "); put(file,diflen,1);
          put_line(file," singular solutions with nonzero slack variable.");
          put(file,"THE SINGULAR SOLUTIONS WITH zz");
          put(file,k,1); put_line(file," /= 0 :");
          put(file,Length_Of(sinsols),natural32(Head_Of(sinsols).n),sinsols);
        end if;
        Clear(sols1); sols1 := regsols;
        put(file,"THE REGULAR SOLUTIONS WITH zz");
        put(file,k,1); put_line(file," /= 0 :");
        put(file,Length_Of(sols1),natural32(Head_Of(sols1).n),sols1);
        Clear(sinsols);
      end if;
      Clear(van_sols);
    end if;
    Clear(tar_sols);
  end Zero_Singular_Split_Filter;

  procedure Filter_and_Split_Solutions
                ( sols : in Solution_List;
                  n,k : in integer32; tol : in double_float;
                  sols0,sols1 : out Solution_List ) is

    tar_sols : Solution_List
             := On_Target_Filter(sols,Create(1.0),tol);
    van_sols : Solution_List := Vanishing_Filter(tar_sols,tol);

  begin
    if k = 0 then
      sols1 := van_sols;
    else
      sols0 := Zero_Component_Filter(van_sols,natural32(n+k),tol);
      sols1 := Free_Component_Filter(van_sols,natural32(n+k),tol);
      Clear(van_sols);  
    end if;
    Clear(tar_sols);
  end Filter_and_Split_Solutions;

  procedure Zero_Singular_Split_Filter
                ( sols : in Solution_List;
                  n,k : in integer32; tolzero,tolsing : in double_float;
                  sols0,sols1 : out Solution_List ) is

    tar_sols : Solution_List
             := On_Target_Filter(sols,Create(1.0),tolzero);
    van_sols : Solution_List := Vanishing_Filter(tar_sols,tolzero);
    regsols,sinsols : Solution_List;

  begin
    if k = 0 then
      sols1 := van_sols;
    else
      sols0 := Zero_Component_Filter(van_sols,natural32(n+k),tolzero);
      sols1 := Free_Component_Filter(van_sols,natural32(n+k),tolzero);
      if not Is_Null(sols1) then
        Silent_Singular_Filter(sols1,tolsing,sinsols,regsols);
        Clear(sols1); sols1 := regsols;
        Clear(sinsols);  
      end if;
      Clear(van_sols);
    end if;
    Clear(tar_sols);
  end Zero_Singular_Split_Filter;

  procedure Silent_Singular_Filter
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
  end Silent_Singular_Filter;

  procedure Reporting_Singular_Filter
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
  end Reporting_Singular_Filter;

end Standard_Solution_Splitters;
