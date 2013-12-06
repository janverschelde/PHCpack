with integer_io;                          use integer_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;           use Standard_Solution_Filters;

procedure Filter_and_Split_Solutions
              ( file : in file_type; sols : in Solution_List;
                n,k : in natural; tol : in double_float;
                sols0,sols1 : out Solution_List ) is

  tar_sols : Solution_List := On_Target_Filter(sols,Create(1.0),tol);
  van_sols : Solution_List := Vanishing_Filter(tar_sols,tol);

begin
  new_line(file);
  put(file,"FILTERED "); put(file,Length_Of(sols),1);
  put(file," computed vectors and found "); 
  put(file,Length_Of(van_sols),1);
  put_line(file," vanishing solutions.");
  if k = 0
   then sols1 := van_sols;
        if not Is_Null(van_sols)
         then put(file,Length_Of(sols1),Head_Of(sols1).n,sols1);
        end if;
   else sols0 := Zero_Component_Filter(van_sols,n+k,tol);
        if Is_Null(sols0)
         then put(file,"NO SOLUTIONS WITH zz");
              put(file,k,1); put_line(file," = 0.");
         else put(file,"THE SOLUTIONS WITH zz");
              put(file,k,1); put_line(file," = 0 :");
              put(file,Length_Of(sols0),Head_Of(sols0).n,sols0);
        end if;
        sols1 := Free_Component_Filter(van_sols,n+k,tol);
        if Is_Null(sols1)
         then put(file,"NO SOLUTIONS WITH zz");
              put(file,k,1); put_line(file," /= 0.");
         else put(file,"THE SOLUTIONS WITH zz");
              put(file,k,1); put_line(file," /= 0 :");
              put(file,Length_Of(sols1),Head_Of(sols1).n,sols1);
        end if;
  end if;
end Filter_and_Split_Solutions;
