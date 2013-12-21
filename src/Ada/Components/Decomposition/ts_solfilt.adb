with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Numbers_io;                        use Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;         use Standard_Solution_Filters;

procedure ts_solfilt is

-- DESCRIPTION :
--   Calls the driver to filter solutions.

  file : file_type;
  sols : Solution_List;

  procedure Driver_for_Solution_Filters
               ( file : in file_type; sols : in Solution_List );

  -- DESCRIPTION :
  --   Creates new solution lists from filtering a given list subject
  --   to the following exclusive criteria :
  --     1) vanishing or spurious;
  --     2) zero or free k-th component;
  --     3) regular or singular;
  --     4) reached a target value.
  --   The criteria are set up with respect to a given tolerance.

  procedure Read_Double_Complex ( c : out Complex_Number ) is

    re,im : double_float;

  begin
    put("  give real part : "); Read_Double_Float(re);
    put("  give imaginary part : "); Read_Double_Float(im);
    c := Create(re,im);
  end Read_Double_Complex;

  procedure Driver_for_Solution_Filters
               ( file : in file_type; sols : in Solution_List ) is

    ans : character;
    tol : double_float;
    target : Complex_Number;
    filtsols : Solution_List;
    k : positive;

  begin
    new_line;
    put_line("MENU for filtering solution lists subject to criterion.");
    put_line("  1. vanishing : error and residual smaller than tolerance");
    put_line("  2. spurious : not vanishing");
    put_line("  3. zero component : a component smaller than tolerance");
    put_line("  4. free component : a component that is nonzero");
    put_line("  5. regular : inverse condition number larger than tolerance");
    put_line("  6. singular : not regular");
    put_line("  7. reached a given target value.");
    put_line("  8. not reached a given target value.");
    put("Type your choice (number between 1 and 8) : ");
    Ask_Alternative(ans,"1234567"); 
    new_line;
    case ans is
      when '1' | '2'
        => put("Give tolerance for threshold on residual : ");
           Read_Double_Float(tol);
           if ans = '1'
            then filtsols := Vanishing_Filter(sols,tol);
            else filtsols := Spurious_Filter(sols,tol);
           end if;
      when '3' | '4'
        => put("Give tolerance to test whether zero or not : ");
           Read_Double_Float(tol);
           put("Give the index to the component : "); Read_Positive(k);
           if ans = '3'
            then filtsols := Zero_Component_Filter(sols,k,tol);
            else filtsols := Free_Component_Filter(sols,k,tol);
           end if;
      when '5' | '6'
        => put("Give tolerance on inverse condition number : ");
           Read_Double_Float(tol);
           if ans = '5'
            then filtsols := Regular_Filter(sols,tol);
            else filtsols := Singular_Filter(sols,tol);
           end if;
      when '7' | '8'
        => put_line("Give a complex target value : ");
           Read_Double_Complex(target);
           put("Give tolerance to test equality : ");
           Read_Double_Float(tol);
           if ans = '7'
            then filtsols := On_Target_Filter(sols,target,tol);
            else filtsols := Off_Target_Filter(sols,target,tol);
           end if;
      when others => null;
    end case;
    if Is_Null(filtsols)
     then put_line("The filtered list is empty");
     else put_line("The filtered list is on the output file.");
          put(file,Length_Of(filtsols),Head_Of(filtsols).n,filtsols);
    end if;
  end Driver_for_Solution_Filters;

begin
  new_line;
  put_line("Filtering solution lists subject to criteria.");
  new_line;
  put_line("WARNING : PROBLEMS WITH SYMBOL TABLE");
  new_line;
  Read(sols);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  Driver_for_Solution_Filters(file,sols);
end ts_solfilt;
