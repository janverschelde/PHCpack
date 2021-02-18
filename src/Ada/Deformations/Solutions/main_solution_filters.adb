with Communications_with_User;          use Communications_with_User;
with Numbers_io;                        use Numbers_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Symbol_Table,Symbol_Table_io;      use Symbol_Table;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;
with Standard_Select_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with DoblDobl_Solution_Filters;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with QuadDobl_Solution_Filters;

package body Main_Solution_Filters is

  procedure Write_Symbols is

    n : constant natural32 := Symbol_Table.Number;

  begin
    put("There are "); put(n,1); put(" symbols :");
    for i in 1..n loop
      put(" "); Symbol_Table_io.put(Symbol_Table.get(i));
    end loop;
    new_line;
  end Write_Symbols;

  function Prompt_Symbol return natural32 is

    res : natural32 := 0;
    name : Symbol;

  begin
    loop
      put("Give the name of a variable : ");
      name := (name'range => ' ');
      Symbol_Table_io.get(name);
      res := Symbol_Table.Get(name);
      exit when (res > 0);
      put("The name ");
      Symbol_Table_io.put(name);
      put_line(" does not occur in the symbol table.");
      Write_Symbols;
      put_line("Please try again ...");
    end loop;
    return res;
  end Prompt_Symbol;

  procedure Read_Double_Double ( f : out double_double ) is
  begin
    f := create(0.0);
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a double double, please try again : ");
      Read_Double_Double(f);
  end Read_Double_Double;

  procedure Read_Quad_Double ( f : out quad_double ) is
  begin
    f := create(0.0);
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a double double, please try again : ");
      Read_Quad_Double(f);
  end Read_Quad_Double;

  procedure Read_Double_Complex
              ( c : out Standard_Complex_Numbers.Complex_Number ) is

    re,im : double_float;

  begin
    put("  give real part : "); Read_Double_Float(re);
    put("  give imaginary part : "); Read_Double_Float(im);
    c := Standard_Complex_Numbers.Create(re,im);
  end Read_Double_Complex;

  procedure Read_DoblDobl_Complex
              ( c : out DoblDobl_Complex_Numbers.Complex_Number ) is

    st_re,st_im : double_float;
    dd_re,dd_im : double_double;

  begin
    put("  give real part : "); Read_Double_Float(st_re);
    dd_re := create(st_re);
    put("  give imaginary part : "); Read_Double_Float(st_im);
    dd_im := create(st_im);
    c := DoblDobl_Complex_Numbers.Create(dd_re,dd_im);
  end Read_DoblDobl_Complex;

  procedure Read_QuadDobl_Complex
              ( c : out QuadDobl_Complex_Numbers.Complex_Number ) is

    st_re,st_im : double_float;
    qd_re,qd_im : quad_double;

  begin
    put("  give real part : "); Read_Double_Float(st_re);
    qd_re := create(st_re);
    put("  give imaginary part : "); Read_Double_Float(st_im);
    qd_im := create(st_im);
    c := QuadDobl_Complex_Numbers.Create(qd_re,qd_im);
  end Read_QuadDobl_Complex;

  function Show_Menu_and_Prompt_Answer return character is

    ans : character;

  begin
    new_line;
    put_line("MENU for filtering solution lists subject to criteria :");
    put_line("  1. reached a value of the continuation parameter, or");
    put_line("  2. target value of continuation parameter not reached;");
    put_line("  3. vanishing : error and residual smaller than tolerance, or");
    put_line("  4. spurious : not vanishing;");
    put_line("  5. regular : inverse condition number > tolerance, or");
    put_line("  6. singular : not regular;");
    put_line("  7. zero component : a component smaller than tolerance, or");
    put_line("  8. free component : a component that is nonzero;");
    put_line("  9. real solutions have imaginary part < tolerance;");
    put_line("  A. given their indices, select some solutions.");
    put("Give a number between 1 and 9, or A, to make your choice : ");
    Ask_Alternative(ans,"123456789A"); 
    return ans;
  end Show_Menu_and_Prompt_Answer;

  procedure Main_Solution_Filter
               ( file : in file_type;
                 sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Standard_Solution_Filters;

    ans : constant character := Show_Menu_and_Prompt_Answer;
    tol : double_float;
    target : Complex_Number;
    filtsols : Solution_List;
    k : natural32;

  begin
    new_line;
    case ans is
      when '1' | '2'
        => put_line("Give a complex target value : ");
           Read_Double_Complex(target);
           put("Give tolerance to test equality : ");
           Read_Double_Float(tol);
           if ans = '1'
            then filtsols := On_Target_Filter(sols,target,tol);
            else filtsols := Off_Target_Filter(sols,target,tol);
           end if;
      when '3' | '4'
        => put("Give tolerance for threshold on residual : ");
           Read_Double_Float(tol);
           if ans = '3'
            then filtsols := Vanishing_Filter(sols,tol);
            else filtsols := Spurious_Filter(sols,tol);
           end if;
      when '5' | '6'
        => put("Give tolerance on inverse condition number : ");
           Read_Double_Float(tol);
           if ans = '5'
            then filtsols := Regular_Filter(sols,tol);
            else filtsols := Singular_Filter(sols,tol);
           end if;
      when '7' | '8'
        => put("Give tolerance to test whether zero or not : ");
           Read_Double_Float(tol);
           Write_Symbols;
           k := Prompt_Symbol;
           if ans = '7'
            then filtsols := Zero_Component_Filter(sols,k,tol);
            else filtsols := Free_Component_Filter(sols,k,tol);
           end if;
      when '9'
        => put("Give tolerance on size of imaginary part : ");
           Read_Double_Float(tol);
           filtsols := Real_Filter(sols,tol);
      when 'A'
        => put("Give number of solutions to select : ");
           get(k);
           declare
             sel : Standard_Natural_Vectors.Vector(1..integer32(k));
           begin
             put("Give "); put(k,1); put(" positive increasing numbers : ");
             get(sel);
             filtsols := Select_Solutions(sols,sel);
           end;
      when others => null;
    end case;
    new_line;
    if Is_Null(filtsols) then
      put_line("The filtered list is empty");
    else
      put_line("The filtered list is on the output file.");
      put(file,Length_Of(filtsols),natural32(Head_Of(filtsols).n),filtsols);
    end if;
  end Main_Solution_Filter;

  procedure Main_Solution_Filter
               ( file : in file_type;
                 sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Filters;

    ans : constant character := Show_Menu_and_Prompt_Answer;
    tol : double_float;
    target : Complex_Number;
    filtsols : Solution_List;
    k : natural32;

  begin
    new_line;
    case ans is
      when '1' | '2'
        => put_line("Give a complex target value : ");
           Read_DoblDobl_Complex(target);
           put("Give tolerance to test equality : ");
           Read_Double_Float(tol);
           if ans = '1'
            then filtsols := On_Target_Filter(sols,target,tol);
            else filtsols := Off_Target_Filter(sols,target,tol);
           end if;
      when '3' | '4'
        => put("Give tolerance for threshold on residual : ");
           Read_Double_Float(tol);
           if ans = '3'
            then filtsols := Vanishing_Filter(sols,tol);
            else filtsols := Spurious_Filter(sols,tol);
           end if;
      when '5' | '6'
        => put("Give tolerance on inverse condition number : ");
           Read_Double_Float(tol);
           if ans = '5'
            then filtsols := Regular_Filter(sols,tol);
            else filtsols := Singular_Filter(sols,tol);
           end if;
      when '7' | '8'
        => put("Give tolerance to test whether zero or not : ");
           Read_Double_Float(tol);
           Write_Symbols;
           k := Prompt_Symbol;
           if ans = '7'
            then filtsols := Zero_Component_Filter(sols,k,tol);
            else filtsols := Free_Component_Filter(sols,k,tol);
           end if;
      when '9'
        => put("Give tolerance on size of imaginary part : ");
           Read_Double_Float(tol);
           filtsols := Real_Filter(sols,tol);
      when 'A'
        => put("Give number of solutions to select : ");
           get(k);
           declare
             sel : Standard_Natural_Vectors.Vector(1..integer32(k));
           begin
             put("Give "); put(k,1); put(" positive increasing numbers : ");
             get(sel);
             filtsols := Select_Solutions(sols,sel);
           end;
      when others => null;
    end case;
    new_line;
    if Is_Null(filtsols) then
      put_line("The filtered list is empty");
    else
      put_line("The filtered list is on the output file.");
      put(file,Length_Of(filtsols),natural32(Head_Of(filtsols).n),filtsols);
    end if;
  end Main_Solution_Filter;

  procedure Main_Solution_Filter
               ( file : in file_type;
                 sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Filters;

    ans : constant character := Show_Menu_and_Prompt_Answer;
    tol : double_float;
    target : Complex_Number;
    filtsols : Solution_List;
    k : natural32;

  begin
    new_line;
    case ans is
      when '1' | '2'
        => put_line("Give a complex target value : ");
           Read_QuadDobl_Complex(target);
           put("Give tolerance to test equality : ");
           Read_Double_Float(tol);
           if ans = '1'
            then filtsols := On_Target_Filter(sols,target,tol);
            else filtsols := Off_Target_Filter(sols,target,tol);
           end if;
      when '3' | '4'
        => put("Give tolerance for threshold on residual : ");
           Read_Double_Float(tol);
           if ans = '3'
            then filtsols := Vanishing_Filter(sols,tol);
            else filtsols := Spurious_Filter(sols,tol);
           end if;
      when '5' | '6'
        => put("Give tolerance on inverse condition number : ");
           Read_Double_Float(tol);
           if ans = '5'
            then filtsols := Regular_Filter(sols,tol);
            else filtsols := Singular_Filter(sols,tol);
           end if;
      when '7' | '8'
        => put("Give tolerance to test whether zero or not : ");
           Read_Double_Float(tol);
           Write_Symbols;
           k := Prompt_Symbol;
           if ans = '7'
            then filtsols := Zero_Component_Filter(sols,k,tol);
            else filtsols := Free_Component_Filter(sols,k,tol);
           end if;
      when '9'
        => put("Give tolerance on size of imaginary part : ");
           Read_Double_Float(tol);
           filtsols := Real_Filter(sols,tol);
      when 'A'
        => put("Give number of solutions to select : ");
           get(k);
           declare
             sel : Standard_Natural_Vectors.Vector(1..integer32(k));
           begin
             put("Give "); put(k,1); put(" positive increasing numbers : ");
             get(sel);
             filtsols := Select_Solutions(sols,sel);
           end;
      when others => null;
    end case;
    new_line;
    if Is_Null(filtsols) then
      put_line("The filtered list is empty");
    else
      put_line("The filtered list is on the output file.");
      put(file,Length_Of(filtsols),natural32(Head_Of(filtsols).n),filtsols);
    end if;
  end Main_Solution_Filter;

  procedure Standard_Solution_Filter
              ( infile,outfile : in file_type; len,dim : in natural32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Standard_Solution_Filters;

    ans : constant character := Show_Menu_and_Prompt_Answer;
    tol : double_float;
    target : Complex_Number;
    ls : Link_to_Solution;
    k,cnt : natural32;

  begin
    new_line;
    case ans is
      when '1' | '2'
        => put_line("Give a complex target value : ");
           Read_Double_Complex(target);
           put("Give tolerance to test equality : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '1'
            then On_Target_Filter(infile,outfile,len,dim,target,tol,cnt);
            else Off_Target_Filter(infile,outfile,len,dim,target,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '1'
	    then put_line(" on target solutions to file.");
	    else put_line(" off target solutions to file.");
           end if;
      when '3' | '4'
        => put("Give tolerance for threshold on residual : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '3'
            then Vanishing_Filter(infile,outfile,len,dim,tol,cnt);
            else Spurious_Filter(infile,outfile,len,dim,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '3'
            then put_line(" vanishing solutions to file.");
            else put_line(" spurious solutions to file.");
           end if;
      when '5' | '6'
        => put("Give tolerance on inverse condition number : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '5'
            then Regular_Filter(infile,outfile,len,dim,tol,cnt);
            else Singular_Filter(infile,outfile,len,dim,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '5'
            then put_line(" regular solutions to file.");
            else put_line(" spurious solutions to file.");
           end if;
      when '7' | '8'
        => put("Give tolerance to test whether zero or not : ");
           Read_Double_Float(tol);
           Read_Next(infile,dim,ls); -- for symbol table !
           Write_Symbols;
           k := Prompt_Symbol;
           new_line;
           put("Scanning the input file ... ");
           if ans = '7'
            then Zero_Component_Filter
                   (infile,outfile,len,dim,k,tol,ls,cnt);
            else Free_Component_Filter
                   (infile,outfile,len,dim,k,tol,ls,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '7'
            then put(" solutions with zero component ");
            else put(" solutions with free component ");
           end if;
           put(k,1); put_line(".");
      when '9'
        => put("Give tolerance on size of imaginary part : ");
           Read_Double_Float(tol);
           put("Scanning the input file ... ");
           Real_Filter(infile,outfile,len,dim,tol,cnt);
           put("wrote "); put(cnt,1);
           put_line(" real solutions to file.");
      when 'A'
        => put("Give number of solutions to select : ");
           get(k);
           declare
             sel : Standard_Natural_Vectors.Vector(1..integer32(k));
           begin
             put("Give "); put(k,1); put(" positive increasing numbers : ");
             get(sel);
             put("Scanning the input file ... ");
             Select_Solutions(infile,outfile,len,dim,sel,cnt);
             put("wrote "); put(cnt,1); put_line(" solutions to file.");
           end;
      when others => null;
    end case;
  end Standard_Solution_Filter;

  procedure DoblDobl_Solution_Filter
              ( infile,outfile : in file_type; len,dim : in natural32 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Filters;

    ans : constant character := Show_Menu_and_Prompt_Answer;
    tol : double_float;
    target : Complex_Number;
    ls : Link_to_Solution;
    k,cnt : natural32;

  begin
    new_line;
    case ans is
      when '1' | '2'
        => put_line("Give a complex target value : ");
           Read_DoblDobl_Complex(target);
           put("Give tolerance to test equality : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '1'
            then On_Target_Filter(infile,outfile,len,dim,target,tol,cnt);
            else Off_Target_Filter(infile,outfile,len,dim,target,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '1'
	    then put_line(" on target solutions to file.");
	    else put_line(" off target solutions to file.");
           end if;
      when '3' | '4'
        => put("Give tolerance for threshold on residual : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '3'
            then Vanishing_Filter(infile,outfile,len,dim,tol,cnt);
            else Spurious_Filter(infile,outfile,len,dim,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '3'
            then put_line(" vanishing solutions to file.");
            else put_line(" spurious solutions to file.");
           end if;
      when '5' | '6'
        => put("Give tolerance on inverse condition number : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '5'
            then Regular_Filter(infile,outfile,len,dim,tol,cnt);
            else Singular_Filter(infile,outfile,len,dim,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '5'
            then put_line(" regular solutions to file.");
            else put_line(" spurious solutions to file.");
           end if;
      when '7' | '8'
        => put("Give tolerance to test whether zero or not : ");
           Read_Double_Float(tol);
           Read_Next(infile,dim,ls); -- for symbol table !
           Write_Symbols;
           k := Prompt_Symbol;
           new_line;
           put("Scanning the input file ... ");
           if ans = '7'
            then Zero_Component_Filter
                   (infile,outfile,len,dim,k,tol,ls,cnt);
            else Free_Component_Filter
                   (infile,outfile,len,dim,k,tol,ls,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '7'
            then put(" solutions with zero component ");
            else put(" solutions with free component ");
           end if;
           put(k,1); put_line(".");
      when '9'
        => put("Give tolerance on size of imaginary part : ");
           Read_Double_Float(tol);
           put("Scanning the input file ... ");
           Real_Filter(infile,outfile,len,dim,tol,cnt);
           put("wrote "); put(cnt,1);
           put_line(" real solutions to file.");
      when 'A'
        => put("Give number of solutions to select : ");
           get(k);
           declare
             sel : Standard_Natural_Vectors.Vector(1..integer32(k));
           begin
             put("Give "); put(k,1); put(" positive increasing numbers : ");
             get(sel);
             put("Scanning the input file ... ");
             Select_Solutions(infile,outfile,len,dim,sel,cnt);
             put("wrote "); put(cnt,1); put_line(" solutions to file.");
           end;
      when others => null;
    end case;
  end DoblDobl_Solution_Filter;

  procedure QuadDobl_Solution_Filter
              ( infile,outfile : in file_type; len,dim : in natural32 ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Filters;

    ans : constant character := Show_Menu_and_Prompt_Answer;
    tol : double_float;
    target : Complex_Number;
    ls : Link_to_Solution;
    k,cnt : natural32;

  begin
    new_line;
    case ans is
      when '1' | '2'
        => put_line("Give a complex target value : ");
           Read_QuadDobl_Complex(target);
           put("Give tolerance to test equality : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '1'
            then On_Target_Filter(infile,outfile,len,dim,target,tol,cnt);
            else Off_Target_Filter(infile,outfile,len,dim,target,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '1'
	    then put_line(" on target solutions to file.");
	    else put_line(" off target solutions to file.");
           end if;
      when '3' | '4'
        => put("Give tolerance for threshold on residual : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '3'
            then Vanishing_Filter(infile,outfile,len,dim,tol,cnt);
            else Spurious_Filter(infile,outfile,len,dim,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '3'
            then put_line(" vanishing solutions to file.");
            else put_line(" spurious solutions to file.");
           end if;
      when '5' | '6'
        => put("Give tolerance on inverse condition number : ");
           Read_Double_Float(tol);
           new_line;
           put("Scanning the input file ... ");
           if ans = '5'
            then Regular_Filter(infile,outfile,len,dim,tol,cnt);
            else Singular_Filter(infile,outfile,len,dim,tol,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '5'
            then put_line(" regular solutions to file.");
            else put_line(" spurious solutions to file.");
           end if;
      when '7' | '8'
        => put("Give tolerance to test whether zero or not : ");
           Read_Double_Float(tol);
           Read_Next(infile,dim,ls); -- for symbol table !
           Write_Symbols;
           k := Prompt_Symbol;
           new_line;
           put("Scanning the input file ... ");
           if ans = '7'
            then Zero_Component_Filter
                   (infile,outfile,len,dim,k,tol,ls,cnt);
            else Free_Component_Filter
                   (infile,outfile,len,dim,k,tol,ls,cnt);
           end if;
           put("wrote "); put(cnt,1);
           if ans = '7'
            then put(" solutions with zero component ");
            else put(" solutions with free component ");
           end if;
           put(k,1); put_line(".");
      when '9'
        => put("Give tolerance on size of imaginary part : ");
           Read_Double_Float(tol);
           put("Scanning the input file ... ");
           Real_Filter(infile,outfile,len,dim,tol,cnt);
           put("wrote "); put(cnt,1);
           put_line(" real solutions to file.");
      when 'A'
        => put("Give number of solutions to select : ");
           get(k);
           declare
             sel : Standard_Natural_Vectors.Vector(1..integer32(k));
           begin
             put("Give "); put(k,1); put(" positive increasing numbers : ");
             get(sel);
             put("Scanning the input file ... ");
             Select_Solutions(infile,outfile,len,dim,sel,cnt);
             put("wrote "); put(cnt,1); put_line(" solutions to file.");
           end;
      when others => null;
    end case;
  end QuadDobl_Solution_Filter;

  procedure Read_Dimensions
              ( infile : out file_type; len,dim : out natural32 ) is

    bannered,fail : boolean;

  begin
    loop
      new_line;
      put_line("Reading the name of the file for the solutions.");
      Read_Name_and_Open_File(infile);
      Standard_Select_Solutions.Prompt_to_Scan_Banner(infile,bannered,fail);
      exit when not fail;
      Close(infile);
      put("Something was wrong with the input file.");
      put_line("  Please try again...");
    end loop;
    Standard_Select_Solutions.Read_Dimensions(infile,len,dim,fail);
    put("Ready to process "); put(len,1);
    put(" solutions of dimension "); put(dim,1); put_line(" ... ");
  exception
    when others => close(infile); raise;
  end Read_Dimensions;

  procedure Main ( infilename,outfilename : in string ) is

    infile,outfile : file_type;
    len,dim : natural32;
    precision : character;

  begin
    new_line;
    put_line("Filtering solution lists subject to criteria.");
    declare
    begin
      Read_Dimensions(infile,len,dim);
    exception
      when others
        => put_line("Please try again...");
           Read_Dimensions(infile,len,dim);
    end;
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    precision := Communications_with_User.Prompt_for_Precision;
    case precision is
      when '0' => Standard_Solution_Filter(infile,outfile,len,dim);
      when '1' => DoblDobl_Solution_Filter(infile,outfile,len,dim);
      when '2' => QuadDobl_Solution_Filter(infile,outfile,len,dim);
      when others => null;
    end case;
  end Main;

end Main_Solution_Filters;
