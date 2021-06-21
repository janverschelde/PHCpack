with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Symbol_Table;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_JacoMats;    use Standard_Complex_Laur_JacoMats;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Double_Linear_Laurent_Solvers;
with Double_Lseries_Polynomials;
with Double_Lseries_Newton_Steps;

package body Main_Laurent_Series_Newton is

  procedure Run_Regular_Newton
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys; 
                x : in Standard_Complex_Vectors.Vector;
                tdx,deg : in integer32; verbose : in boolean ) is

    use Double_Lseries_Polynomials;

    neq : constant integer32 := p'last;
    dim : constant integer32 := neq + 1;
    nvr : constant integer32 := neq; -- number of variables minus t
    tv : constant Table_Vector(neq) := Make_Table_Vector(p,dim,nvr,tdx,deg);
    jp : constant Jaco_Mat(1..neq,1..dim) := Create(p);
    tva : constant Table_Vector_Array(1..neq)
        := Make_Table_Vector_Array(jp,tdx,deg);
    xlead : Standard_Integer_Vectors.Vector(1..nvr);
    xcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    dxnrm,pxnrm : double_float;
    numit : integer32 := 0;

  begin
    Double_Lseries_Newton_Steps.Make_Series(x,deg,xlead,xcffs);
    Double_Lseries_Newton_Steps.Set_Leading_Exponents(xlead);
    Double_Lseries_Newton_Steps.Run_Newton_Steps
      (deg,tv,tva,xlead,xcffs,dxnrm,pxnrm,numit,verbose=>verbose);
    new_line(file);
    Double_Linear_Laurent_Solvers.Write(file,xlead,xcffs,"x");
    put(file,"Maximum |dx| :"); put(file,dxnrm,3);
    put(file,"  residual :"); put(file,pxnrm,3); new_line(file);
  end Run_Regular_Newton;

  procedure Run_Singular_Newton
              ( file : in file_type;
                p,x : in Standard_Complex_Laur_Systems.Laur_Sys; 
                tdx,deg : in integer32; verbose : in boolean ) is

    use Double_Lseries_Polynomials;

    neq : constant integer32 := p'last;
    dim : constant integer32 := neq + 1;
    nvr : constant integer32 := neq; -- number of variables minus t
    tv : constant Table_Vector(neq) := Make_Table_Vector(p,dim,nvr,tdx,deg);
    jp : constant Jaco_Mat(1..neq,1..dim) := Create(p);
    tva : constant Table_Vector_Array(1..neq)
        := Make_Table_Vector_Array(jp,tdx,deg);
    xlead : Standard_Integer_Vectors.Vector(1..nvr);
    xcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    dxnrm,pxnrm : double_float;
    numit : integer32 := 0;

  begin
    Double_Lseries_Newton_Steps.Make_Series(x,deg,xlead,xcffs);
    Double_Lseries_Newton_Steps.Set_Leading_Exponents(xlead);
    Double_Lseries_Newton_Steps.Run_Newton_Steps
      (deg,tv,tva,xlead,xcffs,dxnrm,pxnrm,numit,verbose=>verbose);
    new_line(file);
    Double_Linear_Laurent_Solvers.Write(file,xlead,xcffs,"x");
    put(file,"Maximum |dx| :"); put(file,dxnrm,3);
    put(file,"  residual :"); put(file,pxnrm,3); new_line(file);
  end Run_Singular_Newton;

  procedure Start_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Laur_Sys;
    nq,nv,idx,deg : integer32 := 0;
    sols,ptr : Solution_List;
    ls : Link_to_Solution;
    verbose : boolean;
    ans : character;

  begin
    if vrb > 0 then
      put_line("-> in main_laurent_series_newton.Start_at_Constant");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      Standard_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      Standard_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    put(outfile,natural32(nq),natural32(nv),lp.all);
    idx := Double_Lseries_Polynomials.tsymbol_index;
    if idx = 0 then
      new_line;
      put_line("Warning: no parameter t found ... will exit now.");
    else
      new_line;
      put("The index of the parameter : "); put(idx,1); new_line;
      put("Give the truncation degree : "); get(deg);
      new_line;
      put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
      verbose := (ans = 'y');
      ptr := sols;
      while not Is_Null(ptr) loop
        ls := Head_Of(ptr);
        Run_Regular_Newton(outfile,lp.all,ls.v,idx,deg,verbose);
        ptr := Tail_Of(ptr);
        exit when Is_Null(ptr);
        new_line;
        put("Continue with the next solution ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end loop;
    end if;
  end Start_at_Constant;

  procedure Start_at_Series
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;

    infile,outfile : file_type;
    lp,lsol : Link_to_Laur_Sys;
    neq,nvr,dim,tdx,deg : integer32 := 0;
    ans : character;
    verbose : boolean;

  begin
    if vrb > 0 then
      put_line("-> in main_laurent_series_newton.Start_at_Series");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading a Laurent polynomial system ...");
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
      close(infile);
    end if;
    neq := lp'last;
    nvr := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Read "); put(neq,1); put(" polynomials in ");
    put(nvr,1); put_line(" variables ...");
    tdx := Double_Lseries_Polynomials.tsymbol_Index;
    put("-> index of t : "); put(tdx,1); new_line;
    if tdx = 0 then
      put_line("The polynomials must contain t as a variable!");
    else
      if outfilename = "" then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(outfile);
      else
        Create_Output_File(outfile,outfilename);
      end if;
      new_line;
      put_line("Reading initial terms of a series ...");
      Symbol_Table.Clear;
      get(lsol);
      dim := lsol'last; 
      new_line;
      put("Read "); put(dim,1); put_line(" polynomials ...");
      deg := Degree(lsol(lsol'first));
      put("The degree of the first polynomial : "); put(deg,1); new_line;
      new_line;
      put("Give the truncation degree : "); get(deg);
      new_line;
      put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
      verbose := (ans = 'y');
      Run_Singular_Newton(outfile,lp.all,lsol.all,tdx,deg,verbose);
    end if;
  end Start_at_Series;

  procedure Run_Laurent_Series_Newton
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    ans : character;

  begin
    if vrb > 0 then
      put("-> in main_laurent_series_newton.");
      put_line("Run_Laurent_Series_Newton ...");
    end if;
    put("Start Newton's method at a constant term ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Start_at_Constant(infilename,outfilename,vrb-1);
     else Start_at_Series(infilename,outfilename,vrb-1);
    end if;
  end Run_Laurent_Series_Newton;

end Main_Laurent_Series_Newton;
