with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;                   use String_Splitters;
with Strings_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Strings;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Poly_Strings;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;
with Verification_of_Solutions;          use Verification_of_Solutions;
with Random_Conditioned_Root_Problems;   use Random_Conditioned_Root_Problems;
with Random_Conditioned_Homotopies;      use Random_Conditioned_Homotopies;
with Varbprec_Complex_Newton_Steps;      use Varbprec_Complex_Newton_Steps;

procedure ts_vmphom is

-- DESCRIPTION :
--   Development of variable precision homotopy evaluation.

  function Write_Homotopy ( fcs,gcs,f,g : string ) return string is

  -- DESCRIPTION :
  --   Returns the string representation of a homotopy between g and f
  --   with corresponding coefficients gcs and fcs.

    res : constant string
        := gcs & "*(" & Strip_Semicolon(g) & ") + " 
         & fcs & "*(" & Strip_Semicolon(f) & ");";

  begin
    return res;
  end Write_Homotopy;

  function Write_Homotopy
             ( gamma,t : Standard_Complex_Numbers.Complex_Number;
               k : natural32;
               f,g : Array_of_Strings ) return Array_of_Strings is

  -- DESCRIPTION :
  --   Returns the string gamma*(1-t)^k * g + t^k * f,
  --   evaluating the expressions gamma*(1-t)^k and t^k,
  --   but not the expressions for f and g.

    use Standard_Complex_Numbers;

    res : Array_of_Strings(f'range);
    ftk : constant Complex_Number := t**integer(k);
    one_min_t : constant Complex_Number := Create(1.0) - t; 
    gtk : constant Complex_Number := gamma*(one_min_t**integer(k));
    gcs : constant string := Strings_and_Numbers.Convert(gtk);
    fcs : constant string := Strings_and_Numbers.Convert(ftk);

  begin
    put("The factor for g : "); put_line(gcs);
    put("The factor for f : "); put_line(fcs);
    for i in res'range loop
      declare
        hom : constant string 
            := Write_Homotopy(fcs,gcs,f(i).all,g(i).all);
      begin
        res(i) := new string'(hom);
      end;
    end loop;
    return res;
  end Write_Homotopy;

  procedure Make_Homotopy
              ( m : in natural32; target,start : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Makes a string representation of an artificial-parameter homotopy
  --   between target and start system, in m variables.

  -- ON ENTRY :
  --   m        the number of variables;
  --   target   string representation of the target system;
  --   start    string representation of the start system;
  --   sols     solutions of the start system.

  -- ON RETURN :
  --   sols     updated start solutions.

    use Standard_Complex_Numbers;

    gamma : constant Complex_Number := Random1;
    t : Complex_Number;
    k : natural32 := 0;
    wanted,maxitr,maxprc : natural32;
    verbose : boolean;
 
  begin
    new_line;
    put("A random gamma : "); put(gamma); new_line; 
    put("Give a complex value for t : "); get(t);
    put("Give a positive natural value for k : "); get(k);
    Menu_to_Set_Parameters(wanted,maxitr,maxprc,verbose);
    Symbol_Table.Init(m);
    declare
      hom : constant Array_of_Strings
          := Write_Homotopy(gamma,t,k,target,start);
      sys : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Complex_Poly_Strings.Parse(m,hom);
    begin
      put_line("The string representation of the homotopy : ");
      for i in hom'range loop
        put_line(hom(i).all);
      end loop;
      put_line("The parsed string : "); put(sys);
      if verbose then
        Verify_Solutions_of_Polynomial_System
          (standard_output,hom,sols,wanted,maxitr,maxprc);
      else
        Verify_Solutions_of_Polynomial_System(hom,sols,wanted,maxitr,maxprc);
      end if;
      new_line;
      put_line("THE SOLUTIONS :");
      put(standard_output,Multprec_Complex_Solutions.Length_Of(sols),m,sols);
    end;
  end Make_Homotopy;

  procedure Test_Homotopy_to_String is

  -- DESCRIPTION :
  --   Prompts the user for a target and start system
  --   and writes the string representation of the standard homotopy.

    sf_file : file_type;
    sf,sg : Link_to_Array_of_strings;
    nq,nv,len : natural32;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading the target polynomial system ...");
    Read_Name_and_Open_File(sf_file);
    get(sf_file,integer(nq),integer(nv),sf);
    close(sf_file);
    new_line;
    put_line("Reading the start polynomial system with solutions ...");
    Multprec_System_and_Solutions_io.get(nq,nv,sg,sols);
    new_line;
    len := Multprec_Complex_Solutions.Length_Of(sols);
    put("Read list of "); put(len,1); put_line(" solutions from file.");
    Make_Homotopy(nv,sf.all,sg.all,sols);
  end Test_Homotopy_to_String;

  procedure Write_Homotopy_to_File
              ( h : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and then writes the homotopy h
  --   to the file.

    file : file_type;

  begin
    new_line;
    Read_Name_and_Create_File(file);
    put_line(file,h);
    close(file);
  end Write_Homotopy_to_File;

  procedure Write_Homotopy_to_File
              ( h : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and then writes the homotopy h
  --   to the file.

    file : file_type;

  begin
    new_line;
    Read_Name_and_Create_File(file);
    put(file,h'last,1); new_line(file);
    put_line(file,h);
    close(file);
  end Write_Homotopy_to_File;

  procedure Standard_Random_Conditioned_Homotopy
              ( hom : in Array_of_Strings; root : in string ) is

  -- DESCRIPTION :
  --   Tests the evaluation of the homotopy given in hom by its string
  --   representation using standard double precision arithmetic.

    use Standard_Complex_Numbers;

    dim : constant integer32 := integer32(hom'last);
    nvr : constant natural32 := natural32(dim);
    htp,fzero,fmid,fone : Standard_Complex_Poly_Systems.Poly_Sys(1..dim);
    strfmid : Array_of_Strings(hom'range);
    t : Complex_Number;
    jfrco,fzrco : double_float;
    loss : integer32;
    ans : character;

  begin
    Symbol_Table.Init(nvr+1);
    htp := Standard_Complex_Poly_Strings.Parse(nvr+1,hom);
    put_line("The homotopy parsed into a polynomial system :");
    put_line(htp);
    new_line;
    put("Write homotopy to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Write_Homotopy_to_File(htp);
    end if;
    Symbol_Table.Clear; -- wipe out symbol table for confusion about t
    t := Create(0.0);
    fzero := Standard_Complex_Poly_SysFun.Eval(htp,t,1); -- t is first !
    put_line("The start system :"); put_line(fzero);
    t := Create(1.0);
    fone := Standard_Complex_Poly_SysFun.Eval(htp,t,1);  -- t is first !
    put_line("The target system :"); put_line(fone);
    t := Create(0.5);
    fmid := Standard_Complex_Poly_SysFun.Eval(htp,t,1);  -- t is first !
    put_line("The system in the middle :"); put_line(fmid);
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(dim));
    strfmid := Standard_Complex_Poly_Strings.Write(fmid);
    Standard_Estimate_Loss_for_Polynomial_System(strfmid,root,jfrco,fzrco,loss);
    put("condition of Jacobian matrix : "); put(jfrco,3); new_line;
    put("condition of polynomial evaluation : "); put(fzrco,3); new_line;
    put("-> estimated loss of decimal places : "); put(loss,1); new_line;
    Standard_Complex_Poly_Systems.Clear(htp);
    Standard_Complex_Poly_Systems.Clear(fzero);
    Standard_Complex_Poly_Systems.Clear(fone);
    Standard_Complex_Poly_Systems.Clear(fmid);
  end Standard_Random_Conditioned_Homotopy;

  procedure Multprec_Random_Conditioned_Homotopy
              ( hom : in Array_of_Strings; root : in string;
                size : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the evaluation of the homotopy given in hom by its string
  --   representation using arbitrary multiprecision arithmetic.

    use Multprec_Complex_Numbers;

    dim : constant integer32 := integer32(hom'last);
    nvr : constant natural32 := natural32(dim);
    htp,fzero,fmid,fone : Multprec_Complex_Poly_Systems.Poly_Sys(1..dim);
    strfmid : Array_of_Strings(hom'range);
    t : Complex_Number;
    jfrco,fzrco : Floating_Number;
    deci : constant natural32
         := Multprec_Floating_Numbers.Size_to_Decimal(size);
    loss : integer32;
    ans : character;
    fnb : Floating_Number;

  begin
    Symbol_Table.Init(nvr+1);
    htp := Multprec_Complex_Poly_Strings.Parse(nvr+1,size,hom);
    put_line("The homotopy parsed into a polynomial system :");
    put_line(htp);
    new_line;
    put("Write homotopy to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Write_Homotopy_to_File(htp);
    end if;
    Symbol_Table.Clear; -- wipe out symbol table for confusion about t
    fnb := Create(0.0); t := Create(fnb); Clear(fnb);
    fzero := Multprec_Complex_Poly_SysFun.Eval(htp,t,1); -- t is first !
    put_line("The start system :"); put_line(fzero);
    fnb := Create(1.0); t := Create(fnb); Clear(fnb);
    fone := Multprec_Complex_Poly_SysFun.Eval(htp,t,1);  -- t is first !
    put_line("The target system :"); put_line(fone);
    fnb := Create(0.5); t := Create(fnb); Clear(fnb);
    fmid := Multprec_Complex_Poly_SysFun.Eval(htp,t,1);  -- t is first !
    put_line("The system in the middle :"); put_line(fmid);
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(dim));
    strfmid := Multprec_Complex_Poly_Strings.Write(fmid);
    Multprec_Estimate_Loss_for_Polynomial_System
      (strfmid,root,deci,jfrco,fzrco,loss);
    put("condition of Jacobian matrix : "); put(jfrco,3); new_line;
    put("condition of polynomial evaluation : "); put(fzrco,3); new_line;
    put("-> estimated loss of decimal places : "); put(loss,1); new_line;
    Multprec_Complex_Poly_Systems.Clear(htp);
    Multprec_Complex_Poly_Systems.Clear(fzero);
    Multprec_Complex_Poly_Systems.Clear(fone);
    Multprec_Complex_Poly_Systems.Clear(fmid);
  end Multprec_Random_Conditioned_Homotopy;

  procedure Test_Random_Conditioned_Homotopy
              ( f : in Array_of_Strings; root : in string ) is

  -- DESCRIPTION :
  --   Puts the system f at the middle of a homotopy,
  --   with its zero store in the string root.

    hom : Array_of_Strings(f'range) := Conditioned_Homotopy(f);
    deci,size : natural32 := 0;

  begin
    put_line("The polynomials in the conditioned homotopy :");
    for i in hom'range loop
      put_line(hom(i).all);
    end loop;
    new_line;
    put("Give number of decimal places (0 for standard double precision) : ");
    get(deci);
    if deci = 0 then
      Standard_Random_Conditioned_Homotopy(hom,root);
    else
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Multprec_Random_Conditioned_Homotopy(hom,root,size);
    end if;
    String_Splitters.Clear(hom);
  end Test_Random_Conditioned_Homotopy;

  procedure Test_Random_Conditioned_Homotopy is

  -- DESCRIPTION :
  --   Sets up a random conditioned root problem and puts this problem
  --   at the middle of a homotopy.

    fhalf : Link_to_Array_of_Strings;
    froot : Link_to_String;

  begin
    Random_Conditioned_Root_Problem(fhalf,froot);
    put_line("The root : "); put_line(froot.all);
    Test_Random_Conditioned_Homotopy(fhalf.all,froot.all);
    String_Splitters.Clear(fhalf);
    String_Splitters.Clear(froot);
  end Test_Random_Conditioned_Homotopy;

  procedure Main is

  -- DESCRIPTION :
  --   Presents the menu to test random conditioned homotopies.

    ans : character;

  begin
    new_line;
    put_line("MENU to test random conditioned homotopies :");
    put_line("  1. write a standard linear homotopy to a string;");
    put_line("  2. place a random conditioned system in the middle.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Test_Homotopy_to_String;
     else Test_Random_Conditioned_Homotopy;
    end if;
  end Main;

begin
  Main;
end ts_vmphom;
