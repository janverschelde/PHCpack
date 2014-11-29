with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;                   use String_Splitters;
with Strings_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Strings;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;
with Verification_of_Solutions;          use Verification_of_Solutions;

procedure ts_vmphom is

-- DESCRIPTION :
--   Development of variable precision homotopy evaluation.

  function Strip_Semicolon ( f : string ) return string is

  -- DESCRIPTION :
  --   Removes the semicolon from the last position in the string f.

  begin
    if f(f'last) = ';'
     then return f(f'first..f'last-1);
     else return f;
    end if;
  end Strip_Semicolon;

  function Strip_First_Plus ( f : string ) return string is

  -- DESCRIPTION :
  --   Removes the plus sign at the first position of the string f.

  begin
    if f(f'first) = '+'
     then return (" " & f(f'first+1..f'last));
     else return f;
    end if;
  end Strip_First_Plus;

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
             ( gamma,t : Complex_Number; k : natural32;
               f,g : Array_of_Strings ) return Array_of_Strings is

  -- DESCRIPTION :
  --   Returns the string gamma*(1-t)^k * g + t^k * f,
  --   evaluating the expressions gamma*(1-t)^k and t^k,
  --   but not the expressions for f and g.

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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a target and start system.

    sf_file,sg_file : file_type;
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
  end Main;

begin
  Main;
end ts_vmphom;
