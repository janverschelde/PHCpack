with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Numbers_io;                        use Numbers_io;
with File_Scanning,String_Splitters;    use File_Scanning,String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_to_Multprec_Convertors;   use Standard_to_Multprec_Convertors;
with Multprec_Complex_Poly_Systems;
--with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
--with Multprec_Complex_Solutions;
--with Multprec_Complex_Solutions_io;     use Multprec_Complex_Solutions_io;
with Standard_Deflate_Singularities;    use Standard_Deflate_Singularities;
--with Multprec_Deflate_Singularities;    use Multprec_Deflate_Singularities;
with Standard_Deflation_Methods;        use Standard_Deflation_Methods;
with DoblDobl_Deflation_Methods;
with QuadDobl_Deflation_Methods;
with Multprec_Deflation_Methods;        use Multprec_Deflation_Methods;
with Drivers_to_Deflate_Singularities;  use Drivers_to_Deflate_Singularities;

procedure ts_deflate is

-- DESCRIPTION :
--   Interactive testing of deflation to deal with singularities.

  function Display_Menu_and_Prompt_Answer return character is

    ans : character;

  begin
    new_line;
    put_line("MENU to test Newton's method with deflation :");
    put_line("  1. compute the deflated systems in standard arithmetic;");
    put_line("  2. compute the deflations in multi-precision arithmetic;");
    put_line("  3. run Newton with deflation in standard arithmetic;");
    put_line("  4. run Newton with deflation in multi-precision arithmetic;");
    put_line("  5. test algorithmic deflation on standard solution list;");
    put_line("  6. test algorithmic deflation on multprec solution list;");
    put_line("  7. run the main driver to apply Newton with deflation;");
    put_line("  8. strip multipliers after deflation for standard numbers;");
    put_line("  9. strip multipliers after deflation in multi-precision.");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to make your choice : ");
    Ask_Alternative(ans,"123456789");
    return ans;
  end Display_Menu_and_Prompt_Answer;

  procedure Standard_Strip_Multipliers
              ( file : file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    nq,nv : natural32 := 0;

  begin
    new_line;
    put("Give the number of original equations : "); get(nq);
    put("Give the number of original variables : "); get(nv);
    declare
      f : constant Poly_Sys(1..integer32(nq)) := Strip_Multipliers(p,nq,nv);
      s : Solution_List;
    begin
      put(file,nq,nv,f);
      if not Is_Null(sols) then
        s := Strip_Multipliers(sols,nv);
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(s),natural32(Head_Of(s).n),s);
      end if;
    end;
  end Standard_Strip_Multipliers;

 -- procedure Multprec_Strip_Multipliers
 --             ( file : file_type;
 --               p : in Multprec_Complex_Poly_Systems.Poly_Sys;
 --               sols : in Multprec_Complex_Solutions.Solution_List ) is

 --   use Multprec_Complex_Poly_Systems;
 --   use Multprec_Complex_Solutions;

 --   nq,nv : natural32 := 0;

 -- begin
 --   new_line;
 --   put("Give the number of original equations : "); get(nq);
 --   put("Give the number of original variables : "); get(nv);
 --   declare
 --     f : constant Poly_Sys(1..integer32(nq)) := Strip_Multipliers(p,nq,nv);
 --     s : Solution_List;
 --   begin
 --     put(file,f);
 --     if not Is_Null(sols) then
 --       s := Strip_Multipliers(sols,nv);
 --       new_line(file);
 --       put_line(file,"THE SOLUTIONS :");
 --       put(file,Length_Of(s),natural32(Head_Of(s).n),s);
 --     end if;
 --   end;
 -- end Multprec_Strip_Multipliers;

  procedure Standard_Algorithmic_Deflation
              ( file : in file_type; name : in string; 
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Driver to the Algorithmic_Deflation_and_Clustering procedure.

    use Standard_Complex_Solutions;

    tolerr : constant double_float := 1.0E-12;
    tolres : constant double_float := 1.0E-12;
    tolrnk : constant double_float := 1.0E-06;
    maxitr : constant natural32 := 3;
    maxdef : natural32 := 0;
    output : boolean;
    ans : character;

  begin
    new_line;
    put("Give maximal number of deflations : "); get(maxdef);
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    new_line;
    Algorithmic_Deflation_and_Clustering
      (file,name,p,sols,output,maxitr,maxdef,tolerr,tolres,tolrnk);
    new_line(file);
    put_line(file,"THE SOLUTIONS after deflation :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Standard_Algorithmic_Deflation;

  procedure Multprec_Algorithmic_Deflation is
  begin
    put_line("in progress...");
  end Multprec_Algorithmic_Deflation;

  procedure Main is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    n : natural32 := 0;
    onfile,found : boolean;
    lp : Link_to_Poly_Sys;
    ans,choice : character;
    sols : Solution_List;
    tol : double_float := 1.0E-6;

  begin
    new_line;
    put_line("Testing operations in standard deflate singularities ...");
    choice := Display_Menu_and_Prompt_Answer;
    new_line;
    put("Is the polynomial system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      put_line("Reading the name of the file of the system.");
      Read_Name_and_Open_File(infile);
      get(infile,lp);
      onfile := true;
    else
      put("Give the dimension of the system : "); get(n);
      Symbol_Table.Init(n);
      lp := new Poly_Sys(1..integer32(n));
      put("Give "); put(n,1); put_line(" polynomials :");
      get(lp.all); skip_line;
      onfile := false;
    end if;
    new_line;
    put_line("Reading the name of the output file.");
    declare
      outfilename : constant string := Read_String;
    begin
      Create_Output_File(outfile,outfilename);
      if ((choice /= '8') and (choice /= '9')) then
        put(outfile,lp'last,1); new_line(outfile);
        put(outfile,lp.all); new_line(outfile);
      end if;
      if ((choice /= '1') and (choice /= '2')) then
        if onfile then
          Scan_and_Skip(infile,"SOLUTIONS",found);
          if found
           then get(infile,sols);
           else new_line;
                put_line("No solutions found on the input file.");
          end if;
          close(infile);
          if ((choice /= '7') and (choice /= '8')) then
            if not Is_Null(sols)
             then put_line(outfile,"The solutions on input :");
                  put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
            end if;
          end if;
        else
          new_line;
          put_line("No solutions found on the input file.");
        end if;
      end if;
      case choice is
        when '1'
          => new_line;
             Multiple_Standard_Deflations(outfile,lp.all);
        when '3'
          => Read_Tolerance(tol);
             Interactive_Symbolic_Deflation (outfile,lp.all,sols,tol);
        when '2' | '4'
          => declare
               dgts,size : natural32;
               mp : Multprec_Complex_Poly_Systems.Poly_Sys(lp'range);
             begin
               new_line;
               put("Give the number of decimal places : ");
               Read_Natural(dgts);
               size := Decimal_to_Size(dgts);
               mp := Convert(lp.all);
               Set_Size(mp,size);
               if ans = '2'
                then Multiple_Multprec_Deflations(outfile,mp,size);
                else Read_Tolerance(tol);
                     Interactive_Symbolic_Deflation
                       (outfile,mp,sols,size,tol);
               end if;
             end;
        when '5'
          => Standard_Algorithmic_Deflation(outfile,outfilename,lp.all,sols);
        when '6' => Multprec_Algorithmic_Deflation;
        when '8' => Standard_Strip_Multipliers(outfile,lp.all,sols);
       -- when '9' => Multprec_Strip_Multipliers(outfile,lp.all,sols);
        when others
          => new_line;
             put_line("calling the main driver ...");
             Deflate_Singularities(outfile,outfilename,lp.all,sols);
      end case;
    end;
  end Main;

begin
  Main;
end ts_deflate;
