with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Complex_Vector_Strings;
with Varbprec_Complex_Newton_Steps;      use Varbprec_Complex_Newton_Steps;

package body Verification_of_Solutions is

  procedure Set_Default_Parameters
              ( wanted,maxitr,maxprc : out natural32;
                verbose : out boolean ) is
  begin
    wanted := 8;
    maxitr := 4;
    maxprc := 256;
    verbose := true;
  end Set_Default_Parameters;

  procedure Write_Parameters
              ( file : in file_type;
                wanted,maxitr,maxprc : in natural32;
                verbose : in boolean ) is
  begin
    put_line(file,"Parameters for variable precision Newton steps :");
    put(file,"  1. number of accurate decimal places wanted : ");
    put(file,wanted); new_line(file);
    put(file,"  2. maximum number of Newton steps in the sequence : ");
    put(file,maxitr); new_line(file);
    put(file,"  3. maximum number of decimal places in the precision : ");
    put(file,maxprc); new_line(file);
    put(file,"  4. intermediate output and diagnostics during steps : ");
    if verbose
     then put_line("yes");
     else put_line("no");
    end if;
  end Write_Parameters;

  procedure Menu_to_set_Parameters
              ( wanted,maxitr,maxprc : out natural32;
                verbose : out boolean ) is

    ans : character;

  begin
    Set_Default_Parameters(wanted,maxitr,maxprc,verbose);
    loop
      new_line;
      Write_Parameters(standard_output,wanted,maxitr,maxprc,verbose);
      put("Type 1, 2, 3, or 4 to change a value, or 0 to exit : ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      case ans is
        when '1' =>
          put("Give number of wanted accurate decimal places : ");
          get(wanted);
        when '2' =>
          put("Give maximum number of Newton steps in the sequence : ");
          get(maxitr);
        when '3' =>
          put("Give maximum number of decimal places in the precision : ");
          get(maxprc);
        when '4' =>
          put("Intermediate output and diagnostics wanted in steps ? (y/n) ");
          Ask_Yes_or_No(ans);
          verbose := (ans = 'y');
        when others => null;
      end case;
    end loop;
  end Menu_to_Set_Parameters;

  function to_strings ( sols : Multprec_Complex_Solutions.Solution_List )
                      return Array_of_Strings is

    use Multprec_Complex_Solutions;
    len : constant natural32 := Length_Of(sols);
    res : Array_of_Strings(1..integer(len));
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      declare
        sv : constant string
           := Multprec_Complex_Vector_Strings.Write(ls.v); 
      begin
        res(i) := new string'(sv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end to_strings;

  procedure Verify ( p : in Array_of_Strings; z : in out Array_of_Strings;
                     wanted,maxitr,maxprec : in natural32;
                     err,rco,res : out Standard_Floating_Vectors.Vector ) is

    loss : integer32;

  begin
    for i in z'range loop
      Newton_Steps_to_Wanted_Accuracy
        (p,z(i),integer32(wanted),maxprec,maxitr,loss,
         err(integer32(i)),rco(integer32(i)),res(integer32(i)));
    end loop;
  end Verify;

  procedure Verify ( file : in file_type;
                     p : in Array_of_Strings; z : in out Array_of_Strings;
                     wanted,maxitr,maxprec : in natural32;
                     err,rco,res : out Standard_Floating_Vectors.Vector ) is

    loss : integer32;

  begin
    for i in z'range loop
      put(file,"solution "); put(file,natural32(i),1); put_line(file," :");
      Newton_Steps_to_Wanted_Accuracy
        (file,p,z(i),integer32(wanted),maxprec,maxitr,loss,
         err(integer32(i)),rco(integer32(i)),res(integer32(i)));
      put_line(file,z(i).all);
    end loop;
  end Verify;

end Verification_of_Solutions;
