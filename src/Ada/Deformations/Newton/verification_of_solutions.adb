with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
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
     then put_line(file,"yes");
     else put_line(file,"no");
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

  procedure Verify_Solutions_of_Polynomial_System
              ( p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector ) is

    loss : integer32;

  begin
    for i in z'range loop
      Newton_Steps_on_Polynomial_System
        (p,z(i),integer32(wanted),maxprec,maxitr,loss,
         err(integer32(i)),rco(integer32(i)),res(integer32(i)));
    end loop;
  end Verify_Solutions_of_Polynomial_System;

  procedure Verify_Solutions_of_Laurent_Polynomials
              ( p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector ) is

    loss : integer32;

  begin
    for i in z'range loop
      Newton_Steps_on_Laurent_Polynomials
        (p,z(i),integer32(wanted),maxprec,maxitr,loss,
         err(integer32(i)),rco(integer32(i)),res(integer32(i)));
    end loop;
  end Verify_Solutions_of_Laurent_Polynomials;

  procedure Verify_Solutions_of_Polynomial_System
              ( file : in file_type;
                p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector ) is

    loss : integer32;

  begin
    for i in z'range loop
      put(file,"solution "); put(file,natural32(i),1); put_line(file," :");
      Newton_Steps_on_Polynomial_System
        (file,p,z(i),integer32(wanted),maxprec,maxitr,loss,
         err(integer32(i)),rco(integer32(i)),res(integer32(i)));
      put_line(file,z(i).all);
    end loop;
  end Verify_Solutions_of_Polynomial_System;

  procedure Verify_Solutions_of_Laurent_Polynomials
              ( file : in file_type;
                p : in Array_of_Strings; z : in out Array_of_Strings;
                wanted,maxitr,maxprec : in natural32;
                err,rco,res : out Standard_Floating_Vectors.Vector ) is

    loss : integer32;

  begin
    for i in z'range loop
     -- declare
     -- begin
      put(file,"solution "); put(file,natural32(i),1); put_line(file," :");
      Newton_Steps_on_Laurent_Polynomials
        (file,p,z(i),integer32(wanted),maxprec,maxitr,loss,
         err(integer32(i)),rco(integer32(i)),res(integer32(i)));
      put_line(file,z(i).all);
     -- exception
     --   when others => 
     --     put("Exception at solution "); put(natural32(i),1); new_line; raise;
     -- end;
    end loop;
  end Verify_Solutions_of_Laurent_Polynomials;

  function to_Solutions
              ( z : Array_of_Strings;
                err,rco,res : Standard_Floating_Vectors.Vector )
              return Multprec_Complex_Solutions.Solution_List is

    result,last : Multprec_Complex_Solutions.Solution_List;

  begin
    for i in z'range loop
      declare
        sz : constant Multprec_Complex_Vectors.Vector
           := Multprec_Complex_Vector_Strings.Parse(z(i).all);
        dim : constant integer32 := sz'last;
        sol : Multprec_Complex_Solutions.Solution(dim); 
      begin
        sol.t := Multprec_Complex_Numbers.Create(integer32(1));
        sol.m := 1;
        for k in 1..dim loop
          sol.v(k) := sz(k);
        end loop;
        sol.err := Multprec_Floating_Numbers.create(err(integer32(i)));
        sol.rco := Multprec_Floating_Numbers.create(rco(integer32(i)));
        sol.res := Multprec_Floating_Numbers.create(res(integer32(i)));
        Multprec_Complex_Solutions.Append(result,last,sol);
      end;
    end loop;
    return result;
  end to_Solutions;

  procedure Verify_Solutions_of_Polynomial_System
              ( p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 ) is

    len : constant natural32 := Multprec_Complex_Solutions.Length_Of(sols);
    strsols : Array_of_Strings(1..integer(len)) := to_strings(sols);
    err,rco,res : Standard_Floating_Vectors.Vector(1..integer32(len));

  begin
    Verify_Solutions_of_Polynomial_System
      (p,strsols,wanted,maxitr,maxprc,err,rco,res);
    Multprec_Complex_Solutions.Clear(sols);
    sols := to_Solutions(strsols,err,rco,res);
  end Verify_Solutions_of_Polynomial_System;

  procedure Verify_Solutions_of_Laurent_Polynomials
              ( p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 ) is

    len : constant natural32 := Multprec_Complex_Solutions.Length_Of(sols);
    strsols : Array_of_Strings(1..integer(len)) := to_strings(sols);
    err,rco,res : Standard_Floating_Vectors.Vector(1..integer32(len));

  begin
    Verify_Solutions_of_Laurent_Polynomials
      (p,strsols,wanted,maxitr,maxprc,err,rco,res);
    Multprec_Complex_Solutions.Clear(sols);
    sols := to_Solutions(strsols,err,rco,res);
  end Verify_Solutions_of_Laurent_Polynomials;

  procedure Verify_Solutions_of_Polynomial_System
              ( file : in file_type; p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 ) is

    len : constant natural32 := Multprec_Complex_Solutions.Length_Of(sols);
    strsols : Array_of_Strings(1..integer(len)) := to_strings(sols);
    err,rco,res : Standard_Floating_Vectors.Vector(1..integer32(len));

  begin
    Verify_Solutions_of_Polynomial_System
      (file,p,strsols,wanted,maxitr,maxprc,err,rco,res);
    Multprec_Complex_Solutions.Clear(sols);
    sols := to_Solutions(strsols,err,rco,res);
  end Verify_Solutions_of_Polynomial_System;

  procedure Verify_Solutions_of_Laurent_Polynomials
              ( file : in file_type; p : in Array_of_Strings;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                wanted,maxitr,maxprc : in natural32 ) is

    len : constant natural32 := Multprec_Complex_Solutions.Length_Of(sols);
    strsols : Array_of_Strings(1..integer(len)) := to_strings(sols);
    err,rco,res : Standard_Floating_Vectors.Vector(1..integer32(len));

  begin
    Verify_Solutions_of_Laurent_Polynomials
      (file,p,strsols,wanted,maxitr,maxprc,err,rco,res);
    Multprec_Complex_Solutions.Clear(sols);
    sols := to_Solutions(strsols,err,rco,res);
 -- exception
 --   when others => 
 --     put_line("Exception in Verify_Solutions_of_Laurent_Polynomials");
 --     raise;
  end Verify_Solutions_of_Laurent_Polynomials;

end Verification_of_Solutions;
