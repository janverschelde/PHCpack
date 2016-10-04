with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Random_Polynomials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Random_Polynomials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Random_Polynomials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Multprec_Random_Polynomials;

package body Random_Polynomial_Systems is

  procedure Save_System_to_File
              ( s : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Save_System_to_File
              ( s : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Save_System_to_File
              ( s : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Save_System_to_File
              ( s : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;
    ans : character;

  begin
    put("Do you want to write the system to file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Read_Name_and_Create_File(file);
      put_line(file,s);
    else
      put_line(s);
    end if;
  end Save_System_to_File;

  procedure Standard_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    s : Standard_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : Standard_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := Standard_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := Standard_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
    lp := new Standard_Complex_Poly_Systems.Poly_Sys'(s);
  end Standard_Generate_and_Show;

  procedure DoblDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    s : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : DoblDobl_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := DoblDobl_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := DoblDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
    lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(s);
  end DoblDobl_Generate_and_Show;

  procedure QuadDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    s : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : QuadDobl_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := QuadDobl_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := QuadDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
    lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(s);
  end QuadDobl_Generate_and_Show;

  procedure Multprec_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    s : Multprec_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    for i in 1..e loop
      declare
        p : Multprec_Complex_Polynomials.Poly;
      begin
        if m = 0
         then p := Multprec_Random_Polynomials.Random_Dense_Poly(n,d,c);
         else p := Multprec_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
        end if;
        if ans = 'y' then
          new_line;
          put("-> p = ");
          if c /= 1 and c /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        s(i) := p;
      end;
    end loop;
    Save_System_to_File(s);
    lp := new Multprec_Complex_Poly_Systems.Poly_Sys'(s);
  end Multprec_Generate_and_Show;

end Random_Polynomial_Systems;
