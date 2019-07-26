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

  function Standard_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..neq);

    use Standard_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : Standard_Complex_Polynomials.Poly;
      begin
        if mct = 0
         then p := Random_Dense_Poly(nvr,deg,ctp);
         else p := Random_Sparse_Poly(nvr,deg,mct,ctp);
        end if;
        if verbose then
          new_line;
          put("-> p = ");
          if ctp /= 1 and ctp /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        res(i) := p;
      end;
    end loop;
    return res;
  end Standard_Generate;

  function DoblDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..neq);

    use DoblDobl_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : DoblDobl_Complex_Polynomials.Poly;
      begin
        if mct = 0
         then p := Random_Dense_Poly(nvr,deg,ctp);
         else p := Random_Sparse_Poly(nvr,deg,mct,ctp);
        end if;
        if verbose then
          new_line;
          put("-> p = ");
          if ctp /= 1 and ctp /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        res(i) := p;
      end;
    end loop;
    return res;
  end DoblDobl_Generate;

  function QuadDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..neq);

    use QuadDobl_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : QuadDobl_Complex_Polynomials.Poly;
      begin
        if mct = 0
         then p := Random_Dense_Poly(nvr,deg,ctp);
         else p := Random_Sparse_Poly(nvr,deg,mct,ctp);
        end if;
        if verbose then
          new_line;
          put("-> p = ");
          if ctp /= 1 and ctp /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        res(i) := p;
      end;
    end loop;
    return res;
  end QuadDobl_Generate;

  function Multprec_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(1..neq);

    use Multprec_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : Multprec_Complex_Polynomials.Poly;
      begin
        if mct = 0
         then p := Random_Dense_Poly(nvr,deg,ctp);
         else p := Random_Sparse_Poly(nvr,deg,mct,ctp);
        end if;
        if verbose then
          new_line;
          put("-> p = ");
          if ctp /= 1 and ctp /= 2
           then put_line(p);
           else put(p); new_line;
          end if;
        end if;
        res(i) := p;
      end;
    end loop;
    return res;
  end Multprec_Generate;

  procedure Standard_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : Standard_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := Standard_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new Standard_Complex_Poly_Systems.Poly_Sys'(sys);
  end Standard_Generate_and_Show;

  procedure DoblDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := DoblDobl_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(sys);
  end DoblDobl_Generate_and_Show;

  procedure QuadDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := QuadDobl_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(sys);
  end QuadDobl_Generate_and_Show;

  procedure Multprec_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : Multprec_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := Multprec_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new Multprec_Complex_Poly_Systems.Poly_Sys'(sys);
  end Multprec_Generate_and_Show;

end Random_Polynomial_Systems;
