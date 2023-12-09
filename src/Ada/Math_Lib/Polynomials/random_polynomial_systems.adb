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
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Polynomials_io;    use TripDobl_Complex_Polynomials_io;
with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with TripDobl_Random_Polynomials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Random_Polynomials;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Polynomials_io;    use PentDobl_Complex_Polynomials_io;
with PentDobl_Complex_Poly_Systems_io;   use PentDobl_Complex_Poly_Systems_io;
with PentDobl_Random_Polynomials;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Polynomials_io;    use OctoDobl_Complex_Polynomials_io;
with OctoDobl_Complex_Poly_Systems_io;   use OctoDobl_Complex_Poly_Systems_io;
with OctoDobl_Random_Polynomials;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Polynomials_io;    use DecaDobl_Complex_Polynomials_io;
with DecaDobl_Complex_Poly_Systems_io;   use DecaDobl_Complex_Poly_Systems_io;
with DecaDobl_Random_Polynomials;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Polynomials_io;    use HexaDobl_Complex_Polynomials_io;
with HexaDobl_Complex_Poly_Systems_io;   use HexaDobl_Complex_Poly_Systems_io;
with HexaDobl_Random_Polynomials;
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
              ( s : in TripDobl_Complex_Poly_Systems.Poly_Sys ) is

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
              ( s : in PentDobl_Complex_Poly_Systems.Poly_Sys ) is

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
              ( s : in OctoDobl_Complex_Poly_Systems.Poly_Sys ) is

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
              ( s : in DecaDobl_Complex_Poly_Systems.Poly_Sys ) is

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
              ( s : in HexaDobl_Complex_Poly_Systems.Poly_Sys ) is

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

  function TripDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Poly_Sys(1..neq);

    use TripDobl_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : TripDobl_Complex_Polynomials.Poly;
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
  end TripDobl_Generate;

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

  function PentDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    res : PentDobl_Complex_Poly_Systems.Poly_Sys(1..neq);

    use PentDobl_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : PentDobl_Complex_Polynomials.Poly;
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
  end PentDobl_Generate;

  function OctoDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..neq);

    use OctoDobl_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : OctoDobl_Complex_Polynomials.Poly;
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
  end OctoDobl_Generate;

  function DecaDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(1..neq);

    use DecaDobl_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : DecaDobl_Complex_Polynomials.Poly;
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
  end DecaDobl_Generate;

  function HexaDobl_Generate
             ( nvr,deg,mct,ctp : natural32; neq : integer32;
               verbose : boolean := false )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(1..neq);

    use HexaDobl_Random_Polynomials;

  begin
    for i in 1..neq loop
      declare
        p : HexaDobl_Complex_Polynomials.Poly;
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
  end HexaDobl_Generate;

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

  procedure TripDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : TripDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := TripDobl_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new TripDobl_Complex_Poly_Systems.Poly_Sys'(sys);
  end TripDobl_Generate_and_Show;

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

  procedure PentDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : PentDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := PentDobl_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new PentDobl_Complex_Poly_Systems.Poly_Sys'(sys);
  end PentDobl_Generate_and_Show;

  procedure OctoDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := OctoDobl_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new OctoDobl_Complex_Poly_Systems.Poly_Sys'(sys);
  end OctoDobl_Generate_and_Show;

  procedure DecaDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : DecaDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := DecaDobl_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new DecaDobl_Complex_Poly_Systems.Poly_Sys'(sys);
  end DecaDobl_Generate_and_Show;

  procedure HexaDobl_Generate_and_Show
              ( n,d,m,c : in natural32; e : in integer32;
                lp : out HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    sys : HexaDobl_Complex_Poly_Systems.Poly_Sys(1..e);
    ans : character;
    vrb : boolean;

  begin
    new_line;
    put("Do you want intermediate output during generation ? (y/n) ");
    Ask_Yes_or_No(ans);
    vrb := (ans = 'y');
    sys := HexaDobl_Generate(n,d,m,c,e,vrb);
    Save_System_to_File(sys);
    lp := new HexaDobl_Complex_Poly_Systems.Poly_Sys'(sys);
  end HexaDobl_Generate_and_Show;

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
