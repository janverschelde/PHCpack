with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;        use QuadDobl_Complex_Numbers_cv;
with Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Numbers_io;                         use Numbers_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with DoblDobl_Complex_Poly_Strings;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Strings;
with DoblDobl_Polynomial_Convertors;     use DoblDobl_Polynomial_Convertors;
with QuadDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Strings;
with QuadDobl_Polynomial_Convertors;     use QuadDobl_Polynomial_Convertors;
with Multprec_Complex_Poly_Strings;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Standard_Homotopy;
with Standard_Coefficient_Homotopy;
with Standard_Laurent_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Laurent_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Laurent_Homotopy;
with Projective_Transformations;         use Projective_Transformations;
with Homogenization;                     use Homogenization;
with Multprec_Homotopy;

package body Drivers_for_Homotopy_Creation is

  procedure Info_on_Precision is

    i : array(1..5) of string(1..65);

  begin
    i(1):="  The  number  of  decimal places determines the precision during";
    i(2):="path following.  For any value less than or equal to 16, standard";
    i(3):="machine  arithmetic will be used.  The expensive nature of multi-";
    i(4):="precision arithmetic severely limits the  dimension  of  problems";
    i(5):="solvable within in acceptable time frame.                        ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Precision;

  procedure Info_on_Power is

    i : array(1..6) of string(1..65);

  begin
    i(1):="  The  homotopy  parameter  k  determines  the   power   of   the";
    i(2):="continuation  parameter  t.  Taking k>1 has as effect that at the";
    i(3):="beginning and at the end of the continuation, changes in t do not";
    i(4):="change  the homotopy as much as with a homotopy that is linear in";
    i(5):="t so that paths are better to follow.  The default value  k=2  is";
    i(6):="recommended.                                                     ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Power;

  procedure Info_on_Constant is

    i : array(1..3) of string(1..65);

  begin
    i(1):="  The homotopy parameter a ensures the regularity of the solution";
    i(2):="paths, i.e.: by choosing a random complex number for a, all paths";
    i(3):="are regular and do not diverge for t<1.                          ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Constant;

  procedure Info_on_Target is

    i : array(1..7) of string(1..65);

  begin
    i(1):="  The target value for the continuation parameter t is by default";
    i(2):="1.   To  create  stepping stones in the continuation stage, it is";
    i(3):="possible to let the continuation stop at t<1, for instance at t =";
    i(4):="0.9  or even at a complex value for t.  The solutions at t<1 will";
    i(5):="serve at start solutions to take up the continuation in  a  later";
    i(6):="stage.   In this stage, the same homotopy parameters k and a must";
    i(7):="be used.                                                         ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Target;

  procedure Info_on_Exit is

    i : constant string 
      :="By typing 0, the current settings are used in the homotopy.";

  begin
    put_line(i);
  end Info_on_Exit;

  procedure Info_on_Projective_Transformation is

    i : array(1..5) of string(1..65);

  begin
    i(1):="  A projective transformation of the homotopy and start solutions";
    i(2):="makes  the equations in the polynomials homogeneous and adds an a";
    i(3):="random hyperplane.   The  vectors  of  the  start  solutions  are";
    i(4):="extended  with an additional unknown.  For solutions at infinity,";
    i(5):="this additional unknown equals zero.                             ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Projective_Transformation;

  procedure Read_Complex ( x : in out Complex_Number ) is

  -- DESCRIPTION :
  --   User friendly reading of a complex number.

    re,im : double_float;

  begin
    put("  Give a real number for real part : ");      Read_Double_Float(re);
    put("  Give a real number for imaginary part : "); Read_Double_Float(im);
    x := Create(re,im);
  end Read_Complex;

  procedure Default_Homotopy_Settings
                ( d,k : out natural32; a,t : out Complex_Number ) is
  begin
    d := 16;
    k := 2;
    a := Random1;
    t := Create(1.0);
  end Default_Homotopy_Settings;

  procedure Default_Homotopy_Settings
                ( d,k : out natural32; a,t : out Complex_Number;
                  proj : out boolean ) is
  begin
    Default_Homotopy_Settings(d,k,a,t);
    proj := false;
  end Default_Homotopy_Settings;

  procedure Show_Banner is

  -- DESCRIPTION :
  --   Shows the banner explaining the parameters in the homotopy.

  begin
    new_line;
    put_line
       ("Homotopy is H(x,t) = a*(1-t)^k * Q(x) + t^k * P(x) = 0, t in [0,1],");
    put_line
       ("      with Q(x) = 0 a start system, and P(x) = 0 the target system.");
  end Show_Banner;

  function Show_Menu_and_Prompt_Choice 
              ( d,k : natural32; a,t : Complex_Number;
                proj : boolean ) return string is

  -- DESCRIPTION :
  --   Shows the menu with the current settings of the parameters
  --   and prompts the user to make choice, returned in a string.

    res : string(1..2) := "  ";

  begin
    new_line;
    put_line("MENU with settings for homotopy :");
    put("  number of decimal places : "); put(d,2);
    if d <= 16 then
      put(" (standard double precision)");
    elsif d <= 32 then
      put(" (double double precision)");
    elsif d <= 64 then
      put(" (quad double precision)");
    end if;
    new_line;
    put("  relaxation parameter k   : "); put(k,2); new_line;
    put("  accessibility constant a : "); put(a); new_line;
    put("  the target value for t   : "); put(t); new_line;
    put("  projective transformation : ");
         if proj then put("yes"); else put("no"); end if; new_line;
    put("Type d, k, a, t, or p to change, preceded by i for info." 
        & "  Type 0 to exit : ");
    Ask_Alternative(res,"dkatp0",'i');
    return res;
  end Show_Menu_and_Prompt_Choice;

  function Show_Menu_and_Prompt_Choice 
              ( d,k : natural32; a,t : Complex_Number ) return string is

  -- DESCRIPTION :
  --   Shows the menu with the current settings of the parameters
  --   and prompts the user to make choice, returned in a string.

    res : string(1..2) := "  ";

  begin
    new_line;
    put_line("MENU with settings for homotopy :");
    put("  number of decimal places : "); put(d,2);
    if d <= 16 then
      put(" (standard double precision)");
    elsif d <= 32 then
      put(" (double double precision)");
    elsif d <= 64 then
      put(" (quad double precision)");
    end if;
    new_line;
    put("  relaxation parameter k   : "); put(k,2); new_line;
    put("  accessibility constant a : "); put(a); new_line;
    put("  the target value for t   : "); put(t); new_line;
    put("Type d, k, a, or t to change, preceded by i for info." 
        & "  Type 0 to exit : ");
    Ask_Alternative(res,"dkat0",'i');
    return res;
  end Show_Menu_and_Prompt_Choice;

  function Show_Menu_and_Prompt_Choice 
              ( qd : boolean; k : natural32;
                a,t : Complex_Number ) return string is

  -- DESCRIPTION :
  --   Shows the menu with the current settings of the parameters
  --   and prompts the user to make choice, returned in a string.

    res : string(1..2) := "  ";

  begin
    new_line;
    put_line("MENU with settings for homotopy :");
    put("  number of decimal places : ");
    if qd 
     then put_line("64 (quad double precision)");
     else put_line("32 (double double precision)");
    end if;
    put("  relaxation parameter k   : "); put(k,2); new_line;
    put("  accessibility constant a : "); put(a); new_line;
    put("  the target value for t   : "); put(t); new_line;
    put("Type k, a, or t to change, preceded by i for info." 
        & "  Type 0 to exit : ");
    Ask_Alternative(res,"kat0",'i');
    return res;
  end Show_Menu_and_Prompt_Choice;

  procedure Write_Homotopy_Values
                ( file : in file_type; d,k : in natural32;
                  a,t : in Complex_Number; proj : in boolean ) is

  -- DESCRIPTION :
  --   Writes the values of the homotopy parameters to file.

  begin
    new_line(file);
    put_line(file,"HOMOTOPY PARAMETERS :");
    put(file,"  d : "); put(file,d,2);
    if d <= 16 then
      put(file," (standard double precision)");
    elsif d <= 32 then
      put(file," (double double precision)");
    elsif d <= 64 then
      put(file," (quad double precision)");
    end if;
    new_line(file);
    put(file,"  k : "); put(file,k,2); new_line(file);
    put(file,"  a : "); put(file,a);   new_line(file);
    put(file,"  t : "); put(file,t);   new_line(file);
    if proj
     then put_line(file,"  projective transformation");
     else put_line(file,"  no projective transformation");
    end if;
  end Write_Homotopy_Values;

  procedure Write_Homotopy_Values
                ( file : in file_type; d,k : in natural32;
                  a,t : in Complex_Number ) is

  -- DESCRIPTION :
  --   Writes the values of the homotopy parameters to file.

  begin
    new_line(file);
    put_line(file,"HOMOTOPY PARAMETERS :");
    put(file,"  d : "); put(file,d,2);
    if d <= 16 then
      put(file," (standard double precision)");
    elsif d <= 32 then
      put(file," (double double precision)");
    elsif d <= 64 then
      put(file," (quad double precision)");
    end if;
    new_line(file);
    put(file,"  k : "); put(file,k,2); new_line(file);
    put(file,"  a : "); put(file,a);   new_line(file);
    put(file,"  t : "); put(file,t);   new_line(file);
  end Write_Homotopy_Values;

  procedure Menu_for_Homotopy_Settings
                ( file : in file_type; d,k : in out natural32;
                  a,t : in out Complex_Number; proj : in out boolean ) is

    ans : character;
    choice : string(1..2) := "  ";

  begin
    Show_Banner;
    loop
      choice := Show_Menu_and_Prompt_Choice(d,k,a,t,proj);
      exit when choice(1) = '0';
      case choice(1) is
        when 'd' => 
          put("Give the number of decimal places (<= 16 is standard) : ");
          Read_Positive(integer(d));
        when 'k' => 
          put("Give a value (1,2, or 3) for the relaxation constant k : ");
          Read_Positive(integer(k));
        when 'a' =>
          put_line("Reading a complex value for the regularity constant a.");
          loop
            Read_Complex(a);
            exit when AbsVal(a) /= 0.0;
            put_line("The value 0 for a is not allowed! Try again.");
          end loop;
        when 't' =>
          put_line("Reading the (complex) target value for t.");
          Read_Complex(t);
        when 'p' =>
          put("Do you want a projective transformation? ");
          Ask_Yes_or_No(ans); proj := (ans ='y');
        when 'i' =>
          new_line;
          case choice(2) is
            when 'd' => Info_on_Precision;
            when 'k' => Info_on_Power;
            when 'a' => Info_on_Constant;
            when 't' => Info_on_Target;
            when 'p' => Info_on_Projective_Transformation;
            when '0' => Info_on_Exit;
            when others => null;
          end case;
        when others => null;
      end case;
    end loop;
    Write_Homotopy_Values(file,d,k,a,t,proj);
  end Menu_for_Homotopy_Settings;

  procedure Menu_for_Homotopy_Settings
                ( file : in file_type; d,k : in out natural32;
                  a,t : in out Complex_Number ) is

    choice : string(1..2) := "  ";

  begin
    Show_Banner;
    loop
      choice := Show_Menu_and_Prompt_Choice(d,k,a,t);
      exit when choice(1) = '0';
      case choice(1) is
        when 'd' => 
          put("Give the number of decimal places (<= 16 is standard) : ");
          Read_Positive(integer(d));
        when 'k' => 
          put("Give a value (1,2, or 3) for the relaxation constant k : ");
          Read_Positive(integer(k));
        when 'a' =>
          put_line("Reading a complex value for the regularity constant a.");
          loop
            Read_Complex(a);
            exit when AbsVal(a) /= 0.0;
            put_line("The value 0 for a is not allowed! Try again.");
          end loop;
        when 't' =>
          put_line("Reading the (complex) target value for t.");
          Read_Complex(t);
        when 'i' =>
          new_line;
          case choice(2) is
            when 'd' => Info_on_Precision;
            when 'k' => Info_on_Power;
            when 'a' => Info_on_Constant;
            when 't' => Info_on_Target;
            when '0' => Info_on_Exit;
            when others => null;
          end case;
        when others => null;
      end case;
    end loop;
    Write_Homotopy_Values(file,d,k,a,t);
  end Menu_for_Homotopy_Settings;

  procedure Menu_for_Homotopy_Settings
                ( file : in file_type;
                  qd : in boolean;  k : in out natural32;
                  a,t : in out Complex_Number ) is

    choice : string(1..2) := "  ";

  begin
    Show_Banner;
    loop
      choice := Show_Menu_and_Prompt_Choice(qd,k,a,t);
      exit when choice(1) = '0';
      case choice(1) is
        when 'k' => 
          put("Give a value (1,2, or 3) for the relaxation constant k : ");
          Read_Positive(integer(k));
        when 'a' =>
          put_line("Reading a complex value for the regularity constant a.");
          loop
            Read_Complex(a);
            exit when AbsVal(a) /= 0.0;
            put_line("The value 0 for a is not allowed! Try again.");
          end loop;
        when 't' =>
          put_line("Reading the (complex) target value for t.");
          Read_Complex(t);
        when 'i' =>
          new_line;
          case choice(2) is
            when 'k' => Info_on_Power;
            when 'a' => Info_on_Constant;
            when 't' => Info_on_Target;
            when '0' => Info_on_Exit;
            when others => null;
          end case;
        when others => null;
      end case;
    end loop;
    if qd
     then Write_Homotopy_Values(file,64,k,a,t);
     else Write_Homotopy_Values(file,32,k,a,t);
    end if;
  end Menu_for_Homotopy_Settings;

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 p,q : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 target : out Complex_Number; verbose : in integer32 := 0 ) is

    k,d : natural32;
    a,t : Complex_Number;

  begin
    if verbose > 0 then
      put("-> in drivers_for_homotopy_creation.");
      put_line("Driver_for_Homotopy_Construction 1 ...");
    end if;
    Default_Homotopy_Settings(d,k,a,t);
    Menu_for_Homotopy_Settings(file,false,k,a,t);
    target := t;
    declare
      dd_a : constant DoblDobl_Complex_Numbers.Complex_Number 
           := Standard_to_DoblDobl_Complex(a);
    begin
      DoblDobl_Homotopy.Create(p,q,k,dd_a);
    end;
  end Driver_for_Homotopy_Construction;

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 p,q : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 target : out Complex_Number; verbose : in integer32 := 0 ) is

    k,d : natural32;
    a,t : Complex_Number;

  begin
    if verbose > 0 then
      put("-> in drivers_for_homotopy_creation.");
      put_line("Driver_for_Homotopy_Construction 2 ...");
    end if;
    Default_Homotopy_Settings(d,k,a,t);
    Menu_for_Homotopy_Settings(file,true,k,a,t);
    target := t;
    declare
      qd_a : constant QuadDobl_Complex_Numbers.Complex_Number 
           := Standard_to_QuadDobl_Complex(a);
    begin
      QuadDobl_Homotopy.Create(p,q,k,qd_a);
    end;
  end Driver_for_Homotopy_Construction;

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type; dp : in natural32;
                 p,q : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 target : out Complex_Number; verbose : in integer32 := 0 ) is

    k,d : natural32;
    a,t : Complex_Number;

  begin
    if verbose > 0 then
      put("-> in drivers_for_homotopy_creation.");
      put_line("Driver_for_Homotopy_Construction 3 ...");
    end if;
    Default_Homotopy_Settings(d,k,a,t); d := dp;
    Menu_for_Homotopy_Settings(file,d,k,a,t);
    target := t;
    declare
      mp_a : constant Multprec_Complex_Numbers.Complex_Number 
           := Multprec_Complex_Number_Tools.Create(a);
    begin
      Multprec_Homotopy.Create(p,q,k,mp_a);
    end;
  end Driver_for_Homotopy_Construction;

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type; p,q : in out Laur_Sys;
                -- qsols : in out Solution_List;
                 target : out Complex_Number; deci : in out natural32;
                 verbose : in integer32 := 0 ) is

    k,d : natural32;
    a,t : Complex_Number;
    prt : boolean;
   -- size : natural;

  begin
    if verbose > 0 then
      put("-> in drivers_for_homotopy_creation.");
      put_line("Driver_for_Homotopy_Construction 4 ...");
    end if;
    Default_Homotopy_Settings(d,k,a,t,prt);
    if deci /= 0   
     then d := deci; -- respect the preset value of the precision
    end if;
    Menu_for_Homotopy_Settings(file,d,k,a,t,prt);
    target := t;
    if d <= 16 then
--      if prt then
--        Projective_Transformation(p);
--        Projective_Transformation(q);
--        Projective_Transformation(qsols);
--        declare
--          pp,sysp : Poly_Sys(p'first..p'last+1);
--        begin
--          pp := Add_Random_Hyperplanes(p,1,true);
--          sysp := Add_Standard_Hyperplanes(q,1);
--          Standard_Laurent_Homotopy.Create(pp,sysp,k,a);
--          Clear(pp); Clear(sysp);
--        end;
--      else
        Standard_Laurent_Homotopy.Create(p,q,k,a);
--      end if;
--    else
--      size := Multprec_Floating_Numbers.Decimal_to_Size(d);
--      declare
--        mp,mq : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);
--        aa : constant Multprec_Complex_Numbers.Complex_Number := Create(a);
--      begin
--        mp := Convert(p);  mq := Convert(q);
--        Set_Size(mp,size); Set_Size(mq,size);
--        Multprec_Homotopy.Create(mp,mq,k,aa);
--      end;
    end if;
    deci := d;
  end Driver_for_Homotopy_Construction;

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 ls : in String_Splitters.Link_to_Array_of_Strings;
                 p,q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols : in out Solution_List; target : out Complex_Number;
                 deci : in out natural32; verbose : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;
    a,t : Complex_Number;
    prt : boolean;
    k,d,size : natural32;

  begin
    if verbose > 0 then
      put("-> in drivers_for_homotopy_creation.");
      put_line("Driver_for_Homotopy_Construction 5 ...");
    end if;
    Default_Homotopy_Settings(d,k,a,t,prt);
    if deci > 0
     then d := deci; -- respect the preset value of the precision
    end if;
    Menu_for_Homotopy_Settings(file,d,k,a,t,prt);
    target := t;
    if d <= 16 then
      if prt then
        Projective_Transformation(p);
        Projective_Transformation(q);
        Projective_Transformation(qsols);
        declare
          pp,sysp : Poly_Sys(p'first..p'last+1);
        begin
          pp := Add_Random_Hyperplanes(p,1,true);
          sysp := Add_Standard_Hyperplanes(q,1);
          Standard_Homotopy.Create(pp,sysp,k,a);
          Standard_Homotopy.Create(sysp,pp,k,a);
          Clear(pp); Clear(sysp);
        end;
      else
        Standard_Homotopy.Create(p,q,k,a);
       -- Standard_Coefficient_Homotopy.Create(q,p,k,a);
        Standard_Coefficient_Homotopy.Create(p,q,k,a);
      end if;
    elsif d <= 32 then
      declare
        nv_p : constant natural32 := Number_of_Unknowns(p(p'first));
        dd_p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range)
             := DoblDobl_Complex_Poly_Strings.Parse(nv_p,ls.all);
        dd_q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(q'range)
             := Standard_Poly_Sys_to_DoblDobl_Complex(q);
         -- assuming the start system has only random coefficients ...
        dd_a : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(a);
      begin
        DoblDobl_Homotopy.Create(dd_p,dd_q,k,dd_a);
      end;
    elsif d <= 64 then
      declare
        nv_p : constant natural32 := Number_of_Unknowns(p(p'first));
        qd_p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range)
             := QuadDobl_Complex_Poly_Strings.Parse(nv_p,ls.all);
        qd_q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(q'range)
             := Standard_Poly_Sys_to_QuadDobl_Complex(q);
          -- assuming the start system has only random coefficients ...
        qd_a : constant QuadDobl_Complex_Numbers.Complex_Number
             := Standard_to_QuadDobl_Complex(a);
      begin
       -- put("creating quaddobl homotopy with k = "); put(k,1); new_line;
        QuadDobl_Homotopy.Create(qd_p,qd_q,k,qd_a);
      end;
    else
      size := Multprec_Floating_Numbers.Decimal_to_Size(d);
      declare
        nv_p : constant natural32 := Number_of_Unknowns(p(p'first));
        mp : constant Multprec_Complex_Poly_Systems.Poly_Sys(p'range)
           := Multprec_Complex_Poly_Strings.Parse(nv_p,size,ls.all);
        mq : Multprec_Complex_Poly_Systems.Poly_Sys(p'range)
           := Convert(q);
         -- assuming the start system has only random coefficients ...
        aa : constant Multprec_Complex_Numbers.Complex_Number := Create(a);
      begin
        Set_Size(mq,size);
        Multprec_Homotopy.Create(mp,mq,k,aa);
      end;
    end if;
    deci := d;
  end Driver_for_Homotopy_Construction;

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 ls : in String_Splitters.Link_to_Array_of_Strings;
                 p,q : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out Solution_List; target : out Complex_Number;
                 deci : in out natural32; verbose : in integer32 := 0 ) is

    a,t : Complex_Number;
    prt : boolean;
    k,d : natural32;

  begin
    if verbose > 0 then
      put("-> in drivers_for_homotopy_creation.");
      put_line("Driver_for_Homotopy_Construction 6 ...");
    end if;
    Default_Homotopy_Settings(d,k,a,t,prt);
    if deci > 0
     then d := deci; -- respect the preset value of the precision
    end if;
    Menu_for_Homotopy_Settings(file,d,k,a,t,prt);
    target := t;
    if d <= 16 then
      Standard_Laurent_Homotopy.Create(p,q,k,a);
    elsif d <= 32 then
      declare
        nv_p : constant natural32 := Number_of_Unknowns(p(p'first));
        dd_p : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range)
             := DoblDobl_Complex_Laur_Strings.Parse(nv_p,ls.all);
        dd_q : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(q'range)
             := Standard_Laur_Sys_to_DoblDobl_Complex(q);
         -- assuming the start system has only random coefficients ...
        dd_a : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(a);
      begin
        DoblDobl_Laurent_Homotopy.Create(dd_p,dd_q,k,dd_a);
      end;
    elsif d <= 64 then
      declare
        nv_p : constant natural32 := Number_of_Unknowns(p(p'first));
        qd_p : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range)
             := QuadDobl_Complex_Laur_Strings.Parse(nv_p,ls.all);
        qd_q : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(q'range)
             := Standard_Laur_Sys_to_QuadDobl_Complex(q);
          -- assuming the start system has only random coefficients ...
        qd_a : constant QuadDobl_Complex_Numbers.Complex_Number
             := Standard_to_QuadDobl_Complex(a);
      begin
       -- put("creating quaddobl homotopy with k = "); put(k,1); new_line;
        QuadDobl_Laurent_Homotopy.Create(qd_p,qd_q,k,qd_a);
      end;
    end if;
    deci := d;
  end Driver_for_Homotopy_Construction;

end Drivers_for_Homotopy_Creation;
