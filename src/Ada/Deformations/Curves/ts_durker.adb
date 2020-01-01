with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Durand_Kerner;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;
with DoblDobl_Durand_Kerner;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with QuadDobl_Durand_Kerner;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Durand_Kerner;
with Hybrid_Durand_Kerner;

procedure ts_durker is

-- DESCRIPTION :
--   Test on the solver for polynomial equations in one variable.

  procedure Standard_Read ( cv : in out Standard_Complex_Vectors.Vector ) is
  begin
    for i in cv'range loop
      put(' '); put(i,1); put(" : ");
      Standard_Complex_Numbers_io.get(cv(i));
    end loop;
  end Standard_Read;

  procedure Multprec_Read ( cv : in out Multprec_Complex_Vectors.Vector ) is
  begin
    for i in cv'range loop
      put(' '); put(i,1); put(" : ");
      Multprec_Complex_Numbers_io.get(cv(i));
    end loop;
  end Multprec_Read;

  procedure Standard_Write
               ( step : in natural32;
                 z,res : in Standard_Complex_Vectors.Vector ) is

    absres : double_float;

  begin
    put("Output after step  "); put(step,1); put_line(" :");
    put_line
   ("------------------------------------------------------------------------");
    put_line
   ("|    APPROXIMATED ROOTS                        |     RESIDUALS         |");
    put_line
   ("------------------------------------------------------------------------");
    for i in z'range loop
      absres := AbsVal(res(i));
      put("| "); put(z(i)); put(" | "); put(absres); put(" |");
      new_line;
    end loop;
    put_line
   ("------------------------------------------------------------------------");
  end Standard_Write;

  procedure DoblDobl_Write
               ( step : in natural32;
                 z,res : in DoblDobl_Complex_Vectors.Vector ) is

    absres : double_float;
    dd_absres : double_double;

  begin
    put("Output after step  "); put(step,1); put_line(" :");
    put("------------------------------------");
    put_line
   ("------------------------------------------------------------------------");
    put("|    APPROXIMATED ROOTS                        ");
    put("                                    ");
    put_line("|     RESIDUALS         |");
    put("------------------------------------");
    put_line
   ("------------------------------------------------------------------------");
    for i in z'range loop
      dd_absres := DoblDobl_Complex_Numbers.AbsVal(res(i));
      absres := to_double(dd_absres);
      put("| ");
      if DoblDobl_Complex_Numbers.REAL_PART(z(i)) >= 0.0
       then put(" ");
      end if;
      if DoblDobl_Complex_Numbers.IMAG_PART(z(i)) >= 0.0
       then put(" ");
      end if;
      put(z(i)); put(" | "); put(absres); put(" |");
      new_line;
    end loop;
    put("------------------------------------");
    put_line
   ("------------------------------------------------------------------------");
  end DoblDobl_Write;

  procedure QuadDobl_Write
               ( step : in natural32;
                 z,res : in QuadDobl_Complex_Vectors.Vector ) is

    absres : double_float;
    qd_absres : quad_double;

  begin
    put("Output after step  "); put(step,1); put_line(" :");
    put("--------------------------------");
    put("--------------------------------------------------------------------");
    put_line
   ("------------------------------------------------------------------------");
    put("|    APPROXIMATED ROOTS                        ");
    put("                                ");
    put("                                                                    ");
    put_line("|     RESIDUALS         |");
    put("--------------------------------");
    put("--------------------------------------------------------------------");
    put_line
   ("------------------------------------------------------------------------");
    for i in z'range loop
      qd_absres := QuadDobl_Complex_Numbers.AbsVal(res(i));
      absres := to_double(qd_absres);
      put("| ");
      if QuadDobl_Complex_Numbers.REAL_PART(z(i)) >= 0.0
       then put(" ");
      end if;
      if QuadDobl_Complex_Numbers.IMAG_PART(z(i)) >= 0.0
       then put(" ");
      end if;
      put(z(i)); put(" | "); put(absres); put(" |");
      new_line;
    end loop;
    put("--------------------------------");
    put("--------------------------------------------------------------------");
    put_line
   ("------------------------------------------------------------------------");
  end QuadDobl_Write;

  procedure Multprec_Write
               ( step : in natural32;
                 z,res : in Multprec_Complex_Vectors.Vector ) is

    absres : Floating_Number;

  begin
    put("Output after step  "); put(step,1); put_line(" :");
    put_line
   ("------------------------------------------------------------------------");
    put_line
   ("|    APPROXIMATED ROOTS                                                |");
    put_line
   ("------------------------------------------------------------------------");
    for i in z'range loop
      put(z(i)); new_line;
    end loop;
    put_line
   ("------------------------------------------------------------------------");
    put_line
   ("|    RESIDUALS                                                         |");
    put_line
   ("------------------------------------------------------------------------");
    for i in z'range loop
      absres := AbsVal(res(i));
      put(absres); new_line;
      Clear(absres);
    end loop;
    put_line
   ("------------------------------------------------------------------------");
  end Multprec_Write;

  procedure stdk is 
    new Standard_Durand_Kerner.Reporting_Durand_Kerner(Standard_Write);
  procedure dddk is 
    new DoblDobl_Durand_Kerner.Reporting_Durand_Kerner(DoblDobl_Write);
  procedure qddk is 
    new QuadDobl_Durand_Kerner.Reporting_Durand_Kerner(QuadDobl_Write);
  procedure mpdk is 
    new Multprec_Durand_Kerner.Reporting_Durand_Kerner(Multprec_Write);
  procedure hbdk is
    new Hybrid_Durand_Kerner.Reporting_Durand_Kerner(Multprec_Write);

  procedure Standard_Root_Finder
              ( p : in Standard_Complex_Vectors.Vector;
                z,res : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Asks the user to give in the maximal number of steps, the accuracy,
  --   and whether intermediate output is wanted.

  -- ON ENTRY :
  --   p        vector or range 0..n with the coefficient of the polynomial;
  --   z        initial approximations of the roots;
  --   res      residuals at the initial approximations.

  -- ON RETURN :
  --   z        refined approximations of the roots of p;
  --   res      residuals at the approximations.

    max,nb : natural32 := 0;
    eps : double_float := 1.0;
    fail : boolean;
    ans : character;
    timer : Timing_Widget;

  begin
    put("Give the maximum number of steps : "); get(max);
    put("Give the required accuracy       : "); get(eps);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' then
      stdk(p,z,res,max,eps,nb,fail);
    else
      Standard_Durand_Kerner.Silent_Durand_Kerner(p,z,res,max,eps,nb,fail);
      Standard_Write(nb,z,res);
    end if;
    tstop(timer);
    if fail then
      put("Failed to reach accuracy ");
      put(eps,3); put(" within "); put(max,1); put_line(" steps.");
    else
      put("Reached accuracy ");
      put(eps,3); put(" after "); put(nb,1); put_line(" steps.");
    end if;
    print_times(Standard_Output,timer,"Durand-Kerner in standard arithmetic");
  end Standard_Root_Finder;

  procedure DoblDobl_Root_Finder
              ( p : in DoblDobl_Complex_Vectors.Vector;
                z,res : in out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Asks the user to give in the maximal number of steps, the accuracy,
  --   and whether intermediate output is wanted.

  -- ON ENTRY :
  --   n        degree of the polynomial;
  --   p        vector or range 0..n with the coefficient of the polynomial;
  --   z        initial approximations of the roots;
  --   res      residuals at the initial approximations.

  -- ON RETURN :
  --   z        refined approximations of the roots of p;
  --   res      residuals at the approximations.

    max,nb : natural32 := 0;
    eps : double_float := 1.0;
    fail : boolean;
    ans : character;
    timer : Timing_Widget;

  begin
    put("Give the maximum number of steps : "); get(max);
    put("Give the required accuracy       : "); get(eps);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' then
      dddk(p,z,res,max,eps,nb,fail);
    else
      DoblDobl_Durand_Kerner.Silent_Durand_Kerner(p,z,res,max,eps,nb,fail);
      DoblDobl_Write(nb,z,res);
    end if;
    tstop(timer);
    if fail then
      put("Failed to reach accuracy ");
      put(eps,3); put(" within "); put(max,1); put_line(" steps.");
    else
      put("Reached accuracy ");
      put(eps,3); put(" after "); put(nb,1); put_line(" steps.");
    end if;
    print_times(Standard_Output,timer,
                "Durand-Kerner in double double arithmetic");
  end DoblDobl_Root_Finder;

  procedure QuadDobl_Root_Finder
             ( p : in QuadDobl_Complex_Vectors.Vector;
               z,res : in out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Asks the user to give in the maximal number of steps, the accuracy,
  --   and whether intermediate output is wanted.

  -- ON ENTRY :
  --   n       degree of the polynomial;
  --   p       vector or range 0..n with the coefficient of the polynomial;
  --   z       initial approximations of the roots;
  --   res     residuals at the initial approximations.

  -- ON RETURN :
  --   z       refined approximations of the roots of p;
  --   res     residuals at the approximations.

    max,nb : natural32 := 0;
    eps : double_float := 1.0;
    fail : boolean;
    ans : character;
    timer : Timing_Widget;

  begin
    put("Give the maximum number of steps : "); get(max);
    put("Give the required accuracy       : "); get(eps);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' then
      qddk(p,z,res,max,eps,nb,fail);
    else
      QuadDobl_Durand_Kerner.Silent_Durand_Kerner(p,z,res,max,eps,nb,fail);
      QuadDobl_Write(nb,z,res);
    end if;
    tstop(timer);
    if fail then
      put("Failed to reach accuracy ");
      put(eps,3); put(" within "); put(max,1); put_line(" steps.");
    else
      put("Reached accuracy ");
      put(eps,3); put(" after "); put(nb,1); put_line(" steps.");
    end if;
    print_times(Standard_Output,timer,
                "Durand-Kerner in quad double arithmetic");
  end QuadDobl_Root_Finder;

  procedure Multprec_Root_Finder
              ( p : Multprec_Complex_Vectors.Vector;
                z,res : in out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Specifications are the same as the Standard_Root_Finder.

    max,nb : natural32 := 0;
    eps : Floating_Number;
    fail : boolean;
    ans : character;
    timer : Timing_Widget;

  begin
    put("Give the maximum number of steps : "); get(max);
    put("Give the required accuracy       : "); get(eps);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' then
      mpdk(p,z,res,max,eps,nb,fail);
    else
      Multprec_Durand_Kerner.Silent_Durand_Kerner(p,z,res,max,eps,nb,fail);
      Multprec_Write(nb,z,res);
    end if;
    tstop(timer);
    if fail then
      put("Failed to reach accuracy ");
      put(eps,3); put(" within "); put(max,1); put_line(" steps.");
    else
      put("Reached accuracy ");
      put(eps,3); put(" after "); put(nb,1); put_line(" steps.");
    end if;
    print_times(Standard_Output,timer,"Multi-precision Durand-Kerner");
  end Multprec_Root_Finder;

  procedure Hybrid_Root_Finder
              ( p : Multprec_Complex_Vectors.Vector;
                z,res : in out Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Specifications are the same as the Multprec_Root_Finder.

    max,nb : natural32 := 0;
    eps : Floating_Number;
    fail : boolean;
    ans : character;
    timer : Timing_Widget;

  begin
    put("Give the maximum number of steps : "); get(max);
    put("Give the required accuracy       : "); get(eps);
    put("Do you want intermediate output ? (y/n) "); Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' then
      hbdk(p,z,res,max,eps,nb,fail);
    else
      Hybrid_Durand_Kerner.Silent_Durand_Kerner(p,z,res,max,eps,nb,fail);
      Multprec_Write(nb,z,res);
    end if;
    tstop(timer);
    if fail then
      put("Failed to reach accuracy ");
      put(eps,3); put(" within "); put(max,1); put_line(" steps.");
    else
      put("Reached accuracy ");
      put(eps,3); put(" after "); put(nb,1); put_line(" steps.");
    end if;
    print_times(Standard_Output,timer,"Hybrid Durand-Kerner");
  end Hybrid_Root_Finder;

  procedure Standard_User_Given_Test ( n : in integer32 ) is

    p : Standard_Complex_Vectors.Vector(0..n);
    z : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    res : Standard_Complex_Vectors.Vector(1..n) := z;

  begin
    put_line("p(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n");
    put("Give "); put(n+1,1);
    put_line(" complex coefficients of p(x) :");
    Standard_Read(p);
    Standard_Root_Finder(p,z,res);
  end Standard_User_Given_Test;

  procedure Multprec_User_Given_Test
              ( n : in integer32; size : in natural32 ) is

    p : Multprec_Complex_Vectors.Vector(0..n);
    z : Multprec_Complex_Vectors.Vector(1..n) := Random_Vector(1,n,size);
    res : Multprec_Complex_Vectors.Vector(1..n) := z;

  begin
    put_line("p(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n");
    put("Give "); put(n+1,1);
    put_line(" complex coefficients of p(x) :");
    Multprec_Read(p);
    Set_Size(p,size);
    Multprec_Root_Finder(p,z,res);
  end Multprec_User_Given_Test;

  procedure Standard_Random_Test ( n : in integer32 ) is

    p : constant Standard_Complex_Vectors.Vector(0..n) := Random_Vector(0,n);
    z : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    res : Standard_Complex_Vectors.Vector(1..n) := z;

  begin
    Standard_Root_Finder(p,z,res);
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test ( n : in integer32 ) is

    p : constant DoblDobl_Complex_Vectors.Vector(0..n) := Random_Vector(0,n);
    z : DoblDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    res : DoblDobl_Complex_Vectors.Vector(1..n) := z;

  begin
    DoblDobl_Root_Finder(p,z,res);
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test ( n : in integer32 ) is

    p : constant QuadDobl_Complex_Vectors.Vector(0..n) := Random_Vector(0,n);
    z : QuadDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    res : QuadDobl_Complex_Vectors.Vector(1..n) := z;

  begin
    QuadDobl_Root_Finder(p,z,res);
  end QuadDobl_Random_Test;

  procedure Multprec_Random_Test
              ( n : in integer32; size : natural32 ) is

    p : constant Multprec_Complex_Vectors.Vector(0..n)
      := Random_Vector(0,n,size);
    z : Multprec_Complex_Vectors.Vector(1..n) := Random_Vector(1,n,size);
    res : Multprec_Complex_Vectors.Vector(1..n) := z;

  begin
    Multprec_Root_Finder(p,z,res);
  end Multprec_Random_Test;

  procedure Hybrid_Random_Test
              ( n : in integer32; size : in natural32 ) is

    p : constant Multprec_Complex_Vectors.Vector(0..n)
      := Random_Vector(0,n,size);
    z : Multprec_Complex_Vectors.Vector(1..n) := Random_Vector(1,n,size);
    res : Multprec_Complex_Vectors.Vector(1..n) := z;

  begin
    Hybrid_Root_Finder(p,z,res);
  end Hybrid_Random_Test;

  procedure Two_Stage_Hybrid_Random_Test
              ( n : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Applies first the standard Durand-Kerner and then refines
  --   with the multi-precision Durand-Kerner.

    mpp : constant Multprec_Complex_Vectors.Vector(0..n)
        := Random_Vector(0,n,size);
    stp : constant Standard_Complex_Vectors.Vector(0..n) := Round(mpp);
    stz : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    stres : Standard_Complex_Vectors.Vector(1..n) := stz;
    mpz,mpres : Multprec_Complex_Vectors.Vector(1..n);

  begin
    new_line;
    put_line("Starting initial approximations in standard arithmetic.");
    Standard_Root_Finder(stp,stz,stres);
    new_line;
    put_line("Refining initial approximations in multi-precision arithmetic.");
    mpz := Create(stz);
    mpres := Create(stres);
    Multprec_Root_Finder(mpp,mpz,mpres);
  end Two_Stage_Hybrid_Random_Test;

  procedure Main is

    n : integer32 := 0;
    size,deci : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("The method of Durand-Kerner to compute all complex roots.");
    loop
      new_line;
      put_line("Choose one of the following options : ");
      put_line("  0. exit this program.");
      put_line("  1. given polynomial with standard numbers.");
      put_line("  2. given polynomial with multi-precision numbers.");
      put_line("  3. random polynomial with standard numbers.");
      put_line("  4. random polynomial with double double numbers.");
      put_line("  5. random polynomial with quad double numbers.");
      put_line("  6. random polynomial with multi-precision numbers.");
      put_line("  7. hybrid standard/multi-precision on random polynomial.");
      put("Type 0, 1, 2, 3, 4, 5, 6, or 7 to select : ");
      Ask_Alternative(ans,"01234567");
      exit when ans = '0';
      new_line;
      put("Give the degree of p(x) : "); get(n);
      if (ans = '2') or (ans = '6') or (ans = '7') then
        put("Give the number of decimal places : "); get(deci);
        size := Decimal_to_Size(deci);
      end if;
      case ans is
        when '1' => Standard_User_Given_Test(n);
        when '2' => Multprec_User_Given_Test(n,size);
        when '3' => Standard_Random_Test(n);
        when '4' => DoblDobl_Random_Test(n);
        when '5' => QuadDobl_Random_Test(n);
        when '6' => Multprec_Random_Test(n,size);
        when '7' => Hybrid_Random_Test(n,size);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_durker;
