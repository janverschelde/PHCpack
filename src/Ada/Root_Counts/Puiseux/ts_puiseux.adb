with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Vectors;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Norms;
with Symbol_Table;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Poly_Systems;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Poly_Laur_Convertors;      use DoblDobl_Poly_Laur_Convertors;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Poly_Laur_Convertors;      use QuadDobl_Poly_Laur_Convertors;
with Random_Polynomial_Systems;          use Random_Polynomial_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Standard_Dense_Series_Vectors;
with Standard_Series_Vector_Functions;
with Standard_Dense_Series_VecVecs;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Series_Vector_Functions;
with DoblDobl_Dense_Series_VecVecs;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Series_Vector_Functions;
with QuadDobl_Dense_Series_VecVecs;
with Regular_Solution_Curves_Series;     use Regular_Solution_Curves_Series;

procedure ts_puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm,
--   for regular solution curves, defined by complete intersections,
--   in Noether position, with sufficiently general coefficients.

  procedure Tropisms_by_Mixed_Cells
              ( file : in file_type; 
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 ) is

  -- DESCRIPTION :
  --   Computes a regular mixed cell configuration for
  --   the supports in sup, with some writing to file.
  --   The mixed volume is in mv.

  begin
    Mixed_Cell_Tropisms(file,sup,mcc,mv);
    new_line(file);
    put(file,"The number of pretropisms : ");
    put(file,Length_Of(mcc),1); new_line(file);
    put(file,"The number of series : ");
    put(file,mv,1); new_line(file);
  end Tropisms_by_Mixed_Cells;

  procedure Tropisms_by_Mixed_Cells
              ( sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Computes a regular mixed cell configuration for
  --   the supports in sup, with some writing to screen if report.
  --   The mixed volume is in mv.

  begin
    Mixed_Cell_Tropisms(report,sup,mcc,mv);
    if report then
      new_line;
      put("The number of pretropisms : ");
      put(Length_Of(mcc),1); new_line;
      put("The number of series : "); put(mv,1); new_line;
    end if;
  end Tropisms_by_Mixed_Cells;

  function Standard_Residual
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                s : Standard_Dense_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns the residual of the evaluation of s at t in p,
  --   or ||p(s(t),t)||, with weights for the exponents in w.

    res : double_float := 0.0;
    x : constant Standard_Complex_Vectors.Vector(s'range)
      := Standard_Series_Vector_Functions.Eval(s,w,t);
    xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);
    y : Standard_Complex_Vectors.Vector(p'range);

  begin
    xt(x'range) := x;
    xt(xt'last) := Standard_Complex_Numbers.Create(t);
    y := Standard_Complex_Laur_SysFun.Eval(p,xt);
    res := Standard_Complex_Vector_Norms.Norm2(y);
    return res;
  end Standard_Residual;

  function DoblDobl_Residual
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : DoblDobl_Dense_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : double_double ) return double_double is

  -- DESCRIPTION :
  --   Returns the residual of the evaluation of s at t in p,
  --   or ||p(s(t),t)||.

    res : double_double := create(0.0);
    x : constant DoblDobl_Complex_Vectors.Vector(s'range)
      := DoblDobl_Series_Vector_Functions.Eval(s,w,t);
    xt : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);
    y : DoblDobl_Complex_Vectors.Vector(p'range);

  begin
    xt(x'range) := x;
    xt(xt'last) := DoblDobl_Complex_Numbers.Create(t);
    y := DoblDobl_Complex_Laur_SysFun.Eval(p,xt);
    res := DoblDobl_Complex_Vector_Norms.Norm2(y);
    return res;
  end DoblDobl_Residual;

  function QuadDobl_Residual
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : QuadDobl_Dense_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : quad_double ) return quad_double is

  -- DESCRIPTION :
  --   Returns the residual of the evaluation of s at t in p,
  --   or ||p(s(t),t)||.

    res : quad_double := create(0.0);
    x : constant QuadDobl_Complex_Vectors.Vector(s'range)
      := QuadDobl_Series_Vector_Functions.Eval(s,w,t);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);
    y : QuadDobl_Complex_Vectors.Vector(p'range);

  begin
    xt(x'range) := x;
    xt(xt'last) := QuadDobl_Complex_Numbers.Create(t);
    y := QuadDobl_Complex_Laur_SysFun.Eval(p,xt);
    res := QuadDobl_Complex_Vector_Norms.Norm2(y);
    return res;
  end QuadDobl_Residual;

  function Standard_Residuals
              ( file : file_type;
                p : Standard_Complex_Laur_Systems.Laur_Sys;
                s : Standard_Dense_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_float ) return double_float is

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, with weights for the exponents in w,
  --   computed in standard double precision.
  --   Output is written to file.
  --   Returns the sum of the residuals.

    sum : double_float := 0.0;
    res : double_float;

  begin
    for k in s'range loop
      res := Standard_Residual(p,s(k).all,w(k).all,t);
      put(file,"Residual at series ");
      put(file,k,1); put(file," : ");
      put(file,res,3); new_line(file);
      sum := sum + res;
    end loop;
    put(file,"Sum of all residuals : ");
    put(file,sum,3); new_line(file);
    return sum;
  end Standard_Residuals;

  function Standard_Residuals
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                s : Standard_Dense_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_float; report : in boolean )
              return double_float is

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, with weights for the exponents in w,
  --   computed in standard double precision.
  --   Output is written to screen if report.
  --   Returns the sum of the residuals.

    sum : double_float := 0.0;
    res : double_float;

  begin
    for k in s'range loop
      res := Standard_Residual(p,s(k).all,w(k).all,t);
      if report then
        put("Residual at series "); put(k,1); put(" : ");
        put(res,3); new_line;
      end if;
      sum := sum + res;
    end loop;
    if report then
      put("Sum of all residuals : ");
      put(sum,3); new_line;
    end if;
    return sum;
  end Standard_Residuals;

  function DoblDobl_Residuals
              ( file : file_type;
                p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : DoblDobl_Dense_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_double ) return double_double is

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, with weights for the exponents in w,
  --   computed in double double precision.
  --   Output is written to file.
  --   Returns the sum of the residuals.

    sum : double_double := create(0.0);
    res : double_double;

  begin
    for k in s'range loop
      res := DoblDobl_Residual(p,s(k).all,w(k).all,t);
      put(file,"Residual at series ");
      put(file,k,1); put(file," : ");
      put(file,res,3); new_line(file);
      sum := sum + res;
    end loop;
    put(file,"Sum of all residuals : ");
    put(file,sum,3); new_line(file);
    return sum;
  end DoblDobl_Residuals;

  function DoblDobl_Residuals
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : DoblDobl_Dense_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_double; report : in boolean )
              return double_double is

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, with weights for the exponents in w,
  --   computed in double double precision.
  --   Output is written to screen if report.
  --   Returns the sum of the residuals.

    sum : double_double := create(0.0);
    res : double_double;

  begin
    for k in s'range loop
      res := DoblDobl_Residual(p,s(k).all,w(k).all,t);
      if report then
        put("Residual at series "); put(k,1); put(" : ");
        put(res,3); new_line;
      end if;
      sum := sum + res;
    end loop;
    if report then
      put("Sum of all residuals : ");
      put(sum,3); new_line;
    end if;
    return sum;
  end DoblDobl_Residuals;

  function QuadDobl_Residuals
              ( file : file_type;
                p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : QuadDobl_Dense_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : quad_double ) return quad_double is

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, in quad double precision.
  --   Output is written to file.
  --   Returns the sum of the residuals.

    sum : quad_double := create(0.0);
    res : quad_double;

  begin
    for k in s'range loop
      res := QuadDobl_Residual(p,s(k).all,w(k).all,t);
      put(file,"Residual at series ");
      put(file,k,1); put(file," : ");
      put(file,res,3); new_line(file);
      sum := sum + res;
    end loop;
    put(file,"Sum of all residuals : ");
    put(file,sum,3); new_line(file);
    return sum;
  end QuadDobl_Residuals;

  function QuadDobl_Residuals
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : QuadDobl_Dense_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : quad_double; report : in boolean )
              return quad_double is

  -- DESCRIPTION :
  --   Computes the residuals of the series in s at t,
  --   evaluated in p, computed in quad double precision.
  --   Output is written to screen if report.
  --   Returns the sum of the residuals.

    sum : quad_double := create(0.0);
    res : quad_double;

  begin
    for k in s'range loop
      res := QuadDobl_Residual(p,s(k).all,w(k).all,t);
      if report then
        put("Residual at series "); put(k,1); put(" : ");
        put(res,3); new_line;
      end if;
      sum := sum + res;
    end loop;
    if report then
      put("Sum of all residuals : ");
      put(sum,3); new_line;
    end if;
    return sum;
  end QuadDobl_Residuals;

  procedure Standard_Test
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in standard double precision.  The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : Standard_Dense_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc);
      r : double_float;
    begin
      s := Series(file,p,mcc,mv,nit);
      put_line(file,"Evaluating series at 0.1 ...");
      r := Standard_Residuals(file,p,s,w,0.1);
    end;
  end Standard_Test;

  procedure DoblDobl_Test
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in double double precision.  The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : DoblDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc);
      t : constant double_double := create(0.1);
      r : double_double;
    begin
      s := Series(file,p,mcc,mv,nit);
      put_line(file,"Evaluating series at 0.1 ...");
      r := DoblDobl_Residuals(file,p,s,w,t);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in quad double precision.  The output is written to file.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : QuadDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc);
      t : constant quad_double := create(0.1);
      r : quad_double;
    begin
      s := Series(file,p,mcc,mv,nit);
      put_line(file,"Evaluating series at 0.1 ...");
      r := QuadDobl_Residuals(file,p,s,w,t);
    end;
  end QuadDobl_Test;

  procedure Standard_Test 
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in standard double precision.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : Standard_Dense_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc);
      r : double_float;
    begin
      s := Series(p,mcc,mv,nit,report);
      if report
       then put_line("Evaluating the series at 0.1 ...");
      end if;
      r := Standard_Residuals(p,s,w,0.1,report);
    end;
  end Standard_Test;

  procedure DoblDobl_Test 
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in double double precision.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : DoblDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc);
      t : constant double_double := create(0.1);
      r : double_double;
    begin
      s := Series(p,mcc,mv,nit,report);
      if report
       then put_line("Evaluating the series at 0.1 ...");
      end if;
      r := DoblDobl_Residuals(p,s,w,t,report);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test 
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the output information
  --   and test the computation of the pretropisms and the series,
  --   in quad double precision.
  --   The output is written to screen if report.

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : QuadDobl_Dense_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc);
      t : constant quad_double := create(0.1);
      r : quad_double;
    begin
      s := Series(p,mcc,mv,nit,report);
      if report
       then put_line("Evaluating the series at 0.1 ...");
      end if;
      r := QuadDobl_Residuals(p,s,w,t,report);
    end;
  end QuadDobl_Test;

  function Prompt_for_Precision return character is

  -- DESCRIPTION :
  --   Displays the menu for the available precision and
  --   returns '0', '1', or '2' for double, double double,
  --   or quad double precision.

    res : character;

  begin
    new_line;
    put_line("MENU to set the precision level :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(res,"012");
    return res;
  end Prompt_for_Precision;

  procedure Standard_Read
              ( lp : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system and returns the system,
  --   along with the number of polynomials in nq
  --   and the number of variables in nv.

    use Standard_Complex_Laurentials;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
  end Standard_Read;

  procedure DoblDobl_Read
              ( lp : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system and returns the system,
  --   along with the number of polynomials in nq
  --   and the number of variables in nv.

    use DoblDobl_Complex_Laurentials;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
  end DoblDobl_Read;

  procedure QuadDobl_Read
              ( lp : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent system and returns the system,
  --   along with the number of polynomials in nq
  --   and the number of variables in nv.

    use QuadDobl_Complex_Laurentials;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
  end QuadDobl_Read;

  procedure Standard_Test
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the kind of output (to file or to screen)
  --   and then computes the power series in standard double precision.
  --   The number of polynomials is in nq and
  --   the number of variables is in nv.

    ans : character;
    report : boolean;
    file : file_type;

  begin
    new_line;
    put("Do you want intermediate output to file ? (y/n) ");
    Ask_Yes_or_No(ans); report := (ans = 'y');
    if report then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
      put(file,natural32(nq),natural32(nv),p);
      Standard_Test(file,p);
    else
      new_line;
      put("Do you want intermediate output to screen ? (y/n) ");
      Ask_Yes_or_No(ans); report := (ans = 'y');
      Standard_Test(p,report);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the kind of output (to file or to screen)
  --   and then computes the power series in double double precision.
  --   The number of polynomials is in nq and
  --   the number of variables is in nv.

    ans : character;
    report : boolean;
    file : file_type;
  
  begin
    new_line;
    put("Do you want intermediate output to file ? (y/n) ");
    Ask_Yes_or_No(ans); report := (ans = 'y');
    if report then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
      put(file,natural32(nq),natural32(nv),p);
      DoblDobl_Test(file,p);
    else
      new_line;
      put("Do you want intermediate output to screen ? (y/n) ");
      Ask_Yes_or_No(ans); report := (ans = 'y');
      DoblDobl_Test(p,report);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the kind of output (to file or to screen)
  --   and then computes the power series in quad double precision.
  --   The number of polynomials is in nq and
  --   the number of variables is in nv.

    ans : character;
    report : boolean;
    file : file_type;
  
  begin
    new_line;
    put("Do you want intermediate output to file ? (y/n) ");
    Ask_Yes_or_No(ans); report := (ans = 'y');
    if report then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
      put(file,natural32(nq),natural32(nv),p);
      QuadDobl_Test(file,p);
    else
      new_line;
      put("Do you want intermediate output to screen ? (y/n) ");
      Ask_Yes_or_No(ans); report := (ans = 'y');
      QuadDobl_Test(p,report);
    end if;
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.
  --   Computations are done in standard double precision.

    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    Standard_Read(lp,nq,nv);
    if nv /= nq+1
     then put(nv,1); put(" /= "); put(nq,1); put(" + 1");
     else Standard_Test(lp.all,nq,nv);
    end if;
  end Standard_Main;
  
  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.
  --   Computations are done in double double precision.

    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    DoblDobl_Read(lp,nq,nv);
    if nv /= nq+1 
     then put(nv,1); put(" /= "); put(nq,1); put(" + 1");
     else DoblDobl_Test(lp.all,nq,nv);
    end if;
  end DoblDobl_Main;
  
  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.
  --   Computations are done in quad double precision.

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    QuadDobl_Read(lp,nq,nv);
    if nv /= nq+1 
     then put(nv,1); put(" /= "); put(nq,1); put(" + 1");
     else QuadDobl_Test(lp.all,nq,nv);
    end if;
  end QuadDobl_Main;

  procedure Standard_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the parameters to generate a system
  --   with random coefficients in standard double precision.
  --   The series are computed for the generated system.

    n,d,m,c : natural32 := 0;
    e : integer32 := 0;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

    use Standard_Complex_Laurentials;

  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    e := integer32(n) - 1;
    Standard_Generate_and_Show(n,d,m,c,e,lp);
    declare
      q : Standard_Complex_Laur_Systems.Laur_Sys(lp'range)
        := Polynomial_to_Laurent_System(lp.all);
      nv : constant integer32 := integer32(Number_of_Unknowns(q(q'first)));
    begin
      Standard_Test(q,q'last,nv);
    end;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the parameters to generate a system
  --   with random coefficients in double double precision.
  --   The series are computed for the generated system.

    n,d,m,c : natural32 := 0;
    e : integer32 := 0;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use DoblDobl_Complex_Laurentials;

  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    e := integer32(n) - 1;
    DoblDobl_Generate_and_Show(n,d,m,c,e,lp);
    declare
      q : DoblDobl_Complex_Laur_Systems.Laur_Sys(lp'range)
        := Polynomial_to_Laurent_System(lp.all);
      nv : constant integer32 := integer32(Number_of_Unknowns(q(q'first)));
    begin
      DoblDobl_Test(q,q'last,nv);
    end;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the parameters to generate a system
  --   with random coefficients in quad double precision.
  --   The series are computed for the generated system.

    n,d,m,c : natural32 := 0;
    e : integer32 := 0;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use QuadDobl_Complex_Laurentials;

  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    e := integer32(n) - 1;
    QuadDobl_Generate_and_Show(n,d,m,c,e,lp);
    declare
      q : QuadDobl_Complex_Laur_Systems.Laur_Sys(lp'range)
        := Polynomial_to_Laurent_System(lp.all);
      nv : constant integer32 := integer32(Number_of_Unknowns(q(q'first)));
    begin
      QuadDobl_Test(q,q'last,nv);
    end;
  end QuadDobl_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then launches the proper main driver.

    prc : constant character := Prompt_for_Precision;
    ans : character;

  begin
    new_line;
    put("Generate a random polynomial system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      case prc is
        when '0' => Standard_Random_Test;
        when '1' => DoblDobl_Random_Test;
        when '2' => QuadDobl_Random_Test;
        when others => null;
      end case;
    else
      case prc is
        when '0' => Standard_Main;
        when '1' => DoblDobl_Main;
        when '2' => QuadDobl_Main;
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_puiseux;
