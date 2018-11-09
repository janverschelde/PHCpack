with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_SysFun;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_SysFun;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_CSeries_Vector_Functions;
with DoblDobl_CSeries_Vector_Functions;
with QuadDobl_CSeries_Vector_Functions;
with Complex_Series_and_Polynomials_io;
with Regular_Solution_Curves_Series;     use Regular_Solution_Curves_Series;

package body Regular_Newton_Puiseux is

  procedure Tropisms_by_Mixed_Cells
              ( file : in file_type; 
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32 ) is
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
                s : Standard_Complex_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : double_float ) return double_float is

    res : double_float := 0.0;
    x : constant Standard_Complex_Vectors.Vector(s'range)
      := Standard_CSeries_Vector_Functions.Eval(s,w,t);
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
                s : DoblDobl_Complex_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : double_double ) return double_double is

    res : double_double := create(0.0);
    x : constant DoblDobl_Complex_Vectors.Vector(s'range)
      := DoblDobl_CSeries_Vector_Functions.Eval(s,w,t);
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
                s : QuadDobl_Complex_Series_Vectors.Vector;
                w : Standard_Integer_Vectors.Vector;
                t : quad_double ) return quad_double is

    res : quad_double := create(0.0);
    x : constant QuadDobl_Complex_Vectors.Vector(s'range)
      := QuadDobl_CSeries_Vector_Functions.Eval(s,w,t);
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
                s : Standard_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_float ) return double_float is

    sum : double_float := 0.0;
    res : double_float;

  begin
    put_line(file,"The system p :"); put(file,p);
    for k in s'range loop
      put(file,"-> at the series "); put(file,k); put_line(file," : ");
      Complex_Series_and_Polynomials_io.put(file,s(k).all);
      put(file,"with tropism "); put(file,w(k)); new_line(file);
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
                s : Standard_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_float; report : in boolean )
              return double_float is

    sum : double_float := 0.0;
    res : double_float;

  begin
    if report
     then put_line("The system p :"); put(p);
    end if;
    for k in s'range loop
      if report then
        put("-> at the series "); put(k); put_line(" : ");
        Complex_Series_and_Polynomials_io.put(s(k).all);
        put("with tropism "); put(w(k)); new_line;
      end if;
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
                s : DoblDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_double ) return double_double is

    sum : double_double := create(0.0);
    res : double_double;

  begin
    put_line(file,"The system p :"); put(file,p);
    for k in s'range loop
      put(file,"-> at the series "); put(file,k); put_line(file," : ");
      Complex_Series_and_Polynomials_io.put(file,s(k).all);
      put(file,"with tropism "); put(file,w(k)); new_line(file);
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
                s : DoblDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : double_double; report : in boolean )
              return double_double is

    sum : double_double := create(0.0);
    res : double_double;

  begin
    if report
     then put_line("The system p :"); put(p);
    end if;
    for k in s'range loop
      if report then
        put("-> at the series "); put(k); put_line(" : ");
        Complex_Series_and_Polynomials_io.put(s(k).all);
        put("with tropism "); put(w(k)); new_line;
      end if;
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
                s : QuadDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : quad_double ) return quad_double is

    sum : quad_double := create(0.0);
    res : quad_double;

  begin
    put_line(file,"The system p :"); put(file,p);
    for k in s'range loop
      put(file,"-> at the series "); put(file,k); put_line(file," : ");
      Complex_Series_and_Polynomials_io.put(file,s(k).all);
      put(file,"with tropism "); put(file,w(k)); new_line(file);
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
                s : QuadDobl_Complex_Series_VecVecs.VecVec;
                w : Standard_Integer_VecVecs.VecVec;
                t : quad_double; report : in boolean )
              return quad_double is

    sum : quad_double := create(0.0);
    res : quad_double;

  begin
    if report
     then put_line("The system p :"); put(p);
    end if;
    for k in s'range loop
      if report then
        put("-> at the series "); put(k); put_line(" : ");
        Complex_Series_and_Polynomials_io.put(s(k).all);
        put("with tropism "); put(w(k)); new_line;
      end if;
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

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;
    maxdeg : constant integer32 := 16;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : Standard_Complex_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc,mv);
      r : double_float;
    begin
      s := Series(file,p,mcc,mv,maxdeg,nit);
      put_line(file,"Evaluating series at 0.1 ...");
      r := Standard_Residuals(file,p,s,w,0.1);
      put(file,"The reported sum of residuals : ");
      put(file,r,3); new_line(file);
    end;
  end Standard_Test;

  procedure DoblDobl_Test
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;
    maxdeg : constant integer32 := 16;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : DoblDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc,mv);
      t : constant double_double := create(0.1);
      r : double_double;
    begin
      s := Series(file,p,mcc,mv,maxdeg,nit);
      put_line(file,"Evaluating series at 0.1 ...");
      r := DoblDobl_Residuals(file,p,s,w,t);
      put(file,"The reported sum of residuals : ");
      put(file,r,3); new_line(file);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;
    maxdeg : constant integer32 := 16;

  begin
    Tropisms_by_Mixed_Cells(file,sup,mcc,mv);
    declare
      s : QuadDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc,mv);
      t : constant quad_double := create(0.1);
      r : quad_double;
    begin
      s := Series(file,p,mcc,mv,maxdeg,nit);
      put_line(file,"Evaluating series at 0.1 ...");
      r := QuadDobl_Residuals(file,p,s,w,t);
      put(file,"The reported sum of residuals : ");
      put(file,r,3); new_line(file);
    end;
  end QuadDobl_Test;

  procedure Standard_Test 
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;
    maxdeg : constant integer32 := 16;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : Standard_Complex_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc,mv);
      r : double_float;
    begin
      s := Series(p,mcc,mv,maxdeg,nit,report);
      if report
       then put_line("Evaluating the series at 0.1 ...");
      end if;
      r := Standard_Residuals(p,s,w,0.1,report);
      if report
       then put("The reported sum of residuals : "); put(r,3); new_line;
      end if;
    end;
  end Standard_Test;

  procedure DoblDobl_Test 
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;
    maxdeg : constant integer32 := 16;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : DoblDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc,mv);
      t : constant double_double := create(0.1);
      r : double_double;
    begin
      s := Series(p,mcc,mv,maxdeg,nit,report);
      if report
       then put_line("Evaluating the series at 0.1 ...");
      end if;
      r := DoblDobl_Residuals(p,s,w,t,report);
      if report
       then put("The reported sum of residuals : "); put(r,3); new_line;
      end if;
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test 
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is

    sup : Array_of_Lists(p'range) := Create(p);
    mcc : Mixed_Subdivision;
    mv : natural32;
    nit : constant integer32 := 7;
    maxdeg : constant integer32 := 16;

  begin
    Tropisms_by_Mixed_Cells(sup,mcc,mv,report);
    declare
      s : QuadDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
      w : constant Standard_Integer_VecVecs.VecVec := Tropisms(mcc,mv);
      t : constant quad_double := create(0.1);
      r : quad_double;
    begin
      s := Series(p,mcc,mv,maxdeg,nit,report);
      if report
       then put_line("Evaluating the series at 0.1 ...");
      end if;
      r := QuadDobl_Residuals(p,s,w,t,report);
      if report
       then put("The reported sum of residuals : "); put(r,3); new_line;
      end if;
    end;
  end QuadDobl_Test;

  procedure Standard_Read
              ( lp : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                nq,nv : out integer32 ) is

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

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    QuadDobl_Read(lp,nq,nv);
    if nv /= nq+1 
     then put(nv,1); put(" /= "); put(nq,1); put(" + 1");
     else QuadDobl_Test(lp.all,nq,nv);
    end if;
  end QuadDobl_Main;

end Regular_Newton_Puiseux;
