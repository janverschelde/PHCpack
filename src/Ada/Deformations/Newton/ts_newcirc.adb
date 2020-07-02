with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Vector_Splitters;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Standard_Coefficient_Circuits;      use Standard_Coefficient_Circuits;
with Standard_Circuit_Makers;            use Standard_Circuit_Makers;
with Standard_Newton_Circuits;
with Standard_Inlined_Newton_Circuits;
with Standard_Refiner_Circuits;

procedure ts_newcirc is

-- DESCRIPTION :
--   Test the development of Newton's method on coefficient circuits.

  procedure LU_Newton_Steps
              ( s : in Link_to_System;
                v : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton steps on the system s starting at the vector v.
 
  -- REQUIRED : The system is square.

    use Standard_Newton_Circuits;

    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    res,rco,err,mixres : double_float;
    ans : character;
    condition,nomixres : boolean;
    abscfs : Link_to_System;
    radv : Standard_Complex_Vectors.Vector(v'range);

  begin
    put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
    condition := (ans = 'y');
    put_line("The vector v :"); put_line(v);
    new_line;
    put("Compute mixed residuals ? (y/n) "); Ask_Yes_or_No(ans);
    nomixres := (ans = 'n');
    if ans = 'y' then
      abscfs := Copy(s);
      AbsVal(abscfs);
    end if;
    loop
      if nomixres then
        if condition
         then LU_Newton_Step(standard_output,s,v,xr,xi,ipvt,res,rco,err,true);
         else LU_Newton_Step(standard_output,s,v,xr,xi,ipvt,info,res,err,true);
        end if;
      else
        if condition then
          LU_Newton_Step(standard_output,s,abscfs,v,radv,xr,xi,ipvt,
                         res,rco,err,mixres,true);
        else
          LU_Newton_Step(standard_output,s,abscfs,v,radv,xr,xi,ipvt,info,
                         res,err,mixres,true);
        end if;
      end if;
      put_line("The vector v :"); put_line(v);
      put("Another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end LU_Newton_Steps;

  procedure Inlined_LU_Newton_Steps
              ( s : in Link_to_System;
                v : in out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton steps on the system s starting at the vector v,
  --   using the inlined LU factorizations.
 
  -- REQUIRED : The system is square.

    use Standard_Inlined_Newton_Circuits;

    ipvt : Standard_Integer_Vectors.Vector(1..s.dim);
    info : integer32;
    vxr : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(v'range)
        := (v'range => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    res,rco,err,mixres : double_float;
    ans : character;
    condition,nomixres : boolean;
    abscfs : Link_to_System;
    radv : Standard_Complex_Vectors.Vector(v'range);

  begin
    put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
    condition := (ans = 'y');
    put_line("The vector v :"); put_line(v);
    new_line;
    put("Compute mixed residuals ? (y/n) "); Ask_Yes_or_No(ans);
    nomixres := (ans = 'n');
    if ans = 'y' then
      abscfs := Copy(s);
      AbsVal(abscfs);
    end if;
    loop
      if nomixres then
        if condition
         then LU_Newton_Step(standard_output,s,v,xr,xi,ipvt,res,rco,err,true);
         else LU_Newton_Step(standard_output,s,v,xr,xi,ipvt,info,res,err,true);
        end if;
      else
        if condition then
          LU_Newton_Step(standard_output,s,abscfs,v,radv,xr,xi,
                         ipvt,res,rco,err,mixres,true);
        else
          LU_Newton_Step(standard_output,s,abscfs,v,radv,xr,xi,
                         ipvt,info,res,err,mixres,true);
        end if;
      end if;
      put_line("The vector v :"); put_line(v);
      put("Another step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
  end Inlined_LU_Newton_Steps;

  procedure Interactive_Run
              ( s : in Link_to_System; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Interactive application of Newton's method on the system of
  --   coefficient circuits s, starting at the solutions in sols.
  --   Pauses frequently, prompting user to continue.

    ptr : Solution_List := sols;
    ls : Link_to_Solution;
    ans : character;
    cnt : integer32 := 0;
    inlined : boolean;

  begin
    put("Inlined linear solver ? (y/n) "); Ask_Yes_or_No(ans);
    inlined := (ans = 'y');
    while not Is_Null(ptr) loop
      ls := Head_Of(ptr); cnt := cnt + 1;
      put("Running Newton's method on solution ");
      put(cnt,1); put_line(" ...");
      if inlined
       then Inlined_LU_Newton_Steps(s,ls.v);
       else LU_Newton_Steps(s,ls.v);
      end if;
      put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      ptr := Tail_Of(ptr);
    end loop;
  end Interactive_Run;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions,
  --   makes a coefficient circuit and then runs Newton's method.

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    nbq,len,dim : integer32 := 0;
    s : Link_to_System;
    ans : character;
    file : file_type;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(p,sols,"THE SOLUTIONS");
    nbq := p'last;
    len := integer32(Length_Of(sols));
    if len = 0 then
      put_line("No solutions found on file.");
    else
      dim := Head_Of(sols).n;
      new_line;
      put("Read system of "); put(nbq,1); put(" polynomials and ");
      put(len,1); put(" solutions in dimension "); put(dim,1); put_line(".");
      s := Make_Coefficient_System(p);
      new_line;
      put("Interactive run ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Interactive_Run(s,sols);
      else
        put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'n' then
          Standard_Refiner_Circuits.Run(s,sols);
        else
          new_line;
          put_line("Reading the name of the output file ...");
          Read_Name_and_Create_File(file);
          Standard_Refiner_Circuits.Run(file,s,sols); --,len);
          Close(file);
        end if;
      end if;
    end if;
  end Main;

begin
  Main;
end ts_newcirc;
