with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Poly_Lists;        use Standard_Complex_Poly_Lists;
with Standard_Complex_Poly_Lists_io;     use Standard_Complex_Poly_Lists_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Complex_Prod_Systems_io;   use Standard_Complex_Prod_Systems_io;
with Lexicographic_Root_Enumeration;     use Lexicographic_Root_Enumeration;
with Set_Structure,Set_Structure_io;
with Main_Set_Structures;
with Standard_Linear_Product_System;
with Random_Product_Start_Systems;
with Standard_Complex_Prod_Planes;

package body Test_Product_Systems is

  procedure Write_Support ( h : in Standard_Complex_Vectors.Vector ) is

    zero : constant Complex_Number := Create(0.0);
    first : boolean := true;

  begin
    put("{");
    for i in 1..h'last loop
      if h(i) /= zero then
        if first then
          first := false;
        else
          put(" ");
        end if;
        put("x"); put(i,1);
      end if;
    end loop;
    put("}");
  end Write_Support;

  function Same_Support
             ( h1,h2 : Standard_Complex_Vectors.Vector)
             return boolean is

    zero : constant Complex_Number := Create(0.0);

  begin
    for i in 1..h1'last loop
      if (h1(i) = zero and h2(i) /= zero)
        or (h1(i) /= zero and h2(i) = zero)
       then return false;
      end if;
    end loop;
    return true;
  end Same_Support;

  function Hyper_Multiplicities return Standard_Natural_VecVecs.VecVec is

    n : constant natural32 := Standard_Linear_Product_System.Dimension;
    res : Standard_Natural_VecVecs.VecVec(1..integer32(n));
    hyp1,hyp2 : Standard_Complex_Vectors.Vector(0..integer32(n));
    nh,m : natural32;

  begin
    for i in 1..n loop
      nh := Standard_Linear_Product_System.Number_of_Hyperplanes(i);
      declare
        vm : Standard_Natural_Vectors.Vector(1..integer32(nh));
      begin
        for j in 1..nh loop -- compute multiplicity of j-th hyperplane
          m := 1;
          hyp1 := Standard_Linear_Product_System.Get_Hyperplane(i,j);
          for k in 1..nh loop
            if k /= j then
              hyp2 := Standard_Linear_Product_System.Get_Hyperplane(i,k);
              if Same_Support(hyp1,hyp2) then
                if k < j
                 then m := 0;       -- occurs already
                 else m := m + 1;   -- increase multiplicity
                end if;
              end if;
            end if;
            exit when (m = 0);
          end loop;
          vm(integer32(j)) := m;
        end loop;
        res(integer32(i)) := new Standard_Natural_Vectors.Vector'(vm);
      end;
    end loop;
    return res;
  end Hyper_Multiplicities;

  procedure Write_Multiplicities ( m : in Standard_Natural_VecVecs.VecVec) is

    n : constant natural32 := Standard_Linear_Product_System.Dimension;

  begin
    for i in 1..integer32(n) loop
      put(i,1); put(" : ");
      for j in 1..m(i)'last loop -- write multiplicity of j-th hyperplane
        put("  m("); put(j,1); put(") = "); put(m(i)(j),1);
      end loop;
      new_line;
    end loop;
  end Write_Multiplicities;

  procedure Test_Read_and_Write is

    n : natural32 := 0;
    p : Prod_Poly;
    ep : Poly;

  begin
    put("Give the number of unknowns : "); get(n);
    Symbol_Table.Init(n);
    put("Give a product of polynomials in "); put(n,1);
    put_line(" unknowns, terminate by semicolon : ");
    get(p);
    put("Your polynomial has "); 
    put(Number_of_Factors(p),1); put_line(" factors :");
    put(p); new_line;
    put("The polynomial written in line format : ");
    put_line(p); new_line;
    ep := Expand(p);
    put_line("The polynomial in its expanded form :");
    put(ep); new_line;
  end Test_Read_and_Write;

  procedure Test_Product_Systems is

    lp : Link_to_Prod_Sys;
    n : integer32 := 0;
    m : natural32 := 0;

  begin
    put_line("Reading a system of product polynomials...");
    get(lp);
    n := lp'last;
    m := Number_of_Unknowns(lp(lp'first));
    put("The system of "); put(n,1);
    put(" equations in "); put(m,1); put_line(" unknowns : ");
    put(natural32(n),m,lp.all);
    declare
      ep : constant Poly_Sys(1..n) := Expand(lp.all);
    begin
      put_line("The expanded polynomial system : "); put(ep);
    end;
  end Test_Product_Systems;

  function Count_All_Solutions_of_Linear_Product_System
             ( output : in boolean ) return natural32 is

    res : natural32 := 0;
    tol : constant double_float := 1.0E-8;
    rc : natural64;

    procedure Report ( s : in Standard_Natural_Vectors.Vector;
                       c : out boolean ) is
    begin
      res := res + 1;
      if output
       then put(s); put(" : + 1 = "); put(res,1); new_line;
      end if;
      c := true;
    end Report;
    procedure Count is 
      new Standard_Linear_Product_System.Enumerate_Solutions(Report);

  begin
    Count(tol,rc);
    put("The number of solutions counted : ");
    Standard_Natural_Numbers_io.put(rc,1); new_line;
    return res;
  end Count_All_Solutions_of_Linear_Product_System;

  procedure Count_All_Solutions_of_Linear_Product_System  is

    r : natural32;
    ans : character;
    timer : Timing_Widget;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y'
     then r := Count_All_Solutions_of_Linear_Product_System(true);
     else r := Count_All_Solutions_of_Linear_Product_System(false);
    end if;
    tstop(timer);
    put("the number of solutions : "); put(r,1); new_line;
    print_times(standard_output,timer,"counting the roots");
  end Count_All_Solutions_of_Linear_Product_System;

  procedure Set_Structure_Linear_Product_Start_System is

    lp,lq : Link_to_Poly_Sys;
    ans : character;

  begin
    put_line("Reading any polynomial system..."); get(lp);
    new_line;
    put_line("MENU to define a supporting set structure : ");
    put_line("  0. let the computer generate a set structure;");
    put_line("  1. give your own set structure.");
    put("Type 0 or 1 to select : "); Ask_Alternative(ans,"01");
    if ans = '0'
     then Random_Product_Start_Systems.Build_Set_Structure(lp.all);
     else Main_Set_Structures.Read_Set_Structure(natural32(lp'last));
    end if;
    new_line;
    put_line("The set structure :");
    Set_Structure_io.put;
    Standard_Linear_Product_System.Init(natural32(lp'last));
    Random_Product_Start_Systems.Build_Random_Product_System
      (natural32(lp'last));
    put("Do you see the expanded system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      lq := new Poly_Sys'(Standard_Linear_Product_System.Polynomial_System);
      put_line("A random linear-product start system : "); put_line(lq.all);
    end if;
    declare
      rq : constant Prod_Sys(lp'range)
         := Standard_Complex_Prod_Planes.Create;
      file : file_type;
    begin
      put_line("The system in product form : ");
      put_line(natural32(rq'last),rq);
      put("Do you want to write the linear-product system to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of the output file...");
        Read_Name_and_Create_File(file);
        put_line(file,natural32(rq'last),rq);
      end if;
      new_line;
      put("Do you want to count the solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' 
       then Count_All_Solutions_of_Linear_Product_System;
      end if;
    end;
  end Set_Structure_Linear_Product_Start_System;

  function Residual ( A : Standard_Complex_Matrices.Matrix;
                      b,x : Standard_Complex_Vectors.Vector )
                    return double_float is

    use Standard_Complex_Vectors,Standard_Complex_Matrices;

    y : constant Vector(A'range(1)) := b - A*x;

  begin
    return Max_Norm(y);
  end Residual;

  function Lex_Count_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector ) return natural32 is

    n : constant natural32 := natural32(d'last);
    cont : boolean := true;
    cnt : natural32 := 0;
    v : Standard_Complex_Vectors.Vector(d'range);
    acc : Standard_Natural_Vectors.Vector(d'range);
    fail : boolean;
    tol : constant double_float := 1.0E-10;
    rcond : double_float;

    procedure Count ( acc : in Standard_Natural_Vectors.Vector;
                      continue : out boolean ) is
    begin
      Standard_Linear_Product_System.Solve(acc,tol,rcond,fail,v);
      if not fail then
        cnt := cnt + 1;
      end if;
      continue := true;
    end Count;
    procedure Enum is new Lexicographic_Enumeration(Count);

  begin
    Enum(1,n,d,acc,cont);
    return cnt;
  end Lex_Count_Solutions_of_Linear_Product_System;

  function Multiplicity
              ( m : Standard_Natural_VecVecs.VecVec;
                a : Standard_Natural_Vectors.Vector ) return natural32 is

    res : natural32 := 1;
    mlt : natural32;

  begin
    for i in a'range loop
      mlt := m(i)(integer32(a(i)));
      if mlt = 0
       then return 0;
       else res := res*mlt;
      end if;
    end loop;
    return res;
  end Multiplicity;

  function Lex_Count_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector;
                m : Standard_Natural_VecVecs.VecVec ) return natural32 is

    n : constant natural32 := natural32(d'last);
    cont : boolean := true;
    cnt : natural32 := 0;
    v : Standard_Complex_Vectors.Vector(d'range);
    acc : Standard_Natural_Vectors.Vector(d'range);
    fail : boolean;
    tol : constant double_float := 1.0E-10;
    rcond : double_float;

    procedure Count ( acc : in Standard_Natural_Vectors.Vector;
                      continue : out boolean ) is

      ma : constant natural32 := Multiplicity(m,acc);

    begin
      if ma /= 0 then
        Standard_Linear_Product_System.Solve(acc,tol,rcond,fail,v);
        if not fail then
          cnt := cnt + ma;
        end if;
      end if;
      continue := true;
    end Count;
    procedure Enum is new Lexicographic_Enumeration(Count);

  begin
    Enum(1,n,d,acc,cont);
    return cnt;
  end Lex_Count_Solutions_of_Linear_Product_System;

  procedure Lex_Enumerate_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector ) is

    n : constant natural32 := natural32(d'last);
    cont : boolean := true;
    cnt : natural32 := 0;
    ind : natural32 := 0;
    acc : Standard_Natural_Vectors.Vector(d'range);
    cp : constant Standard_Natural_Vectors.Vector
       := Consecutive_Products(d);
    A : Standard_Complex_Matrices.Matrix(d'range,d'range);
    b,v : Standard_Complex_Vectors.Vector(d'range);
    fail : boolean;
    tol : constant double_float := 1.0E-10;
    rcond,nrm : double_float;
    timer : Timing_Widget;

    procedure Write ( acc : in Standard_Natural_Vectors.Vector;
                      continue : out boolean ) is
    begin
      ind := ind + 1;
      Standard_Linear_Product_System.Linear_System(acc,fail,A,b);
      if not fail then
        Standard_Linear_Product_System.Solve(acc,tol,rcond,fail,v);
        if not fail then
          cnt := cnt + 1; put(cnt,3); put(" : ");
          put(acc); put(" : rco ="); 
          put(rcond,3); put(" : res =");
          nrm := Residual(A,b,v);
          put(nrm,3);
        else
          put(" -- : "); put(acc);
          put(" : rco ="); put(rcond,3);
          put(" : --- ");
        end if;
        put(" : map = "); put(Root_Map(n,ind,d,cp)); new_line;
      end if;
      continue := true;
    end Write;
    procedure Enum is new Lexicographic_Enumeration(Write);

  begin
    tstart(timer);
    Enum(1,n,d,acc,cont);
    tstop(timer);
    put("Counted "); put(cnt,1); put_line(" solutions.");
    new_line;
    print_times(standard_output,timer,"Lexicographic Root Counter");
  end Lex_Enumerate_Solutions_of_Linear_Product_System;

  procedure Enumerate_All_Solutions_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector ) is

    n : constant natural32 := natural32(d'last);
   -- cont : boolean := true;
    cnt : natural32 := 0;
    ind : natural32 := 0;
    cp : constant Standard_Natural_Vectors.Vector
       := Consecutive_Products(d);
    A : Standard_Complex_Matrices.Matrix(d'range,d'range);
    b,v : Standard_Complex_Vectors.Vector(d'range);
    fail : boolean;
    tol : constant double_float := 1.0E-10;
    rcond,nrm : double_float;
    rc : natural64;
    timer : Timing_Widget;

    procedure Write ( acc : in Standard_Natural_Vectors.Vector;
                      continue : out boolean ) is
    begin
      ind := ind + 1;
      Standard_Linear_Product_System.Linear_System(acc,fail,A,b);
      if not fail then
        Standard_Linear_Product_System.Solve(acc,tol,rcond,fail,v);
        if not fail then
          cnt := cnt + 1; 
          Standard_Natural_Numbers_io.put(cnt,3); put(" : ");
          put(acc); put(" : rco ="); 
          put(rcond,3); put(" : res =");
          nrm := Residual(A,b,v);
          put(nrm,3);
        else
          put(" -- : "); put(acc);
          put(" : rco ="); put(rcond,3);
          put(" : --- ");
        end if;
        put(" : map = "); put(Root_Map(n,ind,d,cp)); new_line;
      end if;
      continue := true;
    end Write;
    procedure Enum is 
      new Standard_Linear_Product_System.Enumerate_Solutions(Write);

  begin
    tstart(timer);
    Enum(tol,rc);
    tstop(timer);
    put("Counted ");
    Standard_Natural_Numbers_io.put(cnt,1);
    put_line(" solutions.");
    put("rc = "); put(rc,1); new_line;
    new_line;
    print_times(standard_output,timer,"enumerating all solutions");
  end Enumerate_All_Solutions_of_Linear_Product_System;

  procedure Test_Solution_of_Linear_Product_System
              ( d : Standard_Natural_Vectors.Vector ) is

    s : Standard_Natural_Vectors.Vector(d'range);
    A : Standard_Complex_Matrices.Matrix(d'range,d'range);
    b : Standard_Complex_Vectors.Vector(d'range);
    tol : constant double_float := 10.0E-10;
    ans : character;
    rcond : double_float;
    fail : boolean;
    v : Standard_Complex_Vectors.Vector(d'range);

  begin
    loop
      for i in s'range loop
        put("Give index for equation "); put(i,1); put(" : ");
        s(i) := 0; get(s(i));
      end loop;
      Standard_Linear_Product_System.Linear_System(s,fail,A,b);
      if fail then
        put_line("Retrieval of linear system failed.");
      else
        put_line("The matrix A : "); put(A);
        put_line("The vector b : "); put_line(b);
        Standard_Linear_Product_System.Solve(s,tol,rcond,fail,v);
        put("estimate for condition number of linear system : ");
        put(rcond,3); new_line;
        if fail then
          put_line("Solution of linear system failed.");
        else
          put_line("Solution of linear system succeeded.");
          put_line("The solution vector : "); put_line(v);
          put("The residual :"); put(Residual(A,b,v),3); new_line;
        end if;
      end if;
      put("Do you want more test ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Solution_of_Linear_Product_System;

  procedure Count_Roots ( q : in Prod_Sys ) is

    d : constant Standard_Natural_Vectors.Vector(q'range)
      := Standard_Complex_Prod_Planes.Degrees(q);
    cnt : natural32;
    ans : character;
    timer : Timing_Widget;

  begin
    put("The degrees of the equations in the system : "); put(d); new_line;
    new_line;
    put_line("MENU to solve a linear-product start system : ");
    put_line("  1. enumerate and count all solutions without multiplicities;");
    put_line("  2.                                   with multiplicities;");
    put_line("  3. interactive test of selected linear systems.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    if ans = '1' then
      put("Do you want intermediate output ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
       -- Lex_Enumerate_Solutions_of_Linear_Product_System(d);
        Enumerate_All_Solutions_of_Linear_Product_System(d);
      else
        tstart(timer);
        cnt := Lex_Count_Solutions_of_Linear_Product_System(d);
       -- cnt := Count_All_Solutions_of_Linear_Product_System(d);
        tstop(timer);
        put("Counted "); put(cnt,1); put_line(" solutions.");
        new_line;
        print_times(standard_output,timer,
                    "Set Structure Bound without multiplicities");
      end if;
    elsif ans = '2' then
      declare
        m : constant Standard_Natural_VecVecs.VecVec(d'range)
          := Hyper_Multiplicities;
      begin
        Write_Multiplicities(m);
        tstart(timer);
        cnt := Lex_Count_Solutions_of_Linear_Product_System(d,m);
        tstop(timer);
        put("Set Structure Bound : "); put(cnt,1); new_line;
        new_line;
        print_times(standard_output,timer,
                    "Set Structure Bound with multiplicities");
      end;
    else
      Test_Solution_of_Linear_Product_System(d);
    end if;
  end Count_Roots;

  procedure Solve_Linear_Product_Start_System is

    file : file_type;
    lq : Link_to_Prod_Sys;
    fail : boolean;
    ans : character;

  begin
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(file);
    get(file,lq);
    new_line;
    put("Do you want to see the system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("The system read from file : "); put_line(lq.all);
    end if;
    Standard_Complex_Prod_Planes.Store(lq.all,fail);
    if fail then
      put_line("Storing the system as a linear-product system failed!");
    else
      put_line("The system is stored as a linear-product system.");
      Count_Roots(lq.all);
    end if;
  end Solve_Linear_Product_Start_System;

  procedure Interactive_Inverted_Enumerator ( q : in Prod_Sys ) is

    s,d : Standard_Natural_Vectors.Vector(q'range);
    tol : constant double_float := 1.0E-8;
    fail : boolean;
    cnt : natural32 := 1;
    ans : character;

  begin
    Standard_Linear_Product_System.Get_First(standard_output,tol,s,fail);
    if fail then
      put_line("The given linear-product system has no solution.");
    else
      put("First root is at "); put(s); new_line;
      d := Standard_Complex_Prod_Planes.Degrees(q);
      loop
        cnt := cnt + 1;
        Standard_Linear_Product_System.Get_Next(standard_output,tol,d,s,fail);
        exit when fail;
        put("Root "); put(cnt,1); put(" : "); put(s);
        new_line;
        put("Continue enumeration ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end loop;
    end if;
    Standard_Linear_Product_System.Get_Clear;
  end Interactive_Inverted_Enumerator;

  function Create ( v : Standard_Complex_Vectors.Vector;
                    rcond : double_float ) return Solution is

    res : Solution(v'last);

  begin
    res.t := Create(0.0);
    res.m := 1;
    res.v := v;
    res.err := 0.0;
    res.rco := rcond;
    res.res := 0.0;
    return res;
  end Create;

  procedure Solve_by_Inverted_Enumerator
               ( file : in file_type;
                 q : in Prod_Sys; rc : in natural32 ) is

    s,d : Standard_Natural_Vectors.Vector(q'range);
    tol : constant double_float := 1.0E-8;
    fail : boolean;
    rcond : double_float;
    v : Standard_Complex_Vectors.Vector(q'range);
    n : constant integer32 := q'last;
    sol : Solution(n);
    cnt : natural32 := 0;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    Write_First(file,rc,natural32(n));
    Standard_Linear_Product_System.Get_First(tol,s,fail);
    if fail then
      put_line("The given linear-product system has no solution.");
    else
      d := Standard_Complex_Prod_Planes.Degrees(q);
      loop
        Standard_Linear_Product_System.Solve(s,tol,rcond,fail,v);
        put("root "); put(cnt+1,1); put(" at"); put(s);
        put(" has rcond = "); put(rcond,3);
        if rcond > tol
         then put_line(" okay");
         else put_line(" BUG !!!"); exit;
        end if;
        sol := Create(v,rcond);
        Write_Next(file,cnt,sol);
        Standard_Linear_Product_System.Get_Next(tol,d,s,fail);
        exit when fail;
      end loop;
    end if;
  end Solve_by_Inverted_Enumerator;

  procedure Run_Inverted_Enumerator
              ( file : in out file_type; q : in Prod_Sys ) is

    s,d : Standard_Natural_Vectors.Vector(q'range);
    tol : constant double_float := 1.0E-8;
    fail : boolean;
    cnt : natural32;
    ans : character;

  begin
    Standard_Linear_Product_System.Get_First(tol,s,fail);
    if fail then
      put_line("The given linear-product system has no solution.");
    else
      cnt := 1;
      put("Root "); put(cnt,3); put(" : "); put(s); new_line;
      d := Standard_Complex_Prod_Planes.Degrees(q);
      loop
        Standard_Linear_Product_System.Get_Next(tol,d,s,fail);
        exit when fail;
        cnt := cnt + 1;
        put("Root "); put(cnt,3); put(" : "); put(s); new_line;
      end loop;
    end if;
    Standard_Linear_Product_System.Get_Clear;
    new_line;
    put("Append the "); put(cnt,1);
    put(" solutions to the input file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Reset(file,append_file);
      Solve_by_Inverted_Enumerator(file,q,cnt);
    end if;
  end Run_Inverted_Enumerator;

  procedure Test_Inverted_Enumerator is

    file : file_type;
    lq : Link_to_Prod_Sys;
    fail : boolean;
    ans : character;

  begin
    put_line("Reading the file name for a product system...");
    Read_Name_and_Open_File(file);
    get(file,lq);
    new_line;
    put("Do you want to see the system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("The system read from file : "); put_line(lq.all);
    end if;
    Standard_Complex_Prod_Planes.Store(lq.all,fail);
    if fail then
      put_line("Storing the system as a linear-product system failed!");
    else
      put_line("The system is stored as a linear-product system.");
      new_line;
      put("Do you want interactive incremental version ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Interactive_Inverted_Enumerator(lq.all);
       else Run_Inverted_Enumerator(file,lq.all);
      end if;
    end if;
  end Test_Inverted_Enumerator;

  procedure Write_Supports ( nf : in Standard_Natural_Vectors.Vector ) is

    hyp : Standard_Complex_Vectors.Link_to_Vector;
    n : constant natural32 := Standard_Linear_Product_System.Dimension;
    m : Standard_Natural_VecVecs.VecVec(1..integer32(n));
    zero : constant Complex_Number := Create(0.0);

  begin
    Set_Structure.Init(nf);
    for i in 1..n loop
      for j in 1..Standard_Linear_Product_System.Number_of_Hyperplanes(i) loop
        hyp := Standard_Linear_Product_System.Get_Hyperplane(i,j);
        Write_Support(hyp.all);
        for k in 1..hyp'last loop
          if hyp(k) /= zero 
           then Set_Structure.Add(i,j,natural32(k));
          end if;
        end loop;
      end loop;
      new_line;
    end loop;
   -- Set_Structure_io.put;
    m := Hyper_Multiplicities;
    Write_Multiplicities(m);
  end Write_Supports;

  procedure Support_of_Linear_Product_System is

    file : file_type;
    lq : Link_to_Prod_Sys;
    fail : boolean;
    ans : character;

  begin
    put_line("Reading the file name for a product system...");
    Read_Name_and_Open_File(file);
    get(file,lq);
    new_line;
    put("Do you want to see the system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("The system read from file : "); put_line(lq.all);
    end if;
    Standard_Complex_Prod_Planes.Store(lq.all,fail);
    if fail then
      put_line("The system is not a linear-product system!");
    else
      declare
        nf : constant Standard_Natural_Vectors.Vector(lq'range)
           := Standard_Complex_Prod_Planes.Degrees(lq.all);
      begin
        put("#factors in the system : "); put(nf); new_line;
        Write_Supports(nf);
      end;
    end if;
  end Support_of_Linear_Product_System;

  function Difference ( x,y : Standard_Complex_Matrices.Matrix )
                      return double_float is

    res : double_float := 0.0;
    xmy : Complex_Number;
    v : double_float;

  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        xmy := x(i,j) - y(i,j);
        v := AbsVal(xmy);
        res := res + v;
      end loop;
    end loop;
    return res;
  end Difference;

  procedure Evaluate_Linear_Product_System is

    n : constant natural32 := Standard_Linear_Product_System.Dimension;
    p : Poly_Sys(1..integer32(n));
    jf : Jaco_Mat(p'range,p'range);
    x,y1,y2,df : Standard_Complex_Vectors.Vector(p'range);
    m1,m2 : Standard_Complex_Matrices.Matrix(p'range,p'range);
    ndf : double_float;
    ans : character;

    use Standard_Complex_Vectors;

  begin
    put_line("expanding the linear-product system ...");
    p := Standard_Linear_Product_System.Polynomial_System;
    jf := Create(p);
    loop
      x := Random_Vector(1,integer32(n));
      y1 := Eval(p,x);
      put_line("Evaluated expanded system at a random point : ");
      put_line(y1);
      y2 := Standard_Complex_Prod_Planes.Eval(x);
      put_line("Evaluated linear-product system at the same point : ");
      put_line(y2);
      df := y1 - y2; ndf := Norm2(df);
      put("difference between the values : "); put(ndf,3); new_line;
      put("continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      m1 := Eval(jf,x);
      put_line("Jacobian evaluated from expanded system at random point : ");
      put(m1);
      m2 := Standard_Complex_Prod_Planes.Jacobian(x);
      put_line("Jacobian evaluated from linear-product system : ");
      put(m2);
      ndf := Difference(m1,m2);
      put("difference between the values : "); put(ndf,3); new_line;
      put("Test another point ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Evaluate_Linear_Product_System;

  procedure Test_Evaluation is

    file : file_type;
    lq : Link_to_Prod_Sys;
    fail : boolean;
    ans : character;

  begin
    put_line("Reading the file name for a product system...");
    Read_Name_and_Open_File(file);
    get(file,lq);
    new_line;
    put("Do you want to see the system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("The system read from file : "); put_line(lq.all);
    end if;
    Standard_Complex_Prod_Planes.Store(lq.all,fail);
    if fail then
      put_line("The system is not a linear-product system!");
    else
      Evaluate_Linear_Product_System;
    end if;
  end Test_Evaluation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test operations on products of polynomials :");
    put_line("  1. test read and write of a product polynomial;");
    put_line("  2. test operations on product systems;");
    put_line("  3. set structure linear-product start system;");
    put_line("  4. functional enumerator to solve a linear-product system;");
    put_line("  5. solve a linear-product system with inverted enumerator;");
    put_line("  6. test supporting set structure of linear-product system;");
    put_line("  7. evaluation of linear-product systems.");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to make your choice : ");
    Ask_Alternative(ans,"1234567");
    new_line;
    case ans is
      when '1' => Test_Read_and_Write;
      when '2' => Test_Product_Systems;
      when '3' => Set_Structure_Linear_Product_Start_System;
      when '4' => Solve_Linear_Product_Start_System;
      when '5' => Test_Inverted_Enumerator;
      when '6' => Support_of_Linear_Product_System;
      when '7' => Test_Evaluation;
      when others => null;
    end case;
  end Main;

end Test_Product_Systems;
