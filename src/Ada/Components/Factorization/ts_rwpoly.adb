with text_io,integer_io;                 use text_io,integer_io;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use standard_Complex_Solutions_io;
with Standard_Random_Polynomials;        use Standard_Random_Polynomials;
with Rewrite_Polynomials;                use Rewrite_Polynomials;

procedure ts_rwpoly is

-- DESCRIPTION :
--   This is a test program to rewrite polynomials of high degree into
--   polynomials of modest degree, at the expense of extra variables.

  procedure Test_Binary ( d : in natural32;
                          deco : in Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION : 
  --   Tests whether the binary decomposition of d matches the value.

    val : natural32 := deco(0);
    acc : natural32 := 1;
    len : natural32 := Binary_Length(d);

  begin
    put("Length : "); put(len,1);
    if len = natural32(deco'last)
     then put(" okay ");
     else put(" BUG! ");
    end if;
    for i in 1..deco'last loop
      acc := acc*2;
      if deco(i) = 1
       then val := val + acc;
      end if;
    end loop;
    put(d,1); 
    if d = val
     then put(" = "); put(val,1); put_line("  okay");
     else put(" <> "); put(val,1); put_line("  bug!!!");
    end if;
  end Test_Binary;

  procedure Binary_Decomposition is

  -- DESCRIPTION :
  --   Interactive test on computing binary decompositions of
  --   natural numbers.

    d : natural32 := 0;
    deco : Standard_Natural_Vectors.Link_to_Vector;
    ans : character;

  begin
    new_line;
    loop
      put("Give a degree : "); get(d);
      if d > 0 then
        put("The binary decomposition of "); put(d,1);
        put_line(" :");
        Binary(0,d,deco); new_line;
        put("The decomposition vector : "); put(deco.all); new_line;
        Test_Binary(d,deco.all);
      end if;
      declare
        bindeco : constant Standard_Natural_Vectors.Vector := Binary(d);
      begin
        put("The decomposition vector : "); put(bindeco); new_line;
        Test_Binary(d,bindeco);
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
      Standard_Natural_Vectors.Clear(deco);
    end loop;
  end Binary_Decomposition;

  procedure Ask_to_Save ( sys : in Poly_Sys ) is

    ans : character;
    file : file_type;

  begin
    put_line("The rewritten system : ");
    put(sys(1..sys'last-1));
    put_line(sys(sys'last));
    new_line;
    put("Do you wish to save the system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(file);
      put(file,natural32(sys'last),1);
      new_line(file);
      put(file,sys(1..sys'last-1));
      put_line(file,sys(sys'last));
    else
      new_line;
      put_line("Bye bye.");
    end if;
  end Ask_to_Save;

  procedure Ask_to_Save ( sys : in Poly_Sys; sols : in Solution_List ) is

    ans : character;
    file : file_type;

  begin
    put_line("The rewritten system : ");
    put(sys(1..sys'last-1));
    put_line(sys(sys'last));
    new_line;
    put("Do you wish to save the system on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(file);
      put(file,natural32(sys'last),1);
      new_line(file);
      put(file,sys(1..sys'last-1));
      put_line(file,sys(sys'last));
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    else
      new_line;
      put_line("Bye bye.");
    end if;
  end Ask_to_Save;

  procedure Rewrite_Given_Polynomial is

    n,len : natural32 := 0;
    p : Poly;

  begin
    new_line;
    put_line("Give the number of unknowns, followed by polynomial : ");
    get(n,p);
    put("Your polynomial : "); put(p); new_line;
    put("The symbol used : "); put(Symbol_Table.get(1)); new_line;
    len := Binary_Length(natural32(Degree(p)));
    Enlarge_Symbol_Table(len+1,Symbol_Table.Get(1));
    Ask_to_Save(Rewrite_Univariate_Polynomial(p));
  end Rewrite_Given_Polynomial;

  function Random_Polynomial ( n,d : natural32 ) return Poly is

  -- DESCRIPTION :
  --   Returns a random univariate polynomial of degree d with n terms.

    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..1 => d);
    t.cf := Random1;
    res := Create(t);
    t.dg(1) := 0;
    t.cf := Random1;
    Add(res,t);
    for i in 3..n loop
      t.cf := Random1;
      t.dg(1) := natural32(Random(0,integer32(d)));
      Add(res,t);
    end loop;
    Clear(t);
    return res;
  end Random_Polynomial;

  procedure Initialize_One_Symbol is

  -- DESCRIPTION :
  --   Initializes the symbol table with the first symbol, which is
  --   used to denote the original variable of the polynomial.

  begin
    Symbol_Table.Init(1);
    declare
      xsb : Symbol;
    begin
      xsb := (xsb'range => ' ');
      xsb(1) := 'x';
      Symbol_Table.Add(xsb);
    end;
  end Initialize_One_Symbol;

  procedure Rewrite_Random_Polynomial is

    n,d : natural32 := 0;
    p : Poly;

  begin
    new_line;
    put("Give the number of terms : "); get(n);
    put("Give the degree of the polynomial : "); get(d);
    p := Random_Polynomial(n,d);
    Initialize_One_Symbol;
    put_line("A random polynomial : "); put_line(p);
    Enlarge_Symbol_Table(Binary_Length(d)+1,Symbol_Table.Get(1));
    Ask_to_Save(Rewrite_Univariate_Polynomial(p));
  end Rewrite_Random_Polynomial;

  procedure Rewrite_Start_Polynomial is

  -- DESCRIPTION :
  --   Rewrites the start polynomial x^d - 1 = 0 into a system.

    d,n : natural32 := 0;
    t : Term;
    p : Poly;

  begin
    new_line;
    put("Give the degree of the polynomial : "); get(d);
    t.dg := new Standard_Natural_Vectors.Vector'(1..1 => d);
    t.cf := Create(1.0);
    p := Create(t);
    t.dg(1) := 0;
    Sub(p,t);
    Initialize_One_Symbol;
    put_line("The start polynomial : "); put_line(p);
    n := Binary_Length(d) + 1;
    Enlarge_Symbol_Table(n,Symbol_Table.Get(1));
    declare
      sys : constant Poly_Sys := Rewrite_Univariate_Polynomial(p);
      one : Complex_Number := Create(1.0);
      sol : Solution(sys'length);
      sols,sols_last : Solution_List;
    begin
      sol.m := 1;
      sol.t := Create(0.0);
      sol.err := 0.0;
      sol.rco := 1.0;
      sol.res := 0.0;
      for i in 1..d loop
        sol.v(sol.v'first) := Root(one,d,i);
        for i in sol.v'first+1..sol.v'last loop
          sol.v(i) := sol.v(i-1)*sol.v(i-1);
        end loop;
        put_line("The function value at the solution :");
        put_line(Eval(sys,sol.v));
        Append(sols,sols_last,sol);
      end loop;
      Ask_to_Save(sys,sols);
    end;
  end Rewrite_Start_Polynomial;

  procedure Check_Multivariate_Rewrite
               ( n : in natural32; p : in Poly;
                 nvr : in Standard_Natural_Vectors.Vector;
                 tnv : in natural32; rp : in Poly ) is

  -- DESCRIPTION :
  --   Performs a random test on the rewritten polynomial,
  --   evaluating the original polynomial and the rewritten
  --   polynomial at a random point.

    x1 : Standard_Complex_Vectors.Vector(1..integer32(n))
       := Random_Vector(1,integer32(n));
    x2 : Standard_Complex_Vectors.Vector(1..integer32(tnv));
    y1 : Complex_Number := Eval(p,x1);
    y2 : Complex_Number;
    ind : integer32 := 0;

  begin
    put_line("Checking the rewrite at a random point...");
    for i in 1..integer32(n) loop
      ind := ind + 1;
      x2(ind) := x1(i);
      for j in 1..nvr(i)-1 loop
        ind := ind + 1;
        x2(ind) := x2(ind-1)*x2(ind-1);
      end loop;
    end loop;
    y2 := Eval(rp,x2);
    put("  y1 : "); put(y1); new_line;
    put("  y2 : "); put(y2); new_line;
  end Check_Multivariate_Rewrite; 

  procedure Rewrite_Multivariate_Polynomial is

    d,n,m : natural32;
    p : Poly;

  begin
    new_line;
    put_line("Generation of a random polynomial : ");
    put("  Give the degree : "); get(d);
    put("  Give the number of variables : "); get(n);
    put("  Give the number of monomials : "); get(m);
    p := Random_Sparse_Poly(n,d,m,0);
    put_line("A random polynomial :"); put_line(p);
    declare
      deg : constant Standard_Natural_Vectors.Vector
          := Multi_Degrees(p);
      nvr : Standard_Natural_Vectors.Vector(deg'range);
      tnv : natural32;
    begin
      put("The multi-degrees : "); put(deg); new_line;
      Number_of_Variables(deg,nvr,tnv);
      Define_Symbol_Table(tnv,nvr);
      declare
        sys : Poly_Sys(1..integer32(tnv));
      begin
        Telescope(sys,tnv,nvr);
        put_line("The telescope :"); put(sys(1..integer32(tnv-n)));
        sys(integer32(tnv-n+1)) := Rewrite_Multivariate_Poly(tnv,p,nvr);
        put_line("The rewritten polynomial :");
        put_line(sys(integer32(tnv-n+1)));
        Check_Multivariate_Rewrite(n,p,nvr,tnv,sys(integer32(tnv-n+1)));
      end;
    end;
    put_line(Rewrite_Multivariate_Polynomial(p));
  end Rewrite_Multivariate_Polynomial;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Rewriting high degree polynomials into lower degree systems.");
    new_line;
    put_line("MENU for interactive testing : ");
    put_line("  1. test binary decomposition of natural numbers;");
    put_line("  2. rewrite a given polynomial in one variable;");
    put_line("  3. rewrite a random polynomial in one variable;");
    put_line("  4. rewrite univariate start polynomial into start system;");
    put_line("  5. rewrite a random multivariate polynomial.");
    put("Type 1, 2, 3, 4, or 5 to select : ");
    Ask_Alternative(ans,"12345");
    case ans is
      when '1' => Binary_Decomposition;
      when '2' => Rewrite_Given_Polynomial;
      when '3' => Rewrite_Random_Polynomial;
      when '4' => Rewrite_Start_Polynomial;
      when '5' => Rewrite_Multivariate_Polynomial;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_rwpoly;
