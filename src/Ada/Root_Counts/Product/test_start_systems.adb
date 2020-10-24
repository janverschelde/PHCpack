with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Lexicographic_Root_Enumeration;     use Lexicographic_Root_Enumeration;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;

package body Test_Start_Systems is

  procedure Enumerate_and_Solve
              ( d : in Standard_Natural_Vectors.Vector;
                c : in Standard_Complex_Vectors.Vector ) is

    n : constant natural32 := natural32(d'last);
    cont : boolean := true;
    cnt : natural32 := 0;
    acc : Standard_Natural_Vectors.Vector(d'range);
    cp : constant Standard_Natural_Vectors.Vector
       := Consecutive_Products(d);
    r,y : Standard_Complex_Vectors.Vector(c'range);
    nrm : double_float;

    procedure Write ( acc : in Standard_Natural_Vectors.Vector;
                      continue : out boolean ) is
    begin
      r := Root(d,acc,c);
      y := Eval(d,c,r);
      nrm := Max_Norm(y);
      cnt := cnt + 1;
      put(cnt,3); put(" : ");
      put(acc); put(" : "); put(nrm);
      put(" : map = "); put(Root_Map(n,cnt,d,cp)); new_line;
      continue := true;
    end Write;
    procedure Enum is new Lexicographic_Enumeration(Write);

  begin
    Enum(1,n,d,acc,cont);
    put("Counted "); put(cnt,1); put_line(" solutions.");
  end Enumerate_and_Solve;

  procedure Create_and_Solve_Start_System ( n : in natural32 ) is

    d,s : Standard_Natural_Vectors.Vector(1..integer32(n));
    c : constant Standard_Complex_Vectors.Vector(1..integer32(n))
      := Random_Vector(1,integer32(n));
    q : Poly_Sys(1..integer32(n));
    td,k : natural32 := 0;
    cp : Standard_Natural_Vectors.Vector(1..integer32(n)-1);
    ans : character;

  begin
    for i in d'range loop
      put("Give degree of equation "); put(i,1);
      d(i) := 0; put(" : "); get(d(i));
    end loop;
    q := Start_System(d,c);
    put_line("The created start system : "); put(q);
    td := Product(d);
    put("The start system has "); put(td,1); put_line(" solutions.");
    cp := Consecutive_Products(d);
    put("Consecutive products : "); put(cp); new_line;
    Enumerate_and_Solve(d,c);
    loop
      put("Give index of root : "); get(k);
      s := Root_Map(n,k,d,cp);
      put("The root map : "); put(s); new_line;
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Create_and_Solve_Start_System;

  procedure Create_and_Solve_Start_System ( p : in Poly_Sys ) is

    c : constant Standard_Complex_Vectors.Vector(p'range)
      := Random_Vector(1,p'last);
    d : constant Standard_Natural_Vectors.Vector(p'range) := Degrees(p);
    td : constant natural32 := Product(d);
    q : constant Poly_Sys(p'range) := Start_System(p,c);
    qsols : Solution_List;
    sum : double_float;

  begin
    put_line("A start system based on the total degree :"); put(q);
    put("The total degree equals "); put(td,1); new_line;
    qsols := Solve(d,c);
    put_line("The residuals : ");
    Write_residuals(Standard_Output,q,qsols,sum);
    put("The sum of all residuals : "); put(sum); new_line;
  end Create_and_Solve_Start_System;

  procedure Main is

    lp : Link_to_Poly_Sys;
    n : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("MENU to create and solve a total-degree based start system.");
    put_line("  1. creation based on given polynomial system; or");
    put_line("  2. user supplies only vector of degrees.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    new_line;
    if ans = '1' then
      get(lp); 
      Create_and_Solve_Start_System(lp.all);
    else
      put("Give the number of equations : "); get(n);
      Create_and_Solve_Start_System(n);
    end if;
  end Main;

end Test_Start_Systems;
