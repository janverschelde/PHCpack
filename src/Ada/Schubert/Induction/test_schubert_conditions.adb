with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
-- with Standard_Random_Numbers;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;       use Multprec_Natural_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs_io;       use Standard_Natural_VecVecs_io;
with Standard_Natural_Matrices_io;      use Standard_Natural_Matrices_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with DoblDobl_Random_Vectors;           use DoblDobl_Random_Vectors;
with DoblDobl_Random_Matrices;          use DoblDobl_Random_Matrices;
with QuadDobl_Random_Vectors;           use QuadDobl_Random_Vectors;
with QuadDobl_Random_Matrices;          use QuadDobl_Random_Matrices;
with Symbol_Table,Symbol_Table_io;      use Symbol_Table;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with Matrix_Indeterminates;             use Matrix_Indeterminates;
with Standard_Complex_Poly_Matrices_io; use Standard_Complex_Poly_Matrices_io;
with DoblDobl_Complex_Poly_Matrices_io; use DoblDobl_Complex_Poly_Matrices_io;
with QuadDobl_Complex_Poly_Matrices_io; use QuadDobl_Complex_Poly_Matrices_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Embed_Polynomials;        use Standard_Embed_Polynomials;
with Standard_Numerical_Rank;
with Brackets_io;                       use Brackets_io;
with Bracket_Monomials_io;              use Bracket_Monomials_io;
with Standard_Bracket_Polynomials_io;   use Standard_Bracket_Polynomials_io;
with Standard_Bracket_Systems;
with Standard_Bracket_Systems_io;       use Standard_Bracket_Systems_io;
with DoblDobl_Bracket_Polynomials;
-- with QuadDobl_Bracket_Polynomials;
with Bracket_Polynomial_Convertors;     use Bracket_Polynomial_Convertors;
-- with Plane_Representations;
with Standard_Matrix_Inversion;
with Symbolic_Schubert_Conditions;      use Symbolic_Schubert_Conditions;
with Checker_Boards,Checker_Moves;      use Checker_Boards,Checker_Moves;
with Checker_Boards_io;                 use Checker_Boards_io;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;
with Checker_Posets_io;
with Intersection_Posets;               use Intersection_Posets;
with Intersection_Posets_io;            use Intersection_Posets_io;
with Checker_Homotopies;
with Numeric_Schubert_Conditions;       use Numeric_Schubert_Conditions;
with Remember_Symbolic_Minors;          use Remember_Symbolic_Minors;
with Setup_Flag_Homotopies;             use Setup_Flag_Homotopies;
with Start_Flag_Homotopies;             use Start_Flag_Homotopies;
with Moving_Flag_Homotopies;            use Moving_Flag_Homotopies;
with Schubert_Posets;
with Flag_Transformations;
with Main_Schubert_Induction;

package body Test_Schubert_Conditions is

  procedure Symbolic_Plane ( n,k : in integer32 ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := General_Localization_Map(n,k);
    dim : constant natural32 := Dimension(locmap);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);

  begin
    put("A general localization pattern for a "); put(k,1);
    put("-plane in "); put(n,1); put_line("-space :");
    put(locmap);
    put("The dimension of the space of "); put(k,1);
    put("-planes is "); put(dim,1); new_line;
    Initialize_Symbols(dim,locmap);
    xpm := Symbolic_Form_of_Plane(n,k,locmap);
    put("The symbolic "); put(k,1); put("-plane in ");
    put(n,1); put_line("-space :"); put(xpm);
  end Symbolic_Plane;

  function Identity ( n,nv : integer32 )
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    one : Term;

  begin
    one.dg := new Standard_Natural_Vectors.Vector'(1..nv => 0);
    one.cf := Create(1.0);
    for i in 1..n loop
      for j in 1..n loop
        if i = j
         then res(i,j) := Create(one);
         else res(i,j) := Null_Poly;
        end if;
      end loop;
    end loop;
    Clear(one);
    return res;
  end Identity;

  procedure Elaborate_Flag_Minors
               ( n,k,f,i : in natural32; fm : in Bracket_Polynomial;
                 A : in Standard_Complex_Matrices.Matrix ) is

    use Standard_Bracket_Systems;

    nq : constant natural32 := Number_of_Equations(n,k,f,i);
    fms : Bracket_System(1..integer32(nq));
    ind : integer32 := 0;
    p : Bracket_Polynomial;

  begin
    Flag_Minor_Polynomials(fm,fms,ind);
    put_line("The flag minors :"); put(fms);
    put_line("The elaborated equations : ");
    for j in fms'range loop
      p := Elaborate_One_Flag_Minor
             (integer32(n),integer32(k),integer32(f),integer32(i),fms(j),A);
      put_line(p);
    end loop;
  end Elaborate_Flag_Minors;

  procedure Impose_Schubert_Condition ( n,k : in natural32 ) is

    f,i,m,r : natural32 := 0;
    fm : Bracket_Polynomial;

  begin
    put("Give the dimension of the meeting plane : "); get(f);
    put("Give the dimension of the intersection  : "); get(i);
    m := k + f;  -- number of columns in the matrix [ X | F ]
    r := m - i;  -- rank of the minors
    put("  X meets F("); put(f,1); put(") in a "); put(i,1);
    put("-plane : Rank([ X | F("); put(f,1); put(") ]) = ");
    put(r,1); new_line;
    put("  => all "); put(r+1,1); put("-by-"); put(r+1,1);
    put(" minors of a "); put(n,1); put("-by-"); put(m,1);
    put_line(" matrix must be zero");
    put("this leads to "); put(Number_of_Equations(n,k,f,i),1);
    put_line(" minor equations ...");
    if ((r+1 > n) or (r+1 > m)) then
      put_line("  trivial condition, no minor equations");
    else
      fm := Flag_Minors(n,k,f,i);
      put_line("The minor equations : "); put(fm); new_line;
      declare
        A : constant Standard_Complex_Matrices.Matrix
                       (1..integer32(n),1..integer32(f))
          := Random_Matrix(n,f);
      begin
        Elaborate_Flag_Minors(n,k,f,i,fm,A);
      end;
    end if;
  end Impose_Schubert_Condition;

  procedure Standard_Test_Minor_Expansions ( n,k : in integer32 ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := General_Localization_Map(n,k);
    dim : constant natural32 := Dimension(locmap);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nq : constant natural32 
       := Remember_Symbolic_Minors.Number_of_Minors(natural32(n),natural32(k));
    mt : constant Standard_Symbolic_Minors(integer32(nq))
       := Create(natural32(n),natural32(k),xpm);
   
  begin
    put("A general localization pattern for a "); put(k,1);
    put("-plane in "); put(n,1); put_line("-space :");
    put(locmap);
    put("The dimension of the space of "); put(k,1);
    put("-planes is "); put(dim,1); new_line;
    Initialize_Symbols(dim,locmap);
    put("The symbolic "); put(k,1); put("-plane in ");
    put(n,1); put_line("-space :"); put(xpm);
    put("The number of maximal minors : "); put(nq,1); new_line;
    put_line("The remember table of symbolic minors : "); Write(mt);
    Query(mt,k);
    Impose_Schubert_Condition(natural32(n),natural32(k));
  end Standard_Test_Minor_Expansions;

  procedure DoblDobl_Test_Minor_Expansions ( n,k : in integer32 ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := General_Localization_Map(n,k);
    dim : constant natural32 := Dimension(locmap);
    xpm : constant DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nq : constant natural32 
       := Remember_Symbolic_Minors.Number_of_Minors(natural32(n),natural32(k));
    mt : constant DoblDobl_Symbolic_Minors(integer32(nq))
       := Create(natural32(n),natural32(k),xpm);
   
  begin
    put("A general localization pattern for a "); put(k,1);
    put("-plane in "); put(n,1); put_line("-space :");
    put(locmap);
    put("The dimension of the space of "); put(k,1);
    put("-planes is "); put(dim,1); new_line;
    Initialize_Symbols(dim,locmap);
    put("The symbolic "); put(k,1); put("-plane in ");
    put(n,1); put_line("-space :"); put(xpm);
    put("The number of maximal minors : "); put(nq,1); new_line;
    put_line("The remember table of symbolic minors : "); Write(mt);
    Query(mt,k);
    Impose_Schubert_Condition(natural32(n),natural32(k));
  end DoblDobl_Test_Minor_Expansions;

  procedure QuadDobl_Test_Minor_Expansions ( n,k : in integer32 ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := General_Localization_Map(n,k);
    dim : constant natural32 := Dimension(locmap);
    xpm : constant QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nq : constant natural32 
       := Remember_Symbolic_Minors.Number_of_Minors(natural32(n),natural32(k));
    mt : constant QuadDobl_Symbolic_Minors(integer32(nq))
       := Create(natural32(n),natural32(k),xpm);
   
  begin
    put("A general localization pattern for a "); put(k,1);
    put("-plane in "); put(n,1); put_line("-space :");
    put(locmap);
    put("The dimension of the space of "); put(k,1);
    put("-planes is "); put(dim,1); new_line;
    Initialize_Symbols(dim,locmap);
    put("The symbolic "); put(k,1); put("-plane in ");
    put(n,1); put_line("-space :"); put(xpm);
    put("The number of maximal minors : "); put(nq,1); new_line;
    put_line("The remember table of symbolic minors : "); Write(mt);
    Query(mt,k);
    Impose_Schubert_Condition(natural32(n),natural32(k));
  end QuadDobl_Test_Minor_Expansions;

  procedure Test_Minor_Expansions ( n,k : integer32 ) is

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision of the coefficients : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test_Minor_Expansions(n,k);
      when '1' => DoblDobl_Test_Minor_Expansions(n,k);
      when '2' => QuadDobl_Test_Minor_Expansions(n,k);
      when others => null;
    end case;
  end Test_Minor_Expansions;

  function Read_Localization_Map ( n,k : integer32 )
             return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(1..n,1..k);
    ans : character;

  begin
    put_line("The default is a generalization localization map.");
    put("-> give your own localization map ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give a "); put(n,1); put("-by-"); put(k,1);
      put_line(" natural matrix : "); get(res);
    else
      res := General_Localization_Map(n,k);
    end if;
    return res;
  end Read_Localization_Map;

  procedure Standard_Flag_Conditions ( n,k : in integer32; b : in Bracket ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Bracket_Systems;

    fm : Bracket_System(b'range);
    flag : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
         := Random_Flag(n); -- := Random_Matrix(n,n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Read_Localization_Map(n,k);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nv : constant natural32 := Dimension(locmap);
    nq : natural32 := Number_of_Equations(natural32(n),b);
    ind : integer32 := 0;
    sys : Poly_Sys(1..integer32(nq));
    ans : character;
    cnteffeqs : natural32;

  begin
    put_line("The symbolic form of the plane : "); put(xpm);
    Explain_Equations(natural32(n),b,nq);
    put("We have "); put(nq,1); put(" equations in ");
    put(nv,1); put_line(" variables.");
    cnteffeqs := Number_of_NotAbove(natural32(n),b);
    put("Number of equations in efficient representation : ");
    put(cnteffeqs,1); put_line(".");
    put("Continue to compute the flag minors ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      fm := Flag_Minors(natural32(n),b);
      put_line("The flag minors as bracket polynomials :"); put(fm);
      declare
        fms : Bracket_System(1..integer32(nq));
      begin
        put("Continue to see the flag minors as row/column pairs ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          fms := Flag_Minor_System(nq,fm);
          put_line("Flag minors as pairs of rows and colums :"); put(fms);
        end if;
      end;
     -- put("The linear system in "); put(k,1); put_line("-brackets :");
      for i in b'range loop
        if fm(i) /= Null_Bracket_Poly then
          nq := Number_of_Equations
                  (natural32(n),natural32(k),b(i),natural32(i));
          declare
            fms : constant Bracket_system(1..integer32(nq))
                := Flag_Minor_System(nq,fm(i));
           -- linsys : Bracket_System(1..nq);
           -- polsys : Poly_Sys(1..nq);
          begin
            for j in fms'range loop
             -- linsys(j) := Elaborate_One_Flag_Minor(n,k,b(i),i,fms(j),flag);
           -- polsys(j) := Elaborate_One_Flag_Minor(n,k,b(i),i,fms(j),xpm,flag);
              ind := ind + 1;
              sys(ind) := Elaborate_One_Flag_Minor
                            (n,k,integer32(b(i)),i,fms(j),xpm,flag);
            end loop;
           -- put_line(linsys);
           -- put_line(polsys);
          end;
        end if;
      end loop;
      put("Continue to see the polynomial equations ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put_line("The polynomial system :"); put_line(sys);
      end if;
    end if;
  end Standard_Flag_Conditions;

  procedure DoblDobl_Flag_Conditions ( n,k : in integer32; b : in Bracket ) is

    use DoblDobl_Complex_Poly_Systems;
    use Standard_Bracket_Systems;

    fm : Bracket_System(b'range);
    flag : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
         := Random_Flag(n); -- := Random_Matrix(n,n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Read_Localization_Map(n,k);
    xpm : constant DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nv : constant natural32 := Dimension(locmap);
    nq : natural32 := Number_of_Equations(natural32(n),b);
    ind : integer32 := 0;
    sys : Poly_Sys(1..integer32(nq));
    ans : character;
    cnteffeqs : natural32;

  begin
    put_line("The symbolic form of the plane : "); put(xpm);
    Explain_Equations(natural32(n),b,nq);
    put("We have "); put(nq,1); put(" equations in ");
    put(nv,1); put_line(" variables.");
    cnteffeqs := Number_of_NotAbove(natural32(n),b);
    put("Number of equations in efficient representation : ");
    put(cnteffeqs,1); put_line(".");
    put("Continue to compute the flag minors ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      fm := Flag_Minors(natural32(n),b);
      put_line("The flag minors as bracket polynomials :"); put(fm);
      declare
        fms : Bracket_System(1..integer32(nq));
      begin
        put("Continue to see the flag minors as row/column pairs ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          fms := Flag_Minor_System(nq,fm);
          put_line("Flag minors as pairs of rows and colums :"); put(fms);
        end if;
      end;
     -- put("The linear system in "); put(k,1); put_line("-brackets :");
      for i in b'range loop
        if fm(i) /= Null_Bracket_Poly then
          nq := Number_of_Equations
                  (natural32(n),natural32(k),b(i),natural32(i));
          declare
            fms : constant Bracket_system(1..integer32(nq))
                := Flag_Minor_System(nq,fm(i));
           -- linsys : Bracket_System(1..nq);
           -- polsys : Poly_Sys(1..nq);
            dd_bp : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
          begin
            for j in fms'range loop
             -- linsys(j) := Elaborate_One_Flag_Minor(n,k,b(i),i,fms(j),flag);
           -- polsys(j) := Elaborate_One_Flag_Minor(n,k,b(i),i,fms(j),xpm,flag);
              ind := ind + 1;
              dd_bp := Convert(fms(j));
              sys(ind) := Elaborate_One_Flag_Minor
                            (n,k,integer32(b(i)),i,dd_bp,xpm,flag);
              DoblDobl_Bracket_Polynomials.Clear(dd_bp);
            end loop;
           -- put_line(linsys);
           -- put_line(polsys);
          end;
        end if;
      end loop;
      put("Continue to see the polynomial equations ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put_line("The polynomial system :"); put_line(sys);
      end if;
    end if;
  end DoblDobl_Flag_Conditions;

  procedure Flag_Conditions ( n,k : in integer32; bm : in Bracket_Monomial ) is

    nb : constant integer32 := integer32(Number_of_Brackets(bm));
    cd : constant Array_of_Brackets(1..nb) := Create(bm);
    b : Bracket(1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Read_Localization_Map(n,k);
    nv : constant natural32 := Dimension(locmap);
    nq : natural32;
    ans : character;
    sum : natural32 := 0;

  begin
    for i in 1..nb loop
      b := cd(i).all;
      Explain_Equations(natural32(n),b,nq);
      put("Bracket "); Brackets_io.put(b);
      put(" imposes "); put(nq,1); put(" equations on ");
      put(nv,1); put_line(" variables.");
      put("-> continue with next bracket ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      sum := sum + nq;
    end loop;
    put(bm); put(" imposes "); put(sum,1);
    put(" equations on "); put(nv,1); put_line(" variables.");
  end Flag_Conditions;

  procedure Test_Schubert_Conditions ( n,k : in integer32 ) is

    lambda : Bracket(1..k);
    m : integer32;
    ans : character;

  begin
    put("Give a bracket : "); get(lambda);
    put("Our "); put(k,1); put("-plane X in "); 
    put(n,1); put("-space is subject to "); put(lambda); put_line(" : ");
    for i in 1..k loop
      m := k + integer32(lambda(i)) - i;
      put("  X meets F("); put(lambda(i),1); put(") in a "); put(i,1); 
      put("-plane : Rank([ X | F("); put(lambda(i),1); put(") ]) = ");
      put(m,1); new_line;
    end loop;
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Flag_Conditions(n,k,lambda);
      when '1' => DoblDobl_Flag_Conditions(n,k,lambda);
      when others => null;
    end case;
  end Test_Schubert_Conditions;

  function Generate_Point
             ( n,k : integer32; b : Bracket;
               flag : Standard_Complex_Matrices.Matrix ) 
             return Standard_Complex_Matrices.Matrix is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Matrices.Matrix(1..n,1..k);

  begin
    for j in b'range loop       -- the j-th column is a random combination
      declare                   -- of the first b(j) vectors of the flag
        rv : Standard_Complex_Vectors.Vector(1..integer32(b(j)));
      begin
        rv := Random_Vector(1,integer32(b(j)));
        for i in 1..n loop
          res(i,j) := Create(0.0);
          for kk in 1..integer32(b(j)) loop
            res(i,j) := res(i,j) + rv(kk)*flag(i,kk);
          end loop;
        end loop;
      end;
    end loop;
    return res;
  end Generate_Point;

  function Generate_Point
             ( n,k : integer32; b : Bracket;
               flag : DoblDobl_Complex_Matrices.Matrix ) 
             return DoblDobl_Complex_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);

  begin
    for j in b'range loop       -- the j-th column is a random combination
      declare                   -- of the first b(j) vectors of the flag
        rv : DoblDobl_Complex_Vectors.Vector(1..integer32(b(j)));
      begin
        rv := Random_Vector(1,integer32(b(j)));
        for i in 1..n loop
          res(i,j) := Create(integer(0));
          for kk in 1..integer32(b(j)) loop
            res(i,j) := res(i,j) + rv(kk)*flag(i,kk);
          end loop;
        end loop;
      end;
    end loop;
    return res;
  end Generate_Point;

  function Generate_Point
             ( n,k : integer32; b : Bracket;
               flag : QuadDobl_Complex_Matrices.Matrix ) 
             return QuadDobl_Complex_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);

  begin
    for j in b'range loop       -- the j-th column is a random combination
      declare                   -- of the first b(j) vectors of the flag
        rv : QuadDobl_Complex_Vectors.Vector(1..integer32(b(j)));
      begin
        rv := Random_Vector(1,integer32(b(j)));
        for i in 1..n loop
          res(i,j) := Create(integer(0));
          for kk in 1..integer32(b(j)) loop
            res(i,j) := res(i,j) + rv(kk)*flag(i,kk);
          end loop;
        end loop;
      end;
    end loop;
    return res;
  end Generate_Point;

  procedure Eliminate ( x : in out Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;

    ind : integer32;
    m : Complex_Number;

  begin
    for k in x'range(2) loop           -- first eliminate top
      for i in k+1..x'last(1) loop
        x(i,k) := x(i,k)/x(k,k);
      end loop;
      x(k,k) := Create(1.0);
      for j in k+1..x'last(2) loop     -- eliminate x(k,j)
        for i in k+1..x'last(1) loop
          x(i,j) := x(i,j) - x(k,j)*x(i,k);
        end loop;
        x(k,j) := Create(0.0);
      end loop;
    end loop;
    ind := x'last(1);
    for k in reverse 2..x'last(2) loop -- then eliminate bottom
      for j in 1..(k-1) loop
        m := x(ind,j)/x(ind,k);
        for i in x'range(1) loop
          x(i,j) := x(i,j) - m*x(i,k);
        end loop;
      end loop;
      ind := ind - 1;
    end loop;
  end Eliminate;

  procedure Divide_Pivots ( x : in out Standard_Complex_Matrices.Matrix;
                            b : in Bracket ) is

    use Standard_Complex_Numbers;

    pivot : Complex_Number;

  begin
    for j in b'range loop
      pivot := x(integer32(b(j)),j);
      for i in 1..integer32(b(j)) loop
        x(i,j) := x(i,j)/pivot;
      end loop;
    end loop;
  end Divide_Pivots;

  procedure Eliminate_Pivots
              ( x : in out Standard_Complex_Matrices.Matrix;
                b : in Bracket ) is

    use Standard_Complex_Numbers;

    fac : Complex_Number;

  begin
    for j in b'range loop
      -- pivot is at x(b(j),j), make zeroes at the right of pivot
      for k in j+1..x'last(2) loop
        fac := x(integer32(b(j)),k);
        for i in 1..integer32(b(j)) loop
          x(i,k) := x(i,k) - fac*x(i,j);
        end loop;
       -- x(integer32(b(j)),k) := Create(0.0);
      end loop;
    end loop;
  end Eliminate_Pivots;

  function Generate_Standard_Point
             ( n,k : integer32; b : Bracket;
               flag : Standard_Complex_Matrices.Matrix ) 
             return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..k)
        := Generate_Point(n,k,b,flag);

  begin
    Eliminate(res);
    return res;
  end Generate_Standard_Point;

  function Full_Localization_Map
             ( n,k : in integer32 ) return Standard_Natural_Matrices.Matrix is 

    res : Standard_Natural_Matrices.Matrix(1..n,1..k);

  begin
    for i in 1..n loop
      for j in 1..k loop
        res(i,j) := 2;
      end loop;
    end loop;
    return res;
  end Full_Localization_Map;

  function Stiefel_Localization_Map
             ( n,k : in integer32; b : in Bracket )
             return Standard_Natural_Matrices.Matrix is 

    res : Standard_Natural_Matrices.Matrix(1..n,1..k);

  begin
    for i in 1..n loop
      for j in 1..k loop
        res(i,j) := 2;
      end loop;
    end loop;
    for i in b'range loop
      for j in 1..k loop
        res(integer32(b(i)),j) := 0;
      end loop;
      res(integer32(b(i)),i) := 1;
      for j in integer32(b(i))+1..n loop
        res(j,i) := 0;
      end loop;
    end loop;
    return res;
  end Stiefel_Localization_Map;

  function Full_Flatten
             ( x : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector is

    n : constant integer32 := x'last(1)*x'last(2);
    res : Standard_Complex_Vectors.Vector(1..n);
    ind : integer32 := 0;

  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        ind := ind + 1;
        res(ind) := x(i,j);
      end loop;
    end loop;
    return res;
  end Full_Flatten;

  function Full_Flatten
             ( x : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Vectors.Vector is

    n : constant integer32 := x'last(1)*x'last(2);
    res : DoblDobl_Complex_Vectors.Vector(1..n);
    ind : integer32 := 0;

  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        ind := ind + 1;
        res(ind) := x(i,j);
      end loop;
    end loop;
    return res;
  end Full_Flatten;

  function Full_Flatten
             ( x : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Vectors.Vector is

    n : constant integer32 := x'last(1)*x'last(2);
    res : QuadDobl_Complex_Vectors.Vector(1..n);
    ind : integer32 := 0;

  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        ind := ind + 1;
        res(ind) := x(i,j);
      end loop;
    end loop;
    return res;
  end Full_Flatten;

  function General_Flatten
             ( x : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector is

    n : constant integer32 := (x'last(1) - x'last(2))*x'last(2);
    res : Standard_Complex_Vectors.Vector(1..n);
    ind : integer32 := 0;

  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        if i > j and i <= x'last(1) - x'last(2) + j
         then ind := ind + 1;
              res(ind) := x(i,j);
        end if;
      end loop;
    end loop;
    return res;
  end General_Flatten;

  procedure Standard_Point_Test_at_Conditions ( n,k : in integer32 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;

    lambda : Bracket(1..k);
    nq : integer32;
    flag : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
         := Random_Matrix(natural32(n),natural32(n));
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Full_Localization_Map(n,k); -- := General_Localization_Map(n,k);
    dim : constant natural32 := Dimension(locmap);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nv : constant natural32 := Dimension(locmap);
    ans : character;

  begin
    new_line;
    put("Give a bracket : "); get(lambda);
    Explain_Equations(natural32(n),lambda,natural32(nq));
    put("The number of variables equals "); put(nv,1); put_line(".");
    declare
      p : constant Poly_Sys(1..nq) := Expand(n,k,nq,lambda,xpm,flag);
      x : constant Standard_Complex_Matrices.Matrix(1..n,1..k)
        := Generate_Point(n,k,lambda,flag);
      v : constant Standard_Complex_Vectors.Vector(1..n*k) := Full_Flatten(x);
      y : constant Standard_Complex_Vectors.Vector(1..nq) := Eval(p,v);
    begin
      put("Do you want to see the polynomial system (y/n) ? ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Initialize_Symbols(dim,locmap);
        put_line(p);
      end if;
      put_line("The residual at a random plane : "); put_line(y);
    end;
  end Standard_Point_Test_at_Conditions;

  procedure DoblDobl_Point_Test_at_Conditions ( n,k : in integer32 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;

    lambda : Bracket(1..k);
    nq : integer32;
    flag : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
         := Random_Matrix(natural32(n),natural32(n));
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Full_Localization_Map(n,k); -- := General_Localization_Map(n,k);
    dim : constant natural32 := Dimension(locmap);
    xpm : constant DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nv : constant natural32 := Dimension(locmap);
    ans : character;

  begin
    new_line;
    put("Give a bracket : "); get(lambda);
    Explain_Equations(natural32(n),lambda,natural32(nq));
    put("The number of variables equals "); put(nv,1); put_line(".");
    declare
      p : constant Poly_Sys(1..nq) := Expand(n,k,nq,lambda,xpm,flag);
      x : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..k)
        := Generate_Point(n,k,lambda,flag);
      v : constant DoblDobl_Complex_Vectors.Vector(1..n*k) := Full_Flatten(x);
      y : constant DoblDobl_Complex_Vectors.Vector(1..nq) := Eval(p,v);
    begin
      put("Do you want to see the polynomial system (y/n) ? ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Initialize_Symbols(dim,locmap);
        put_line(p);
      end if;
      put_line("The residual at a random plane : "); put_line(y);
    end;
  end DoblDobl_Point_Test_at_Conditions;

  procedure QuadDobl_Point_Test_at_Conditions ( n,k : in integer32 ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;

    lambda : Bracket(1..k);
    nq : integer32;
    flag : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
         := Random_Matrix(natural32(n),natural32(n));
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Full_Localization_Map(n,k); -- := General_Localization_Map(n,k);
    dim : constant natural32 := Dimension(locmap);
    xpm : constant QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    nv : constant natural32 := Dimension(locmap);
    ans : character;

  begin
    new_line;
    put("Give a bracket : "); get(lambda);
    Explain_Equations(natural32(n),lambda,natural32(nq));
    put("The number of variables equals "); put(nv,1); put_line(".");
    declare
      p : constant Poly_Sys(1..nq) := Expand(n,k,nq,lambda,xpm,flag);
      x : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..k)
        := Generate_Point(n,k,lambda,flag);
      v : constant QuadDobl_Complex_Vectors.Vector(1..n*k) := Full_Flatten(x);
      y : constant QuadDobl_Complex_Vectors.Vector(1..nq) := Eval(p,v);
    begin
      put("Do you want to see the polynomial system (y/n) ? ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Initialize_Symbols(dim,locmap);
        put_line(p);
      end if;
      put_line("The residual at a random plane : "); put_line(y);
    end;
  end QuadDobl_Point_Test_at_Conditions;

  procedure Point_Test_at_Conditions ( n,k : in integer32 ) is

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Point_Test_at_Conditions(n,k);
      when '1' => DoblDobl_Point_Test_at_Conditions(n,k);
      when '2' => QuadDobl_Point_Test_at_Conditions(n,k);
      when others => null;
    end case;
  end Point_Test_at_Conditions;

  procedure Truncate_Triangular_Part
              ( A : in out Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
  
  begin
    for j in A'range(2) loop
     -- for i in A'first(1)..j-1 loop
      A(j,j) := Create(1.0);
      for i in j+1..A'last(1) loop
        A(i,j) := Create(0.0);
      end loop;
    end loop;
  end Truncate_Triangular_Part;

 -- procedure Add_Random_Last_Columns
 --             ( A : in out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Adds a random multiple of the last column to the other columns
  --   to make the matrix no longer upper triangular.

 --   use Standard_Complex_Numbers;

 --   rnd : Complex_Number;

 -- begin
 --   for k in reverse A'first(2)+1..A'last(2) loop
      -- add multiple of k-th column to previous columns
 --     for j in A'first(2)..k-1 loop
 --       rnd := Standard_Random_Numbers.Random1;
 --       for i in A'range(1) loop
 --         A(i,j) := A(i,j) + rnd*A(i,k);
 --       end loop;
 --     end loop;
 --   end loop;
 -- end Add_Random_Last_Columns;

  procedure Point_Test_at_Minimal_Conditions ( n,k : in integer32 ) is

    lambda : Bracket(1..k);
    nq : integer32;
    flag : Standard_Complex_Matrices.Matrix(1..n,1..n)
         := Random_Matrix(natural32(n),natural32(n));
    trans : Standard_Complex_Matrices.Matrix(1..n,1..n)
          := Random_Matrix(natural32(n),natural32(n));
    locmap : Standard_Natural_Matrices.Matrix(1..n,1..k);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    x : Standard_Complex_Matrices.Matrix(1..n,1..k);
    nv : natural32;
    ans : character;

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;

  begin
    new_line;
    put("Give a bracket : "); get(lambda);
    locmap := Stiefel_Localization_Map(n,k,lambda);
    put_line("The localization map :"); put(locmap);
    nv := Dimension(locmap);
    Initialize_Symbols(nv,locmap);
    xpm := Symbolic_Form_of_Plane(n,k,locmap);
    put_line("The symbolic form of the solution plane :"); put(xpm);
    nq := integer32(Number_of_NotAbove(natural32(n),lambda));
    put("The number of equations equals "); put(nq,1); put_line(".");
    put("The number of variables equals "); put(nv,1); put_line(".");
    Truncate_Triangular_Part(flag);  -- to fix localization map
    Truncate_Triangular_Part(trans); -- upper triangular to multiply flag with
    x := Generate_Point(n,k,lambda,flag);
    Divide_Pivots(x,lambda);
    Eliminate_Pivots(x,lambda);
    put_line("The solution plane : "); put(x,3);
   -- Add_Random_Last_Columns(flag);   -- so no longer upper triangular
    flag := flag*trans;
   -- Explain_Equations(natural32(n),lambda,natural32(nq));
    declare
      p : constant Poly_Sys(1..nq)
        := Minimal_Expand(n,k,nq,lambda,xpm,flag);
       -- := Expand(n,k,nq,lambda,xpm,flag);
      v : constant Standard_Complex_Vectors.Vector(1..n*k) := Full_Flatten(x);
      y : constant Standard_Complex_Vectors.Vector(1..nq) := Eval(p,v);
    begin
      put("Do you want to see the polynomial system (y/n) ? ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put_line(p);
      end if;
      put_line("The residual at a random plane : "); put_line(y);
    end;
  end Point_Test_at_Minimal_Conditions;

  procedure Cheater_Homotopy ( n,k,nq : in integer32; b : in Bracket ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;

    flag : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
         := Random_Matrix(natural32(n),natural32(n));
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := General_Localization_Map(n,k);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
        := Symbolic_Form_of_Plane(n,k,locmap);
    p : constant Poly_Sys(1..nq) := Expand(n,k,nq,b,xpm,flag);
    x : constant Standard_Complex_Matrices.Matrix(1..n,1..k)
      := Generate_Standard_Point(n,k,b,flag);
    v : constant Standard_Complex_Vectors.Vector(1..(n-k)*k)
      := General_Flatten(x);
    y : constant Standard_Complex_Vectors.Vector(1..nq) := Eval(p,v);
    file : file_type;
    s : Solution(v'last);
    sols : Solution_List;

  begin
    put_line("After elimination : "); put(x,3);
    put_line("The residual at the random point : "); put_line(y);
    s.t := Create(0.0); s.m := 1; s.v := v;
    s.err := 0.0; s.rco := 1.0; s.res := 0.0;
    Add(sols,s);
    Initialize_Symbols(natural32(s.n),locmap);
    new_line;
    skip_line;
    put_line("Reading the name of a file to write the system on...");
    skip_line;
    Read_Name_and_Create_File(file);
    put_line(file,p);
    new_line(file);
    put_line(file,"THE SOLUTIONS ");
    put(file,1,natural32(s.n),sols);
    close(file);
    new_line;
    put_line("See the output file for system and its solution");
  end Cheater_Homotopy;

  procedure Test_Cheater_Homotopy ( n,k : in integer32 ) is

    lambda : Bracket(1..k);
    nq : integer32;

  begin
    new_line;
    put("Give a bracket : "); get(lambda);
    Explain_Equations(natural32(n),lambda,natural32(nq));
    Cheater_Homotopy(n,k,nq,lambda);
  end Test_Cheater_Homotopy;

  procedure Generalizing_Moving_Flags ( n,dim : in integer32 ) is

    use Standard_Complex_Poly_Matrices;

    m : constant integer32 := integer32(Number_of_Moves(natural32(n)));
    all_moves : constant Standard_Natural_VecVecs.VecVec(1..m)
              := Checker_Posets.Generalizing_Moves(n);
    ans : character;
    p,p1 : Standard_Natural_Vectors.Vector(1..n);
    mf,t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    sf,tmp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    acc : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n) := Identity(n,dim);
    f,a : integer32;
    ind : natural32;
    tsb : Symbol;

  begin
    put(all_moves);
    put("Do you wish to see all the boards ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      for i in all_moves'first..all_moves'last-1 loop
        p := all_moves(i).all; mf := Moving_Flag(p);
        f := Falling_Checker(p); a := Ascending_Checker(p,f);
        p1 := all_moves(i+1).all;
        t := Transformation(n,integer32(p1(f)));
        Write_Permutation(p,mf,t);
        put("After Move #"); put(i,1); 
        put(" (a,f) = ("); put(a,1); put(","); put(f,1);
        put(") : ");
        tsb := Matrix_Indeterminates.X_ij(natural32(n+1-a),natural32(f));
        ind := Symbol_Table.get(tsb);
        Symbol_Table_io.put(tsb); new_line;
        Standard_Complex_Poly_Matrices.Clear(sf);
        sf := Symbolic_Transformation(dim,integer32(ind),t);
        put_line("The symbolic form of the transformation : "); put(sf);
        tmp := acc*sf;
        Copy(tmp,acc); Clear(tmp);
        put_line("The transformed moving flag : "); put(acc);
      end loop;
      put_line("The final configuration : ");
      p := all_moves(all_moves'last).all; mf := Moving_Flag(p);
      Write_Permutation(p,mf);
    end if;
  end Generalizing_Moving_Flags;

  procedure Symbolic_Moving_Flags ( n : in integer32 ) is

    p : Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    mf : Standard_Natural_Matrices.Matrix(1..n,1..n) := Moving_Flag(p);
    dim : constant natural32 := Dimension(mf);
    t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    sf : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    ans : character;
    r,d : integer32;
    ind : natural32;
    tsb : Symbol;
    
  begin
    new_line;
    Write_Permutation(p,mf);
    Initialize_Symbols(dim,mf);
    sf := Symbolic_Form_of_Plane(n,n,mf);
    put_line("The symbolic form of the full flag :"); put(sf);
    put("Do you want to see all moves ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Generalizing_Moving_Flags(n,integer32(dim));
    end if;
    loop
      new_line;
      put("continue interactively ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      new_line;
      put_line("Reading a permutation ...");
      Read_Permutation(p);
      d := Descending_Checker(p);
      r := Rising_Checker(p,d);
      put("(r,d) = ("); put(r,1); put(","); put(d,1); put(") : ");
      tsb := Matrix_Indeterminates.X_ij(natural32(n+1-r),natural32(d));
      ind := Symbol_Table.get(tsb);
      Symbol_Table_io.put(tsb); new_line;
      t := Transformation(n,integer32(p(d)));
      put_line("The pattern for the transformation : "); put(t);
      Standard_Complex_Poly_Matrices.Clear(sf);
      sf := Symbolic_Transformation(integer32(dim),integer32(ind),t);
      mf := Moving_Flag(p);
      Write_Permutation(p,mf);
      put_line("The symbolic form of the transformation :"); put(sf);
    end loop;
  end Symbolic_Moving_Flags;

 -- procedure Compare_Patterns ( p,q : in Standard_Natural_Matrices.Matrix ) is
--
  -- DESCRIPTION :
  --   On input are two row patterns, p is the previous and q is the current.
--
 --   order : natural := 0;
 --   cnt : natural := 0;
--
 -- begin
 --   for i in p'range(1) loop
 --     for j in p'range(2) loop
 --        if p(i,j) = 2 then
 --          cnt := cnt + 1;
 --          if p(i,j) /= q(i,j) then
 --            order := 1;
 --            for k in i+1..p'last(1) loop
 --              if p(k,j) = 2
 --               then order := 2;
 --              end if;
 --            end loop;
 --          end if;
 --        end if;
 --        exit when (order > 1);
 --     end loop;
 --     exit when (order > 1);
 --   end loop;
 --   if order = 0 then
 --     put_line("patterns have identical free variables");
 --   elsif order = 1 then
 --     put("free variables shifted at "); put(cnt,1); new_line;
 --   else
 --     put("permuted free variables at "); put(cnt,1); new_line;
 --   end if;
 -- end Compare_Patterns;

  procedure Symbolic_Localization_Patterns ( n,k : in integer32 ) is

    rows,cols : Standard_Natural_Vectors.Vector(1..k);
    ip : constant Standard_Natural_Vectors.Vector(1..n)
       := Identity_Permutation(natural32(n));
    mf : constant Standard_Complex_Matrices.Matrix(1..n,1..n) := One_Flag(n);
    nf : Standard_Complex_Matrices.Matrix(1..n,1..n);
    ps : Poset;
    cnt : integer32 := 0;

    procedure Display_Transition_Equations
                ( a,b : in Standard_Complex_Poly_Matrices.Matrix ) is

    -- DESCRIPTION :
    --   Displays the transition equations between a and b.

    begin
      for i in a'range(1) loop
        for j in a'range(2) loop
         -- if Degree(a(i,j)) > 0 or Degree(b(i,j)) > 0 then
          if not Equal(a(i,j),b(i,j)) then
            put("transition condition at ("); put(i,1);
            put(","); put(j,1); put(") : ");
            put(a(i,j)); put(" = "); put(b(i,j)); new_line;
          end if;
        end loop;
      end loop;
    end Display_Transition_Equations;

    procedure Show_Transition
               -- ( pvt,
                ( pvx : in Standard_Complex_Poly_Matrices.Matrix;
                  invct,x : in Standard_Complex_Poly_Matrices.Matrix ) is

    -- DESCRIPTION :
    --   Shows the transition between localization patterns.

    -- ON ENTRY :
    --   pvt      n-by-n matrix represents symbolic transformation
    --            of the previous move;
    --   pvx      n-by-k matrix is localization pattern for k-plane
    --            of the previous move;
    --   invct    n-by-n matrix is inverse of the symbolic transformation
    --            of the current move;
    --   cx       n-by-k matrix is localization pattern for k-plane
    --            of the current move.

      use Standard_Complex_Numbers;
      use Standard_Complex_Poly_Matrices;

     -- valpvt : constant Matrix(pvt'range(1),pvt'range(2))
     --        := Evaluate_Transformation(pvt,Create(1.0));
      valict : constant Matrix(invct'range(1),invct'range(2))
             := Evaluate_Transformation(invct,Create(1.0));
     -- mvx : Matrix(pvx'range(1),pvx'range(2)) := valpvt*pvx;
     -- cmx : Matrix(pvx'range(1),pvx'range(2)) := valict*mvx;
      cmx : constant Matrix(pvx'range(1),pvx'range(2)) := valict*pvx;

    begin
      put_line("transition between two moves :");
     -- put_line("Value of previous transformation at t = 1 : ");
     -- put(valpvt);
     -- put_line("Form of solution plane after previous move : ");
     -- put(mvx);
      put_line("The current inverse transformation : ");
      put(invct);
      put_line("Value of current inverse transformation at t = 1 : ");
      put(valict);
      put_line("Solution after previous move times inverse transformation :");
      put(cmx);
      put_line("The current form of the solution plane : "); put(x);
      Display_Transition_Equations(cmx,x);
    end Show_Transition;

    procedure Show_Move
                ( nds : in Array_of_Nodes; i : in integer32;
                  p,q,rp,cp : in Standard_Natural_Vectors.Vector;
                  pvt,pvx : in out Standard_Complex_Poly_Matrices.Matrix ) is

    -- DESCRIPTION :
    --   Shows the status of the checker board at the i-th move.

    -- ON ENTRY :
    --   nds      array of nodes along a path;
    --   i        index of the current node;
    --   p        permutation at current i-th node;
    --   q        next permutation for the next node;
    --   rp       row positions of white checkers of i-th node;
    --   cp       column positions of white checkers of i-th node;
    --   pvt      n-by-n matrix of transformation of previous move;
    --   pvx      n-by-k matrix represents localization of previous move.

    -- ON RETURN :
    --   pvt      transformation of current move;
    --   pvx      symbolic form of solution plane in current move.

      use Standard_Complex_Numbers;
      use Standard_Complex_Matrices;
      use Standard_Complex_Poly_Matrices;

      f : constant integer32 := Checker_Moves.Falling_Checker(p);
      a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
      t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
        := Checker_Localization_Patterns.Transformation(n,integer32(q(f)));
      pf : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
         := Checker_Localization_Patterns.Moving_Flag(p);
      map : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
          := Checker_Localization_Patterns.Column_Pattern(n,k,p,rp,cp);
      dim : constant integer32
          := integer32(Checker_Localization_Patterns.Degree_of_Freedom(map));
      x,xt,xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
      st,invst,sf : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
      nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
      gamma : Complex_Number;

    begin
      put("Checker board at node "); put(i,1); put(" : "); 
      Checker_Posets_io.Write_Node(nds(i).all); new_line;
      Checker_Boards_io.Write_Permutation(standard_output,p,rp,cp,pf,t,map);
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Initialize_Homotopy_Symbols(natural32(dim),map);
      x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,map);
     -- put_line("Symbolic form of the plane : "); put(x);
      put("Order of symbols :");
      for j in 1..Symbol_Table.Number loop
        put(" "); Symbol_Table_io.put(Symbol_Table.get(j));
      end loop;
      new_line;
      xt := Standard_Embed_Polynomials.Add_Variables(x,1);
      put("falling checker : "); put(f,1);
      put(", ascending checker : "); put(a,1); 
      put(" -> gamma coordinates : ("); put(n+1-a,1);
      put(","); put(f,1); put_line(")");
      gamma := mf(n+1-a,f);
     -- put("gamma : "); put(gamma); new_line;
      nt := Numeric_Transformation(t,gamma);
      st := Symbolic_Transformation(dim+1,dim+1,gamma,t);
      invst := Inverse_Symbolic_Transformation(dim+1,dim+1,gamma,t);
     -- put_line("the symbolic transformation :"); put(st);
     -- put_line("the inverse symbolic transformation :"); put(invst);
     -- put_line("the product of transformation and its inverse :");
     -- put(st*invst);
      sf := Moving_Flag(nf,st);
     -- put_line("the moving flag : "); put(sf);
      xp := sf*xt;
      put_line("moving coordinates of the plane : "); put(xp);
      nf := nf*nt;
      if i > 1 then
       -- Show_Transition(pvt,pvx,invst,x);
        Show_Transition(pvx,invst,x);
      end if;
      Copy(st,pvt); Copy(x,pvx);
    end Show_Move;

    procedure Show_Patterns ( nds : in Array_of_Nodes; cont : out boolean ) is

      ans : character;
      p,q : Standard_Natural_Vectors.Vector(1..n);
      rp,cp,rq,cq : Standard_Natural_Vectors.Vector(1..k);
      stm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
      xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
      ind,homtp,ctr : integer32;
      stay_child : boolean;

    begin
      cnt := cnt + 1;
      put("Path "); put(cnt,1); put(" starts at ");
      Checker_Posets_io.Write_Node(nds(nds'first).all); new_line;
      put("Do you want to see all patterns along this path ? (y/n) ");
      Ask_Yes_or_No(ans);
      nf := Identity(n);
      if ans = 'y' then
        ind := ps.black'last; p := ps.black(ind).all;
        rp := nds(nds'first).rows; cp := nds(nds'first).cols;
        for i in nds'first..nds'last-1 loop
          q := ps.black(ind-1).all;
          rq := nds(i+1).rows; cq := nds(i+1).cols;
          stay_child := Checker_Posets.Is_Stay_Child(nds(i+1).all,nds(i).all);
          Checker_Homotopies.Define_Generalizing_Homotopy
            (standard_output,n,q,rq,cq,stay_child,homtp,ctr);
          Show_Move(nds,i,p,q,rp,cp,stm,xpm);
         -- put("continue display ? (y/n) ");
         -- Ask_Yes_or_No(ans);
         -- exit when (ans /= 'y');
          p := q; rp := rq; cp := cq; ind := ind - 1;
        end loop;
      end if;
      put("Do you want to see more paths ? (y/n) ");
      Ask_Yes_or_No(ans);
      cont := (ans = 'y');
    end Show_Patterns;
    procedure Show_Paths is new Enumerate_Paths_in_Poset(Show_Patterns);

  begin
    Main_Schubert_Induction.Read_Intersection_Conditions(ip,rows,cols);
    ps := Create(n,rows,cols);
    Show_Paths(ps);
  end Symbolic_Localization_Patterns;

  procedure Define_Moving_Flag_Homotopy ( n,k : in integer32 ) is

    use Standard_Complex_Poly_Matrices;

    p : Standard_Natural_Vectors.Vector(1..n);
    b : Board(1..n,1..n);
    mf,t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    xp : Standard_Natural_Matrices.Matrix(1..k,1..n);
    locmap : Standard_Natural_Matrices.Matrix(1..n,1..k);
    dim : natural32;
    d : integer32;
    xpm,xpt,tx : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    tpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    rows,cols : Standard_Natural_Vectors.Vector(1..k);
    happy : boolean;

  begin
    new_line;
    loop
      put_line("Reading a permutation ...");
      Read_Permutation(p);
      put_line("Reading rows and columns for white checkers ...");
      Read_Permutation(rows); Read_Permutation(cols);
      Check_Happiness(p,rows,cols,happy);
      exit when happy;
      put_line("The white checkers are not happy.  Please try again ...");
    end loop;
    b := Configuration(p);
    Place_White(b,rows,cols);
    mf := Moving_Flag(p);
    xp := Row_Pattern(n,k,p,rows,cols);
    Write(b,mf,p,rows,cols);
    put("The localization pattern for the ");
    put(k,1); put_line("-plane : "); put(xp);
    locmap := Column_Pattern(n,k,p,rows,cols);
    dim := Dimension(locmap);
    Initialize_Homotopy_Symbols(dim,locmap);
    xpm := Symbolic_Form_of_Plane(n,k,locmap);
    put("The symbolic "); put(k,1); put("-plane in ");
    put(n,1); put_line("-space :"); put(xpm);
    d := Descending_Checker(p);
    t := Transformation(n,integer32(p(d)));
    put_line("The pattern for the transformation : "); put(t);
    tpm := Symbolic_Transformation(integer32(dim)+1,integer32(dim)+1,t);
    put_line("The symbolic form of the transformation :"); put(tpm);
    xpt := Add_Variables(xpm,1);
    tx := tpm*xpt;
    put("The "); put(k,1); put_line("-plane in moving coordinates :");
    put(tx);
  end Define_Moving_Flag_Homotopy;

  procedure Test_One_Flag_Homotopy ( n,k : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;

    vf : constant Standard_Complex_Matrices.Matrix(1..n,1..n) := Random_Flag(n);
    mf : constant Standard_Complex_Matrices.Matrix(1..n,1..n) := Random_Flag(n);
    p,q : Standard_Natural_Vectors.Vector(1..n);
    rows,cols,cond : Standard_Natural_Vectors.Vector(1..k);
    happy : boolean;
    b : Board(1..n,1..n);
    f : Standard_Natural_Matrices.Matrix(1..n,1..n);

  begin
    new_line;
    put_line("Testing one flag homotopy ...");
    loop
      put_line("Reading a permutation ...");
      Read_Permutation(p);
      put_line("Reading rows and columns for white checkers ...");
      Read_Permutation(rows); Read_Permutation(cols);
      Check_Happiness(p,rows,cols,happy);
      exit when happy;
      put_line("The white checkers are not happy.  Please try again ...");
    end loop;
    b := Configuration(p);
    Place_White(b,rows,cols);
    f := Moving_Flag(p);
    Write(b,f,p,rows,cols);
    put_line("Reading the permutation at the parent ...");
    Read_Permutation(q);
    put_line("Reading an intersection condition ...");
    Read_Permutation(cond);
    declare
      c : constant Bracket(1..k) := Bracket(cond);
      locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
             := Column_Pattern(n,k,p,rows,cols);
      dim : constant integer32 := integer32(Dimension(locmap));
      nf : constant Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
      s : Link_to_Poly_Sys;
      x : Standard_Complex_Vectors.Vector(1..dim);
      xt : Standard_Complex_Vectors.Vector(1..dim+1);
      fail : boolean;
      res : double_float;
    begin
      Initialize_Homotopy_Symbols(natural32(dim),locmap);
      One_Flag_Homotopy(standard_output,n,k,q,p,rows,cols,c,vf,mf,nf,s);
      put_line("The one flag homotopy : "); put(s.all);
      Start_Solution(s.all,fail,x,res);
      if fail then
        put_line("No start solution found.");
      else
        put_line("The start solution is "); put_line(x);
        put("The residual : "); put(res,3); new_line;
        xt(x'range) := x; xt(xt'last) := Create(0.0);
        declare
          y : constant Standard_Complex_Vectors.Vector(s'range)
            := Eval(s.all,xt);
        begin
          put_line("The value of the start solution : "); put_line(y);
        end;
      end if;
    end;
  end Test_One_Flag_Homotopy;

  procedure Write_Coefficients ( ps : in Poset ) is
  begin
    put("Littlewood-Richardson coefficient at root : ");
    put(ps.white(ps.white'first).coeff); new_line;
    put("Littlewood-Richardson coefficients at leaves : ");
    declare
      nd : Link_to_Node := ps.white(ps.white'last);
    begin
      while nd /= null loop
        put("+"); put(nd.coeff);
        put("*("); put(nd.rows); put(",");
                   put(nd.cols); put(")");
        nd := nd.next_sibling;
      end loop;
      new_line;
    end;
  end Write_Coefficients;

  procedure Create_Schubert_Poset
               ( n,k,nb : in integer32; b : Bracket_Monomial ) is

    ips : constant Intersection_Poset(nb-1)
        := Schubert_Posets.Specialize(n,k,b);
    tmp : Poset_List;
    pnd : Link_to_Poset_Node;
    nbsols : Natural_Number;

  begin
    if ips.level < ips.m then
      put_line("Schubert problem has no solutions.");
    else
      put("Writing coefficients at root level ");
      put(ips.nodes'first,1); put_line(" :");
      tmp := ips.nodes(ips.nodes'first);
      while not Is_Null(tmp) loop
        pnd := Head_Of(tmp);
        Write_Coefficients(pnd.ps);
        tmp := Tail_Of(tmp);
      end loop;
      put("Writing coefficients at leaf level ");
      put(ips.level,1); put_line(" :");
      tmp := ips.nodes(ips.level);
      while not Is_Null(tmp) loop
        pnd := Head_Of(tmp);
        Write_Coefficients(pnd.ps);
        tmp := Tail_Of(tmp);
      end loop;
      nbsols := Schubert_Posets.Root_Count_at_Leaves(ips);
      put("Number of solutions : "); put(nbsols); new_line;
      Write_Expansion(ips);
      put_line("Counting roots from the leaves up ...");
      Schubert_Posets.Count_Roots(ips);
      declare
        ps : constant Poset := Head_Of(ips.nodes(ips.nodes'first)).ps;
        ns : constant Natural_Number := ps.white(ps.white'first).coeff;
      begin
        put("Coefficients at root of intersection poset : ");
        put(ns); 
        if Equal(ns,nbsols)
         then put_line("  okay");
         else put(" /= "); put(nbsols); put_line("!");
        end if;
      end;
      Write_Expansion(ips);
    end if;
  end Create_Schubert_Poset;

  procedure Define_Schubert_Systems ( n,k : in integer32 ) is

    b : Bracket_Monomial
      := Main_Schubert_Induction.Prompt_for_Bracket_Monomial;
    nb : constant natural32 := Number_of_Brackets(b);
    ans : character;

  begin
    put("-> "); put(nb,1);
    put(" intersection conditions : "); put(b); new_line;
    put("See interpretation of conditions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Flag_Conditions(n,k,b);
    end if;
    Create_Schubert_Poset(n,k,integer32(nb),b);
    Clear(b);
  end Define_Schubert_Systems;

  procedure Verify_a_Solution
              ( n,k : in integer32; cond : in Bracket;
                flag : in Standard_Complex_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                sol : in Solution; tol : in double_float ) is

    solplane : constant Standard_Complex_Matrices.Matrix(1..n,1..k)
             := Checker_Localization_Patterns.Map(locmap,sol.v);

  begin
    put_line("The solution as a matrix : ");
    put(solplane);
    for L in cond'range loop
      declare
        nbcols : constant integer32 := k + integer32(cond(L));
        A : Standard_Complex_Matrices.Matrix(1..n,1..nbcols);
        rnk,expected : natural32;
      begin
        for i in 1..n loop
          for j in 1..k loop
            A(i,j) := solplane(i,j);
          end loop;
          for j in 1..integer32(cond(L)) loop
            A(i,k+j) := flag(i,j);
          end loop;
        end loop;
        rnk := Standard_Numerical_Rank.Numerical_Rank(A,tol);
        put("Rank at condition "); put(L,1); put(" : ");
        put(rnk,1); 
        expected := natural32(k + integer32(cond(L)) - L);
        if rnk = expected then
          put(" = "); put(expected,1);
          put_line(", the expected rank.");
        else
          put(" /= "); put(expected,1);
          put_line(", the expected rank, bug!?");
        end if;
      end;
    end loop;
  end Verify_a_Solution;

  procedure Verify_Solutions
              ( n,k : in integer32; cond : in Bracket;
                flag : in Standard_Complex_Matrices.Matrix;
                sols : in Solution_List; tol : in double_float ) is

    p : constant Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    rows,cols : Standard_Natural_Vectors.Vector(1..k);
    locmap : Standard_Natural_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in cond'range loop
      rows(i) := cond(i);
      cols(i) := cond(i);
    end loop;
    locmap := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    put("The general localization map for a "); put(k,1);
    put("-plane in "); put(n,1); put_line("-space :");
    put(locmap);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Verify_a_Solution(n,k,cond,flag,locmap,ls.all,tol);
      tmp := Tail_Of(tmp);
    end loop;
  end Verify_Solutions;

  procedure Evaluate_Solutions
              ( n,k : in integer32; cond : in Bracket;
                flag : in Standard_Complex_Matrices.Matrix;
                sols : in Solution_List ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;

    sys : constant Poly_Sys := Expanded_Polynomial_Equations(n,k,cond,flag);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        y : constant Standard_Complex_Vectors.Vector := Eval(sys,ls.v);
      begin
        put_line("The value of the solution : "); put_line(y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Evaluate_Solutions;

  procedure Verify_Solutions_of_Schubert_Problem ( n,k : in integer32 ) is

    sols : Solution_List; 
    flagfile : file_type;
    flag : Standard_Complex_Matrices.Matrix(1..n,1..n);
    cond : Bracket(1..k);
    moved : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
          := Setup_Flag_Homotopies.Moved_Flag(n);
    inverse_moved : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
                  := Standard_Matrix_Inversion.Inverse(moved);
    idemat : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           :=  Setup_Flag_Homotopies.Identity(n);
    ranflag : Standard_Complex_Matrices.Matrix(1..n,1..n);
    rsd : double_float;
    tol : constant double_float := 1.0E-6;

    use Standard_Complex_Matrices;

  begin
    new_line;
    put_line("Reading the file name for the coordinates of the flag ...");
    Read_Name_and_Open_File(flagfile);
    new_line;
    put("Reading a "); put(n,1); put("-by-"); put(n,1);
    put_line(" matrix ...");
    get(flagfile,flag);
    close(flagfile);
    put_line("The matrix : "); put(flag);
    new_line;
    Read(sols);
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put("Give "); put(k,1); put(" increasing integers : ");
    get(cond);
    put("The bracket condition : "); put(cond); new_line;
    new_line;
    put_line("Verifying the bracket condition for the fixed flag ...");
    Verify_Solutions(n,k,cond,flag,sols,tol);
    new_line;
    put_line("The coordinates for the moved flag :"); put(moved);
    put_line("Verifying the bracket condition for the moved flag ...");
    Verify_Solutions(n,k,cond,inverse_moved,sols,tol);
    Flag_Transformations.Move_to_Generic_Flag(n,ranflag,rsd);
    put("The residual : "); put(rsd); new_line;
    new_line;
    put_line("Verifying the condition for the identity flag ...");
    Verify_Solutions(n,k,cond,idemat,sols,tol);
    new_line;
    put_line("Verifying the condition for the fixed flag ...");
    Verify_Solutions(n,k,cond,flag,sols,tol);
    new_line;
    put_line("Verifying the condition for the transformed random flag ...");
    Verify_Solutions(n,k,cond,ranflag,sols,tol);
    Evaluate_Solutions(n,k,cond,ranflag,sols);
  end Verify_Solutions_of_Schubert_Problem;

  procedure Main is

    n,k : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing symbolic & numeric Schubert intersection conditions...");
    new_line;
    put("Give the ambient dimension n : "); get(n);
    put("Give k, the dimension of the planes : "); get(k);
    new_line;
    put_line("MENU to test the operations on Schubert conditions :");
    put_line("  1. test minor expansions of one condition;");
    put_line("  2. test all conditions imposed by one bracket;");
    put_line("  3. generate polynomial system and test one point;");
    put_line("  4. set up a cheater homotopy for one point;");
    put_line("  5. symbolic forms of transformating the moving flag;");
    put_line("  6. symbolic forms of localization patterns for solution;");
    put_line("  7. define homotopy with a moving flag;");
    put_line("  8. test one flag homotopy for one move;");
    put_line("  9. test moving flag continuation;");
    put_line("  A. define all systems for a general Schubert problem;");
    put_line("  B. numerically verifying solutions of flag transformation.");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, 9, A, or B to make your choice : ");
    Ask_Alternative(ans,"123456789AB");
    case ans is
      when '1' => Test_Minor_Expansions(n,k);
      when '2' => Symbolic_Plane(n,k);
                  Test_Schubert_Conditions(n,k);
      when '3' =>
        new_line;
        put("Use the more efficient representation ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then Point_Test_at_Minimal_Conditions(n,k);
         else Point_Test_at_Conditions(n,k);
        end if;
      when '4' => Test_Cheater_Homotopy(n,k);
      when '5' => Symbolic_Moving_Flags(n);
      when '6' => Symbolic_Localization_Patterns(n,k);
      when '7' => Define_Moving_Flag_Homotopy(n,k);
      when '8' => Test_One_Flag_Homotopy(n,k);
      when '9' => Main_Schubert_Induction.Run_Moving_Flag_Continuation(n,k);
      when 'A' => Define_Schubert_Systems(n,k);
      when 'B' => Verify_Solutions_of_Schubert_Problem(n,k);
      when others => null;
    end case;
  end Main;

end Test_Schubert_Conditions;
