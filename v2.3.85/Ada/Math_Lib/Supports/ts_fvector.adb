with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;        use Standard_Integer_VecVecs_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Standard_Random_VecVecs;            use Standard_Random_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Face_Cardinalities;                 use Face_Cardinalities;

procedure ts_fvector is

-- DESCRIPTION :
--   Computes the f-vector of a polytope and checks the Euler-Poincare formula.

  procedure Write ( f : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the f-vector on screen.

    n : constant integer32 := f'last;

  begin
    put_line("The f-vector : ");
    put(" f(-1)     : "); put(f(-1)); new_line;
    put(" #vertices : "); put(f(0)); new_line;
    put(" #edges    : "); put(f(1)); new_line;
    for i in 2..(n-1) loop
      put(" #"); put(i,1); put("-faces  : "); put(f(i)); new_line;
    end loop;
    put(" f("); put(n,1); put(")      : "); put(f(n)); new_line;
  end Write;

  function Euler_Poincare ( f : Standard_Integer_Vectors.Vector )
                          return integer32 is

  -- DESCRIPTION :
  --   Computes the alternating sum: sum_{i in f'range} (-1)^(i) f(i).

    sum : integer32 := 0;
    pos : boolean := false;

  begin
    for i in f'range loop
      if pos
       then sum := sum + f(i);
       else sum := sum - f(i);
      end if;
      pos := not pos;
    end loop;
    return sum;
  end Euler_Poincare;

  function Euler_Poincare ( flb_pts : Array_of_Lists ) return integer32 is

  -- DESCRIPTION :
  --   Computes the alternating sum from the labels to the faces.

    sum,sign,len : integer32;

  begin
    sum := -1;
    sign := +1;
    put("- 1");
    for i in flb_pts'range loop
      sign := -sign;
      len := integer32(Length_Of(flb_pts(i)));
      if sign > 0 then
        put(" + ");
        sum := sum + len;
        put(len,1);
      else
        put(" - ");
        sum := sum - len;
        put(len,1);
      end if;
    end loop;
    if len > 0 then
      sign := -sign;
      if sign > 0
       then sum := sum+1;
       else sum := sum-1;
      end if;
    end if;
    put(" = "); put(sum,1);
    return sum;
  end Euler_Poincare;

  procedure Integer_Interactive_Testing is

    n,m,sum : integer32 := 0;

  begin
    put("Give the dimension n : "); get(n);
    put("Give the number of points that span the polytope : "); get(m);
    declare
      pts : Standard_Integer_VecVecs.VecVec(1..m);
      f : Standard_Integer_Vectors.Vector(-1..n);
      flb_pts : Array_of_Lists(0..n);
    begin
      put("Give "); put(m,1); put(" "); put(n,1);
      put_line("-dimensional integer vectors :");
      get(natural32(n),pts);
      put_line("Counting vertices, edges , .., facets, ...");
      f := fvector(pts); sum := Euler_Poincare(f);
      Write(f); 
      put("The result of the Euler-Poincare formula : "); put(sum,1);
      if sum /= 0 then
        put_line("   BUG DISCOVERED !!!");
        put_line("The labels to the faces :");
        flb_pts := Face_Labels(pts);
        put(flb_pts);
        sum := Euler_Poincare(flb_pts);
        if sum = 0
         then put_line("  Recomputation is OK...");
         else put_line("  BUG confirmed.");
        end if;
      else
        put_line("   OK.");
      end if;
    end;
  end Integer_Interactive_Testing;

  procedure Floating_Interactive_Testing is

    n,m,sum : integer32 := 0;
   -- tol : constant double_float := 10.0**(-8);

  begin
    put("Give the dimension n : "); get(n);
    put("Give the number of points that span the polytope : "); get(m);
    declare
      pts : Standard_Floating_VecVecs.VecVec(1..m);
      f : Standard_Integer_Vectors.Vector(-1..n);
      flb_pts : Array_of_Lists(0..n);
    begin
      put("Give "); put(m,1); put(" "); put(n,1);
      put_line("-dimensional floating point vectors :");
      get(natural32(n),pts);
      put_line("Counting vertices, edges , .., facets, ...");
      f := fvector(pts); sum := Euler_Poincare(f);
      Write(f);
      put("The result of the Euler-Poincare formula : "); put(sum,1);
      if sum /= 0 then
        put_line("   BUG DISCOVERED !!!");
        put_line("The labels to the faces : ");
        flb_pts := Face_Labels(pts);
        put(flb_pts);
        sum := Euler_Poincare(flb_pts);
        if sum = 0
         then put_line("  Recomputation is OK...");
         else put_line("  BUG confirmed.");
        end if;
      else
        put_line("   OK.");
      end if;
    end;
  end Floating_Interactive_Testing;

  function Floating_Random_Polytope 
             ( n,m : natural32 ) return Standard_Floating_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns m randomly chosen n-dimensional floating-point vectors.

    res : constant Standard_Floating_VecVecs.VecVec(1..integer32(m))
        := Random_VecVec(n,m);

  begin
    return res;
  end Floating_Random_Polytope;

  function Integer_Random_Polytope
             ( n,m : natural32; lower,upper : integer32 )
             return Standard_Integer_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns m randomly chosen n-dimensional vectors,
  --   with integer entries between lower and upper.

    res : Standard_Integer_VecVecs.VecVec(1..integer32(m));
    done : boolean;

    function Is_In ( i : integer32 ) return boolean is

    -- DESCRIPTION :
    --   Returns true if the ith vector already occurs in res(1..i-1).

      use Standard_Integer_Vectors;

    begin
      for j in 1..i-1 loop
        if Equal(res(j).all,res(i).all)
         then return true;
        end if;
      end loop;
      return false;
    end Is_In;

  begin
    for i in 1..integer32(m) loop
      res(i) := new Standard_Integer_Vectors.Vector(1..integer32(n));
      done := false;
      while not done loop
        for j in 1..integer32(n) loop
          res(i)(j) := Random(lower,upper);
        end loop;
        done := not Is_In(i);
      end loop;
    end loop;
    return res;
  end Integer_Random_Polytope;

  procedure Floating_Automatic_Testing is

    n,m,times,cnt,sum : integer32 := 0;
   -- tol : constant double_float := 10.0**(-8);
    bug : boolean := false;

  begin
    put("Give the dimension n : "); get(n);
    put("Give the number of points that span the polytope : "); get(m);
    put("Give the number of testing cycles : "); get(times);
    declare
      pts : Standard_Floating_VecVecs.VecVec(1..m);
      f : Standard_Integer_Vectors.Vector(-1..n);
      flb_pts : Array_of_Lists(0..n);
    begin
      for i in 1..times loop
        cnt := i;
        pts := Floating_Random_Polytope(natural32(n),natural32(m));
        f := fvector(pts); sum := Euler_Poincare(f);
        Write(f);
        put("The result of the Euler-Poincare formula : "); put(sum,1);
        if sum /= 0 then
          put_line("   BUG DISCOVERED !!!"); bug := true;
          put_line("The generated random configuration is");
          put(pts);
          put_line("The labels to the faces are :");
          flb_pts := Face_Labels(pts);
          put(flb_pts);
          sum := Euler_Poincare(flb_pts);
          if sum = 0
           then put_line("  Recomputation is OK...");
           else put_line("  BUG confirmed.");
          end if;
        else
          put_line("   OK.");  bug := false;
        end if;
        Standard_Floating_VecVecs.Clear(pts);
        exit when bug;
      end loop;
    end;
    if not bug then
      put("No bugs found, with "); put(times,1);
      put_line(" generated cases tested.");
      put("Dimension : "); put(n,1);
      put(" and #points : "); put(m,1); put_line(".");
    else
      put("Bug found at case "); put(cnt,1); put_line(".");
    end if;
  end Floating_Automatic_Testing;

  procedure Integer_Automatic_Testing is

    n,m,times,cnt,sum : integer32 := 0;
   -- tol : constant double_float := 10.0**(-8);
    bug : boolean := false;
    lower,upper : integer32 := 0;

  begin
    put("Give the dimension n : "); get(n);
    put("Give the number of points that span the polytope : "); get(m);
    put("Give lower bound on the entries : "); get(lower);
    put("Give upper bound on the entries : "); get(upper);
    put("Give the number of testing cycles : "); get(times);
    declare
      pts : Standard_Integer_VecVecs.VecVec(1..m);
      f : Standard_Integer_Vectors.Vector(-1..n);
      flb_pts : Array_of_Lists(0..n);
    begin
      for i in 1..times loop
        cnt := i;
        pts := Integer_Random_Polytope(natural32(n),natural32(m),lower,upper);
        f := fvector(pts); sum := Euler_Poincare(f);
        Write(f);
        put("The result of the Euler-Poincare formula : "); put(sum,1);
        if sum /= 0 then
          put_line("   BUG DISCOVERED !!!"); bug := true;
          put_line("The generated random configuration is");
          put(pts);
          put_line("The labels to the faces are ");
          flb_pts := Face_Labels(pts);
          put(flb_pts);
          sum := Euler_Poincare(flb_pts);
          if sum = 0
           then put_line("  Recomputation is OK...");
           else put_line("  BUG confirmed");
          end if;
        else
          put_line("   OK.");  bug := false;
        end if;
        Standard_Integer_VecVecs.Clear(pts);
        exit when bug;
      end loop;
    end;
    if not bug then
      put("No bugs found, with "); put(times,1);
      put_line(" generated cases tested.");
      put("Dimension : "); put(n,1);
      put(" and #points : "); put(m,1); put_line(".");
    else
      put("Bug found at case "); put(cnt,1); put_line(".");
    end if;
  end Integer_Automatic_Testing;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing the face enumerators by computing f-vectors.");
    loop
      new_line;
      put_line("Choose one of the following :                  ");
      put_line("  0. Exit this program.                        ");
      put_line("  1. f-vector of given integer polytope.       ");
      put_line("  2. f-vector of given floating polytope.      ");
      put_line("  3. f-vector of random integer polytope.      ");
      put_line("  4. f-vector of random floating polytope.     ");
      put("Type 0,1,2,3, or 4 to select : "); get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Integer_Interactive_Testing;
        when '2' => Floating_Interactive_Testing;
        when '3' => Integer_Automatic_Testing;
        when '4' => Floating_Automatic_Testing;
        when others => null;
      end case;
    end loop;
  end Main;

begin 
  Main;
end ts_fvector;
