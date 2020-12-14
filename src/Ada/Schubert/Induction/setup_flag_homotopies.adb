with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Natural_Matrices;         use Standard_Natural_Matrices;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with DoblDobl_Random_Numbers;           use DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;           use QuadDobl_Random_Numbers;
with Standard_Random_Matrices;
with DoblDobl_Random_Matrices;
with QuadDobl_Random_Matrices;
with Symbol_Table;
with Matrix_Indeterminates;
with Checker_Moves;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;
with Symbolic_Schubert_Conditions;

package body Setup_Flag_Homotopies is

  function Random_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    dim : constant natural32 := natural32(n);
   -- res : Standard_Complex_Matrices.Matrix(1..n,1..n);
    res : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Standard_Random_Matrices.Random_Orthogonal_Matrix(dim,dim);

  begin
   -- for j in 1..n loop
   --   for i in 1..n loop
   --     res(i,j) := Random1;
   --   end loop;
   --   for i in 1..(n-j) loop
   --     res(i,j) := Random1;
   --   end loop; 
   --   res(n-j+1,j) := Standard_Complex_Numbers.Create(1.0);
   --   for i in (n-j+2)..n loop
   --     res(i,j) := Standard_Complex_Numbers.Create(0.0);
   --   end loop;
   -- end loop;
    return res;
  end Random_Flag;

  function Random_Flag
             ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix is

    dim : constant natural32 := natural32(n);
    res : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
        := DoblDobl_Random_Matrices.Random_Orthogonal_Matrix(dim,dim);

  begin
    return res;
  end Random_Flag;

  function Random_Flag
             ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix is

    dim : constant natural32 := natural32(n);
    res : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
        := QuadDobl_Random_Matrices.Random_Orthogonal_Matrix(dim,dim);

  begin
    return res;
  end Random_Flag;

  function One_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for j in 1..n loop
      for i in 1..(n-j+1) loop
        res(i,j) := Standard_Complex_Numbers.Create(1.0);
      end loop; 
      for i in (n-j+2)..n loop
        res(i,j) := Standard_Complex_Numbers.Create(0.0);
      end loop;
    end loop;
    return res;
  end One_Flag;

  function One_Flag
             ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);

  begin
    for j in 1..n loop
      for i in 1..(n-j+1) loop
        res(i,j) := DoblDobl_Complex_Numbers.Create(one);
      end loop; 
      for i in (n-j+2)..n loop
        res(i,j) := DoblDobl_Complex_Numbers.Create(zero);
      end loop;
    end loop;
    return res;
  end One_Flag;

  function One_Flag
             ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    for j in 1..n loop
      for i in 1..(n-j+1) loop
        res(i,j) := QuadDobl_Complex_Numbers.Create(one);
      end loop; 
      for i in (n-j+2)..n loop
        res(i,j) := QuadDobl_Complex_Numbers.Create(zero);
      end loop;
    end loop;
    return res;
  end One_Flag;

  function Identity
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := Standard_Complex_Numbers.Create(1.0);
         else res(i,j) := Standard_Complex_Numbers.Create(0.0);
        end if;
      end loop;
    end loop;
    return res;
  end Identity;

  function Identity
             ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := DoblDobl_Complex_Numbers.Create(one);
         else res(i,j) := DoblDobl_Complex_Numbers.Create(zero);
        end if;
      end loop;
    end loop;
    return res;
  end Identity;

  function Identity
             ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if i = j
         then res(i,j) := QuadDobl_Complex_Numbers.Create(one);
         else res(i,j) := QuadDobl_Complex_Numbers.Create(zero);
        end if;
      end loop;
    end loop;
    return res;
  end Identity;

  function Moved_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

    use Standard_Complex_Numbers;

  begin
    for i in 1..n loop
      for j in 1..n-i+1 loop
        if i mod 2 = 0
         then res(i,j) := Create(-1.0);
         else res(i,j) := Create(1.0);
        end if;
      end loop;
      for j in (n-i+2)..n loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    return res;
  end Moved_Flag;

  function Moved_Flag
             ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    minone : constant double_double := create(-1.0);

    use DoblDobl_Complex_Numbers;

  begin
    for i in 1..n loop
      for j in 1..n-i+1 loop
        if i mod 2 = 0
         then res(i,j) := Create(minone);
         else res(i,j) := Create(one);
        end if;
      end loop;
      for j in (n-i+2)..n loop
        res(i,j) := Create(zero);
      end loop;
    end loop;
    return res;
  end Moved_Flag;

  function Moved_Flag
             ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    minone : constant quad_double := create(-1.0);

    use QuadDobl_Complex_Numbers;

  begin
    for i in 1..n loop
      for j in 1..n-i+1 loop
        if i mod 2 = 0
         then res(i,j) := Create(minone);
         else res(i,j) := Create(one);
        end if;
      end loop;
      for j in (n-i+2)..n loop
        res(i,j) := Create(zero);
      end loop;
    end loop;
    return res;
  end Moved_Flag;

  procedure Write_Standard_Moving_Flag
              ( file : in file_type;
                flag : in Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;

    realcff : double_float;
    intcff : integer32;

  begin
    for i in flag'range(1) loop
      for j in flag'range(2) loop
        realcff := REAL_PART(flag(i,j));
        intcff := integer32(realcff);
        put(file,intcff,3);
      end loop;
      new_line(file);
    end loop;
  end Write_Standard_Moving_Flag;

  procedure Write_DoblDobl_Moving_Flag
              ( file : in file_type;
                flag : in DoblDobl_Complex_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;

    realcff : double_double;
    intcff : integer32;

  begin
    for i in flag'range(1) loop
      for j in flag'range(2) loop
        realcff := REAL_PART(flag(i,j));
        intcff := integer32(hi_part(realcff));
        put(file,intcff,3);
      end loop;
      new_line(file);
    end loop;
  end Write_DoblDobl_Moving_Flag;

  procedure Write_QuadDobl_Moving_Flag
              ( file : in file_type;
                flag : in QuadDobl_Complex_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;

    realcff : quad_double;
    intcff : integer32;

  begin
    for i in flag'range(1) loop
      for j in flag'range(2) loop
        realcff := REAL_PART(flag(i,j));
        intcff := integer32(hihi_part(realcff));
        put(file,intcff,3);
      end loop;
      new_line(file);
    end loop;
  end Write_QuadDobl_Moving_Flag;

  function Numeric_Transformation
              ( t : Standard_Natural_Matrices.Matrix;
                g : Standard_Complex_Numbers.Complex_Number )
              return Standard_Complex_Matrices.Matrix is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Matrices.Matrix(t'range(1),t'range(2));

  begin
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Create(0.0);
        elsif t(i,j) = 1 then
          res(i,j) := Create(1.0);
        else
          res(i,j) := g;
        end if;
      end loop;
    end loop;
    return res;
  end Numeric_Transformation;

  function Numeric_Transformation
              ( t : Standard_Natural_Matrices.Matrix;
                g : DoblDobl_Complex_Numbers.Complex_Number )
              return DoblDobl_Complex_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Matrices.Matrix(t'range(1),t'range(2));
    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

  begin
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Create(zero);
        elsif t(i,j) = 1 then
          res(i,j) := Create(one);
        else
          res(i,j) := g;
        end if;
      end loop;
    end loop;
    return res;
  end Numeric_Transformation;

  function Numeric_Transformation
              ( t : Standard_Natural_Matrices.Matrix;
                g : QuadDobl_Complex_Numbers.Complex_Number )
              return QuadDobl_Complex_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Matrices.Matrix(t'range(1),t'range(2));
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Create(zero);
        elsif t(i,j) = 1 then
          res(i,j) := Create(one);
        else
          res(i,j) := g;
        end if;
      end loop;
    end loop;
    return res;
  end Numeric_Transformation;

  function Numeric_Transformation
              ( t : Standard_Natural_Matrices.Matrix )
              return Standard_Complex_Matrices.Matrix is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Matrices.Matrix(t'range(1),t'range(2));

  begin
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Create(0.0);
        elsif t(i,j) = 1 then
          res(i,j) := Create(1.0);
          if i > t'first(1) then
            if t(i-1,j) = 2
             then res(i,j) := Create(-1.0);
            end if;
          end if;
        else
          res(i,j) := Create(1.0);
          if i = t'first(1)  -- special case, see i > t'first(1) above
           then res(i+1,j) := Create(-1.0);
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Numeric_Transformation;

  function Numeric_Transformation
              ( t : Standard_Natural_Matrices.Matrix )
              return DoblDobl_Complex_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Matrices.Matrix(t'range(1),t'range(2));
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    minone : constant double_double := create(-1.0);

  begin
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Create(zero);
        elsif t(i,j) = 1 then
          res(i,j) := Create(one);
          if i > t'first(1) then
            if t(i-1,j) = 2
             then res(i,j) := Create(minone);
            end if;
          end if;
        else
          res(i,j) := Create(one);
          if i = t'first(1)  -- special case, see i > t'first(1) above
           then res(i+1,j) := Create(minone);
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Numeric_Transformation;

  function Numeric_Transformation
              ( t : Standard_Natural_Matrices.Matrix )
              return QuadDobl_Complex_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Matrices.Matrix(t'range(1),t'range(2));
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    minone : constant quad_double := create(-1.0);

  begin
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Create(zero);
        elsif t(i,j) = 1 then
          res(i,j) := Create(one);
          if i > t'first(1) then
            if t(i-1,j) = 2
             then res(i,j) := Create(minone);
            end if;
          end if;
        else
          res(i,j) := Create(one);
          if i = t'first(1)  -- special case, see i > t'first(1) above
           then res(i+1,j) := Create(minone);
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Numeric_Transformation;

  procedure Add_t_Symbol is

    sb : Symbol_Table.Symbol;

  begin
   -- put("#symbols in table : "); put(Symbol_Table.Number,1); new_line;
    Symbol_Table.Enlarge(1);
    sb := (sb'range => ' ');
    sb(sb'first) := 't';
    Symbol_Table.Add(sb);
   -- put("#symbols in table : "); put(Symbol_Table.Number,1); new_line;
  end Add_t_Symbol;

  procedure Initialize_Homotopy_Symbols
             ( dim : in natural32;
               locmap : in Standard_Natural_Matrices.Matrix ) is
  begin
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
    Add_t_Symbol;
  end Initialize_Homotopy_Symbols;

  procedure Insert_Scaling_Symbol ( i,j : in natural32 ) is

    sb : Symbol_Table.Symbol := Matrix_Indeterminates.X_ij(i,j);
    nb : constant natural32 := Symbol_Table.Number;

  begin
    sb(sb'first) := 's';
    Symbol_Table.Replace(nb,sb); -- replace t by sij
    Add_t_Symbol; -- add t again to the symbol table
  end Insert_Scaling_Symbol;

  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Complex_Numbers.Create(-1.0); --Create(1.0);

  begin
    return Symbolic_Transformation(n,v,gamma,t);
  end Symbolic_Transformation;

  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    minone : constant double_double := create(-1.0);
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers.Create(minone);

  begin
    return Symbolic_Transformation(n,v,gamma,t);
  end Symbolic_Transformation;

  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    minone : constant quad_double := create(-1.0);
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers.Create(minone);

  begin
    return Symbolic_Transformation(n,v,gamma,t);
  end Symbolic_Transformation;

  function Symbolic_Transformation
             ( n,v : integer32;
               gamma : Standard_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    trm : Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.cf := Create(1.0);
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif t(i,j) = 1 then
          res(i,j) := Create(trm);
        else
          trm.dg(v) := 1;
          trm.cf := gamma;
          res(i,j) := Create(trm);
          trm.dg(v) := 0;
          trm.cf := Create(1.0);
        end if;
      end loop;
    end loop;
    Clear(trm);
    return res;
  end Symbolic_Transformation;

  function Symbolic_Transformation
             ( n,v : integer32;
               gamma : DoblDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    one : constant double_double := create(1.0);
    trm : Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.cf := Create(one);
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif t(i,j) = 1 then
          res(i,j) := Create(trm);
        else
          trm.dg(v) := 1;
          trm.cf := gamma;
          res(i,j) := Create(trm);
          trm.dg(v) := 0;
          trm.cf := Create(one);
        end if;
      end loop;
    end loop;
    Clear(trm);
    return res;
  end Symbolic_Transformation;

  function Symbolic_Transformation
             ( n,v : integer32;
               gamma : QuadDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    one : constant quad_double := create(1.0);
    trm : Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.cf := Create(one);
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif t(i,j) = 1 then
          res(i,j) := Create(trm);
        else
          trm.dg(v) := 1;
          trm.cf := gamma;
          res(i,j) := Create(trm);
          trm.dg(v) := 0;
          trm.cf := Create(one);
        end if;
      end loop;
    end loop;
    Clear(trm);
    return res;
  end Symbolic_Transformation;

  function Inverse_Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    gamma : constant Complex_Number := Create(1.0);

  begin
    return Inverse_Symbolic_Transformation(n,v,gamma,t);
  end Inverse_Symbolic_Transformation;

  function Inverse_Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    one : constant double_double := create(1.0);
    gamma : constant Complex_Number := Create(one);

  begin
    return Inverse_Symbolic_Transformation(n,v,gamma,t);
  end Inverse_Symbolic_Transformation;

  function Inverse_Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    one : constant quad_double := create(1.0);
    gamma : constant Complex_Number := Create(one);

  begin
    return Inverse_Symbolic_Transformation(n,v,gamma,t);
  end Inverse_Symbolic_Transformation;

  function Inverse_Symbolic_Transformation
             ( n,v : integer32;
               gamma : Standard_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    trm : Term;
    i2,j2 : integer32;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.cf := Create(1.0);
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif t(i,j) = 1 then
          res(i,j) := Create(trm);
        else
          i2 := i; j2 := j;
          res(i,j) := Null_Poly;
        end if;
      end loop;
    end loop;
    trm.dg(v) := 1;
    trm.cf := -gamma;
    res(i2+1,j2+1) := Create(trm);
    Clear(trm);
    return res;
  end Inverse_Symbolic_Transformation;

  function Inverse_Symbolic_Transformation
             ( n,v : integer32;
               gamma : DoblDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    one : constant double_double := create(1.0);
    trm : Term;
    i2,j2 : integer32;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.cf := Create(one);
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif t(i,j) = 1 then
          res(i,j) := Create(trm);
        else
          i2 := i; j2 := j;
          res(i,j) := Null_Poly;
        end if;
      end loop;
    end loop;
    trm.dg(v) := 1;
    trm.cf := -gamma;
    res(i2+1,j2+1) := Create(trm);
    Clear(trm);
    return res;
  end Inverse_Symbolic_Transformation;

  function Inverse_Symbolic_Transformation
             ( n,v : integer32;
               gamma : QuadDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    one : constant quad_double := create(1.0);
    trm : Term;
    i2,j2 : integer32;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.cf := Create(one);
    for i in t'range(1) loop
      for j in t'range(2) loop
        if t(i,j) = 0 then
          res(i,j) := Null_Poly;
        elsif t(i,j) = 1 then
          res(i,j) := Create(trm);
        else
          i2 := i; j2 := j;
          res(i,j) := Null_Poly;
        end if;
      end loop;
    end loop;
    trm.dg(v) := 1;
    trm.cf := -gamma;
    res(i2+1,j2+1) := Create(trm);
    Clear(trm);
    return res;
  end Inverse_Symbolic_Transformation;

  function Evaluate_Transformation
             ( t : Standard_Complex_Poly_Matrices.Matrix;
               v : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    n : integer32;
    trm,nrt : Term;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if t(i,j) = Null_Poly then
          res(i,j) := Null_Poly;
        else
          trm := Head(t(i,j));
          n := trm.dg'last;
          nrt.dg := new Standard_Natural_Vectors.Vector'(1..n-1 => 0);
          if trm.dg(n) = 0
           then nrt.cf := trm.cf;
           else nrt.cf := trm.cf*v;
          end if;
          res(i,j) := Create(nrt);
          Clear(nrt);
        end if;
      end loop;
    end loop;
    return res;
  end Evaluate_Transformation;

  function Evaluate_Transformation
             ( t : DoblDobl_Complex_Poly_Matrices.Matrix;
               v : DoblDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    n : integer32;
    trm,nrt : Term;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if t(i,j) = Null_Poly then
          res(i,j) := Null_Poly;
        else
          trm := Head(t(i,j));
          n := trm.dg'last;
          nrt.dg := new Standard_Natural_Vectors.Vector'(1..n-1 => 0);
          if trm.dg(n) = 0
           then nrt.cf := trm.cf;
           else nrt.cf := trm.cf*v;
          end if;
          res(i,j) := Create(nrt);
          Clear(nrt);
        end if;
      end loop;
    end loop;
    return res;
  end Evaluate_Transformation;

  function Evaluate_Transformation
             ( t : QuadDobl_Complex_Poly_Matrices.Matrix;
               v : QuadDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Poly_Matrices.Matrix(t'range(1),t'range(2));
    n : integer32;
    trm,nrt : Term;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if t(i,j) = Null_Poly then
          res(i,j) := Null_Poly;
        else
          trm := Head(t(i,j));
          n := trm.dg'last;
          nrt.dg := new Standard_Natural_Vectors.Vector'(1..n-1 => 0);
          if trm.dg(n) = 0
           then nrt.cf := trm.cf;
           else nrt.cf := trm.cf*v;
          end if;
          res(i,j) := Create(nrt);
          Clear(nrt);
        end if;
      end loop;
    end loop;
    return res;
  end Evaluate_Transformation;

  function Moving_Flag
             ( f : Standard_Complex_Matrices.Matrix;
               t : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Poly_Matrices.Matrix(f'range(1),t'range(2));
    zero : constant Complex_Number := Create(0.0);
    one : constant Complex_Number := Create(1.0);
    acc : Poly;

  begin
    for i in f'range(1) loop
      for j in t'range(2) loop
        res(i,j) := Null_Poly;
        for k in f'range(2) loop
          if f(i,k) /= zero and t(k,j) /= Null_Poly then
            if f(i,k) = one then
              Add(res(i,j),t(k,j));
            else
              acc := f(i,k)*t(k,j);
              Add(res(i,j),acc);
              Clear(acc);
            end if;
          end if;
        end loop;
      end loop;
    end loop;
    return res;
  end Moving_Flag;

  function Moving_Flag
             ( f : DoblDobl_Complex_Matrices.Matrix;
               t : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Poly_Matrices.Matrix(f'range(1),t'range(2));
    dd_zero : constant double_double := create(0.0);
    dd_one : constant double_double := create(1.0);
    zero : constant Complex_Number := Create(dd_zero);
    one : constant Complex_Number := Create(dd_one);
    acc : Poly;

  begin
    for i in f'range(1) loop
      for j in t'range(2) loop
        res(i,j) := Null_Poly;
        for k in f'range(2) loop
          if f(i,k) /= zero and t(k,j) /= Null_Poly then
            if f(i,k) = one then
              Add(res(i,j),t(k,j));
            else
              acc := f(i,k)*t(k,j);
              Add(res(i,j),acc);
              Clear(acc);
            end if;
          end if;
        end loop;
      end loop;
    end loop;
    return res;
  end Moving_Flag;

  function Moving_Flag
             ( f : QuadDobl_Complex_Matrices.Matrix;
               t : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Poly_Matrices.Matrix(f'range(1),t'range(2));
    qd_zero : constant quad_double := create(0.0);
    qd_one : constant quad_double := create(1.0);
    zero : constant Complex_Number := Create(qd_zero);
    one : constant Complex_Number := Create(qd_one);
    acc : Poly;

  begin
    for i in f'range(1) loop
      for j in t'range(2) loop
        res(i,j) := Null_Poly;
        for k in f'range(2) loop
          if f(i,k) /= zero and t(k,j) /= Null_Poly then
            if f(i,k) = one then
              Add(res(i,j),t(k,j));
            else
              acc := f(i,k)*t(k,j);
              Add(res(i,j),acc);
              Clear(acc);
            end if;
          end if;
        end loop;
      end loop;
    end loop;
    return res;
  end Moving_Flag;

  procedure Move_to_Opposite_Flag
             ( f : in out Standard_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in Standard_Complex_Matrices.Matrix ) is

    n : constant integer32 := q'last;
    fc : constant integer32 := Checker_Moves.Falling_Checker(p);
    ac : constant integer32 := Checker_Moves.Ascending_Checker(p,fc);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(fc)));
    gamma : constant Standard_Complex_Numbers.Complex_Number := mf(n+1-ac,fc);
    nt : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Numeric_Transformation(t,gamma);

    use Standard_Complex_Matrices;

  begin
    f := f*nt;
  end Move_to_Opposite_Flag;

  procedure Move_to_Opposite_Flag
             ( f : in out DoblDobl_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in DoblDobl_Complex_Matrices.Matrix ) is

    n : constant integer32 := q'last;
    fc : constant integer32 := Checker_Moves.Falling_Checker(p);
    ac : constant integer32 := Checker_Moves.Ascending_Checker(p,fc);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(fc)));
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number := mf(n+1-ac,fc);
    nt : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
       := Numeric_Transformation(t,gamma);

    use DoblDobl_Complex_Matrices;

  begin
    f := f*nt;
  end Move_to_Opposite_Flag;

  procedure Move_to_Opposite_Flag
             ( f : in out QuadDobl_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in QuadDobl_Complex_Matrices.Matrix ) is

    n : constant integer32 := q'last;
    fc : constant integer32 := Checker_Moves.Falling_Checker(p);
    ac : constant integer32 := Checker_Moves.Ascending_Checker(p,fc);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(fc)));
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number := mf(n+1-ac,fc);
    nt : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
       := Numeric_Transformation(t,gamma);

    use QuadDobl_Complex_Matrices;

  begin
    f := f*nt;
  end Move_to_Opposite_Flag;

  function Symbolic_Plane
             ( n,k : integer32; p,rows,cols : Standard_Natural_Vectors.Vector )
             return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    xmp : Standard_Natural_Matrices.Matrix(1..n,1..k);

  begin
    xmp := Column_Pattern(n,k,p,rows,cols);
    res := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,xmp);
    return res;
  end Symbolic_Plane;

  function Symbolic_Plane
             ( n,k : integer32; p,rows,cols : Standard_Natural_Vectors.Vector )
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    res : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    xmp : Standard_Natural_Matrices.Matrix(1..n,1..k);

  begin
    xmp := Column_Pattern(n,k,p,rows,cols);
    res := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,xmp);
    return res;
  end Symbolic_Plane;

  function Symbolic_Plane
             ( n,k : integer32; p,rows,cols : Standard_Natural_Vectors.Vector )
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    res : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    xmp : Standard_Natural_Matrices.Matrix(1..n,1..k);

  begin
    xmp := Column_Pattern(n,k,p,rows,cols);
    res := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,xmp);
    return res;
  end Symbolic_Plane;

  function Filter_Zero_Equations
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    res : Poly_Sys(p'range);
    cnt : integer32 := p'first-1;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        cnt := cnt + 1;
        res(cnt) := p(i);
      end if;
    end loop;
    return res(res'first..cnt);
  end Filter_Zero_Equations;

  function Filter_Zero_Equations
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    res : Poly_Sys(p'range);
    cnt : integer32 := p'first-1;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        cnt := cnt + 1;
        res(cnt) := p(i);
      end if;
    end loop;
    return res(res'first..cnt);
  end Filter_Zero_Equations;

  function Filter_Zero_Equations
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    res : Poly_Sys(p'range);
    cnt : integer32 := p'first-1;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        cnt := cnt + 1;
        res(cnt) := p(i);
      end if;
    end loop;
    return res(res'first..cnt);
  end Filter_Zero_Equations;

  function Square ( n : integer32;
                    p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    res : Poly_Sys(1..n);
    acc : Poly;

  begin
   -- put("in Square, p'last = "); put(p'last,1);
   -- put("  n = "); put(n,1); new_line;
    for i in res'range loop
      Copy(p(i),res(i));
      for j in n+1..p'last loop
        acc := Random1*p(j);
        Add(res(i),acc);
        Clear(acc);
      end loop;
    end loop;
    return res;
  end Square;

  function Square ( n : integer32;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    res : Poly_Sys(1..n);
    acc : Poly;

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in n+1..p'last loop
        acc := Random1*p(j);
        Add(res(i),acc);
        Clear(acc);
      end loop;
    end loop;
    return res;
  end Square;

  function Square ( n : integer32;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    res : Poly_Sys(1..n);
    acc : Poly;

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in n+1..p'last loop
        acc := Random1*p(j);
        Add(res(i),acc);
        Clear(acc);
      end loop;
    end loop;
    return res;
  end Square;

  function Concatenate
             ( s : Standard_Complex_Poly_Systems.Array_of_Poly_Sys )
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys is

    use Standard_Complex_Poly_Systems;

    res : Link_to_Poly_Sys;
    cnt : integer32 := 0;
    ind : integer32 := 0;

  begin
    for i in s'range loop
      cnt := cnt + s(i)'last;
    end loop;
    res := new Poly_Sys(1..cnt);
    for i in s'range loop
      for j in s(i)'range loop
        ind := ind + 1;
        res(ind) := s(i)(j);
      end loop;    
    end loop;
    return res;
  end Concatenate;

  function Concatenate
             ( s : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;

    res : Link_to_Poly_Sys;
    cnt : integer32 := 0;
    ind : integer32 := 0;

  begin
    for i in s'range loop
      cnt := cnt + s(i)'last;
    end loop;
    res := new Poly_Sys(1..cnt);
    for i in s'range loop
      for j in s(i)'range loop
        ind := ind + 1;
        res(ind) := s(i)(j);
      end loop;    
    end loop;
    return res;
  end Concatenate;

  function Concatenate
             ( s : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;

    res : Link_to_Poly_Sys;
    cnt : integer32 := 0;
    ind : integer32 := 0;

  begin
    for i in s'range loop
      cnt := cnt + s(i)'last;
    end loop;
    res := new Poly_Sys(1..cnt);
    for i in s'range loop
      for j in s(i)'range loop
        ind := ind + 1;
        res(ind) := s(i)(j);
      end loop;    
    end loop;
    return res;
  end Concatenate;

  procedure Append
              ( s : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                p : in Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Poly_Systems;

    res : Link_to_Poly_Sys;

  begin
    res := new Poly_Sys(s'first..s'last+1); 
    res(s'range) := s.all;
    res(res'last) := p;
    s := res;
  end Append;

  procedure Append
              ( s : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                p : in DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Poly_Systems;

    res : Link_to_Poly_Sys;

  begin
    res := new Poly_Sys(s'first..s'last+1); 
    res(s'range) := s.all;
    res(res'last) := p;
    s := res;
  end Append;

  procedure Append
              ( s : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                p : in QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Poly_Systems;

    res : Link_to_Poly_Sys;

  begin
    res := new Poly_Sys(s'first..s'last+1); 
    res(s'range) := s.all;
    res(res'last) := p;
    s := res;
  end Append;

end Setup_Flag_Homotopies;
