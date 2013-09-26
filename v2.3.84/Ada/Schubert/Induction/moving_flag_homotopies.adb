with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Natural_Matrices;         use Standard_Natural_Matrices;
--with Standard_Natural_Matrices_io;      use Standard_Natural_Matrices_io;
--with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Poly_Matrices_io; use Standard_Complex_Poly_Matrices_io;

with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Norms_Equals;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Symbol_Table;
with Matrix_Indeterminates;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Embed_Polynomials;
with Planes_and_Polynomials;
with Checker_Moves;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;
with Symbolic_Schubert_Conditions;
with Numeric_Schubert_Conditions;

package body Moving_Flag_Homotopies is

  function Random_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for j in 1..n loop
      for i in 1..(n-j) loop
        res(i,j) := Random1;
      end loop; 
      res(n-j+1,j) := Create(1.0);
      for i in (n-j+2)..n loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    return res;
  end Random_Flag;

  function One_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for j in 1..n loop
      for i in 1..(n-j+1) loop
        res(i,j) := Create(1.0);
      end loop; 
      for i in (n-j+2)..n loop
        res(i,j) := Create(0.0);
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
         then res(i,j) := Create(1.0);
         else res(i,j) := Create(0.0);
        end if;
      end loop;
    end loop;
    return res;
  end Identity;

  function Numeric_Transformation
              ( t : Standard_Natural_Matrices.Matrix; g : Complex_Number )
              return Standard_Complex_Matrices.Matrix is

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

  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    gamma : constant Complex_Number := Create(1.0);

  begin
    return Symbolic_Transformation(n,v,gamma,t);
  end Symbolic_Transformation;

  function Symbolic_Transformation
             ( n,v : integer32; gamma : Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

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

  function Inverse_Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

    gamma : constant Complex_Number := Create(1.0);

  begin
    return Inverse_Symbolic_Transformation(n,v,gamma,t);
  end Inverse_Symbolic_Transformation;

  function Inverse_Symbolic_Transformation
             ( n,v : integer32; gamma : Complex_Number;
               t : Standard_Natural_Matrices.Matrix ) 
             return Standard_Complex_Poly_Matrices.Matrix is

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

  function Evaluate_Transformation
             ( t : Standard_Complex_Poly_Matrices.Matrix; v : Complex_Number )
             return Standard_Complex_Poly_Matrices.Matrix is

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

  function Moving_Flag
             ( f : Standard_Complex_Matrices.Matrix;
               t : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix is

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

  procedure Move_to_Opposite_Flag
             ( f : in out Standard_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in Standard_Complex_Matrices.Matrix ) is

    n : constant integer32 := q'last;
    fc : constant integer32 := Checker_Moves.Falling_Checker(p);
    ac : constant integer32 := Checker_Moves.Ascending_Checker(p,fc);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(fc)));
    gamma : constant Complex_Number := mf(n+1-ac,fc);
    nt : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Numeric_Transformation(t,gamma);

    use Standard_Complex_Matrices;

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

  function Filter_Zero_Equations ( p : Poly_Sys ) return Poly_Sys is

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

  function Square ( n : integer32; p : Poly_Sys ) return Poly_Sys is

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

  procedure One_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys ) is

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    dim : constant natural32 := Degree_of_Freedom(locmap);
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),ic));
    sc : Poly_Sys(1..nq);

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := Standard_Embed_Polynomials.Add_Variables(x,1);
   -- nt := Numeric_Transformation(t,gamma);
    st := Symbolic_Transformation(integer32(dim)+1,integer32(dim)+1,gamma,t);
    sf := Moving_Flag(nf,st);
    xp := sf*xt;
    sc := Numeric_Schubert_Conditions.Expand(n,k,nq,ic,xp,vf);
   -- nf := nf*nt;
    h := new Poly_Sys'(Filter_Zero_Equations(sc));
  end One_Flag_Homotopy;

  procedure One_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys ) is

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
          -- := Column_Pattern(n,k,p,rows,cols); -- change from p
           := Column_Pattern(n,k,q,rows,cols);   -- to current q
    dim : constant integer32 := integer32(Degree_of_Freedom(locmap));
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),ic));
    sc : Poly_Sys(1..nq);

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Matrices;

  begin
   -- put(file,"dimension : "); put(file,dim,1); new_line(file);
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := Standard_Embed_Polynomials.Add_Variables(x,1);
   -- put_line(file,"Localization pattern at current node :");
   -- put(file,Column_Pattern(n,k,q,rows,cols));
    put_line(file,"The symbolic plane : "); put(file,xt);
    st := Symbolic_Transformation(dim+1,dim+1,gamma,t);
   -- nt := Numeric_Transformation(t,gamma);
   -- put_line(file,"The numeric transformation : "); put(file,nt,3);
   -- put_line(file,"The symbolic transformation : "); put(file,st);
    sf := Moving_Flag(nf,st);
    put_line(file,"Symbolic form of the moving flag : "); put(file,sf);
    xp := sf*xt;
    put_line(file,"moving coordinates : "); put(file,xp);
    sc := Numeric_Schubert_Conditions.Expand(n,k,nq,ic,xp,vf);
   -- nf := nf*nt;
    h := new Poly_Sys'(Filter_Zero_Equations(sc));
  end One_Flag_Homotopy;

  function Concatenate ( s : Array_of_Poly_Sys ) return Link_to_Poly_Sys is

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

  procedure Flag_Conditions
             ( n,k : in integer32;
               x : in Standard_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Link_to_Poly_Sys ) is

    s : Array_of_Poly_Sys(ic'range);

  begin
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Flag_Conditions;

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Link_to_Poly_Sys ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Flag_Conditions;

  procedure Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               f : out Link_to_Poly_Sys ) is

   -- fc : constant natural := Checker_Moves.Falling_Checker(p);
   -- ac : constant natural := Checker_Moves.Ascending_Checker(p,fc);
   -- t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
   --   := Transformation(n,q(fc));
   -- gamma : constant Complex_Number := Create(0.0); -- start solution !
    nt : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
      -- := Numeric_Transformation(t,gamma);
       := Identity(n);
    st : Standard_Complex_Matrices.Matrix(1..n,1..n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
   -- dim : constant natural := Degree_of_Freedom(locmap);
    x,xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    st := nf*nt;
    put_line(file,"Moving flag in intersection conditions : ");
    put(file,st,2);
    xp := Moving_Flag(st,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,xp,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Flag_Conditions;

  procedure Moving_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys ) is

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    dim : constant natural32 := Degree_of_Freedom(locmap);
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := Standard_Embed_Polynomials.Add_Variables(x,1);
   -- nt := Numeric_Transformation(t,gamma);
    st := Symbolic_Transformation(integer32(dim)+1,integer32(dim)+1,gamma,t);
    sf := Moving_Flag(nf,st);
    xp := sf*xt;
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,xp,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
   -- nf := nf*nt;
    h := Concatenate(s);
  end Moving_Flag_Homotopy;

  procedure Moving_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys ) is

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols); -- change from p
          -- rows and cols must be with p
          -- := Column_Pattern(n,k,p,rows,cols); -- change from p
          -- := Column_Pattern(n,k,q,rows,cols);   -- to current q
    dim : constant natural32 := Degree_of_Freedom(locmap);
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : Standard_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Matrices;

  begin
    put(file,"dimension : "); put(file,dim,1); new_line(file);
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := Standard_Embed_Polynomials.Add_Variables(x,1);
   -- put_line(file,"Localization pattern at current node :");
   -- put(file,Column_Pattern(n,k,q,rows,cols));
   -- put_line(file,"The symbolic plane : "); put(file,xt);
    st := Symbolic_Transformation(integer32(dim)+1,integer32(dim)+1,gamma,t);
   -- nt := Numeric_Transformation(t,gamma);
   -- put_line(file,"The numeric transformation : "); put(file,nt,3);
   -- put_line(file,"The symbolic transformation : "); put(file,st);
    sf := Moving_Flag(nf,st);
   -- put_line(file,"Symbolic form of the moving flag : "); put(file,sf);
    xp := sf*xt;
   -- put_line(file,"moving coordinates : "); put(file,xp);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,xp,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    h := Concatenate(s);
   -- nf := nf*nt;
  end Moving_Flag_Homotopy;

  function Inconsistent ( p : Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if degree(p(i)) = 0
         then return true;
        end if;
      end if;
    end loop;
    return false;
  end Inconsistent;

  function Linear_Equations ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);
    cnt : integer32 := res'first-1;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if Degree(p(i)) = 1 then
          cnt := cnt + 1;
          res(cnt) := p(i);
        end if;
      end if;
    end loop;
    return res(res'first..cnt);
  end Linear_Equations;

  procedure Coefficients
             ( p : in Poly_Sys; 
               A : out Standard_Complex_Matrices.Matrix;
               b : out Standard_Complex_Vectors.Vector ) is

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    h : Standard_Complex_Vectors.Vector(0..integer32(n));

  begin
    for i in p'range loop
      h := Planes_and_Polynomials.Polynomial(p(i));
      b(i) := -h(0);
      for j in A'range(2) loop
        A(i,j) := h(j);
      end loop;
    end loop;
  end Coefficients;

  procedure Solve ( A : in out Standard_Complex_Matrices.Matrix;
                    b : in Standard_Complex_Vectors.Vector;
                    x : out Standard_Complex_Vectors.Vector;
                    res : out double_float ) is

  -- DESCRIPTION :
  --   Applies singular value decomposition to A to solve A*x = b.

    AA : constant Standard_Complex_Matrices.Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    m : constant integer32 := Min0(n+1,p);
    s : Standard_Complex_Vectors.Vector(1..m);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Standard_Complex_Matrices.Matrix(1..n,1..n);
    v : Standard_Complex_Matrices.Matrix(1..p,1..p);
    info : integer32;
    r : Standard_Complex_Vectors.Vector(b'range);
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

  begin
    SVD(A,n,p,s,e,u,v,11,info);
    x := Solve(u,v,s,b);
    r := b - AA*x;
   -- put_line("the singular values : "); put_line(s);
   -- put_line("The residual :"); put_line(r);
    res := Standard_Complex_Norms_Equals.Norm2(r);
  end Solve;

  procedure First_Solution
             ( f : in Poly_Sys; fail : out boolean;
               x : out Standard_Complex_Vectors.Vector;
               res : out double_float ) is

    n : constant integer32 := x'last;

  begin
    fail := true;
    x := (x'range => Create(0.0));
    if not Inconsistent(f) then
      declare
        s : constant Poly_Sys := Linear_Equations(f);
        m : constant integer32 := s'last;
        A : Standard_Complex_Matrices.Matrix(1..m,1..n);
        b : Standard_Complex_Vectors.Vector(1..m);
      begin
        if s'last >= n then
          Coefficients(s,A,b);
         -- put_line("The coefficient matrix : "); put(A,3);
         -- put_line("The righthandside vector : "); put_line(b);
          Solve(A,b,x,res);
          fail := (res > 1.0E-8);
        end if;
      end;
    end if;
  end First_Solution;

  procedure Start_Solution
             ( h : in Poly_Sys; fail : out boolean;
               x : out Standard_Complex_Vectors.Vector;
               res : out double_float ) is

    n : constant integer32 := x'last;
    h0 : Poly_Sys(h'range) := Eval(h,Create(0.0),n+1);
      -- h0 is system h where t is substituted by zero

  begin
    First_Solution(h0,fail,x,res);
    Clear(h0);
  end Start_Solution;

  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix is
 
    res : Standard_Complex_Poly_Matrices.Matrix(start'range(1),start'range(2));
    z : Complex_Number;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nv+1 => 0);
    for i in start'range(1) loop
      for j in start'range(2) loop
        z := target(i,j);
        if ((REAL_PART(z) = 0.0) and (IMAG_PART(z) = 0.0)) then
          res(i,j) := Null_Poly;
        else
          t.cf := z;
          t.dg(nv+1) := 1;
          res(i,j) := Create(t);
        end if;
        z := start(i,j);
        if ((REAL_PART(z) /= 0.0) or (IMAG_PART(z) /= 0.0)) then
          t.cf := z;
          t.dg(nv+1) := 0;
          Add(res(i,j),t);
          t.dg(nv+1) := 1;
          t.cf := -z;
          Add(res(i,j),t);
        end if;
      end loop;
    end loop;
    Clear(t);
    return res;
  end Cheater_Homotopy_Flag;

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               cond : in Bracket;
               start,target : in Standard_Complex_Matrices.Matrix;
               f : out Link_to_Poly_Sys ) is

    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),cond));
    sc : Poly_Sys(1..n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    nv : constant integer32 := integer32(Degree_of_Freedom(locmap));
    x,xt : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    flag : constant Standard_Complex_Poly_Matrices.Matrix(1..n,1..n)
         := Cheater_Homotopy_Flag(nv,start,target);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := Standard_Embed_Polynomials.Add_Variables(x,1);
    sc := Numeric_Schubert_Conditions.Expand(n,k,nq,cond,xt,flag);
    f := new Poly_Sys'(Filter_Zero_Equations(sc));
  end Flag_Conditions;

end Moving_Flag_Homotopies;
