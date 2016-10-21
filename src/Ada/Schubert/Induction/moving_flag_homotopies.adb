with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Matrices_io;      use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Matrices_io;      use QuadDobl_Complex_Matrices_io;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Standard_Complex_Poly_Matrices_io; use Standard_Complex_Poly_Matrices_io;
with DoblDobl_Complex_Poly_Matrices_io; use DoblDobl_Complex_Poly_Matrices_io;
with QuadDobl_Complex_Poly_Matrices_io; use QuadDobl_Complex_Poly_Matrices_io;
with Standard_Embed_Polynomials;
with DoblDobl_Embed_Polynomials;
with QuadDobl_Embed_Polynomials;
with Planes_and_Polynomials;
with Checker_Moves;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;
with Symbolic_Schubert_Conditions;
with Numeric_Schubert_Conditions;
with Setup_Flag_Homotopies;             use Setup_Flag_Homotopies;

package body Moving_Flag_Homotopies is

  procedure One_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Standard_Complex_Numbers.Complex_Number := mf(n+1-a,f);
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
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    dim : constant natural32 := Degree_of_Freedom(locmap);
    st,sf : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),ic));
    sc : Poly_Sys(1..nq);

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := DoblDobl_Embed_Polynomials.Add_Variables(x,1);
   -- nt := Numeric_Transformation(t,gamma);
    st := Symbolic_Transformation(integer32(dim)+1,integer32(dim)+1,gamma,t);
    sf := Moving_Flag(nf,st);
    xp := sf*xt;
    sc := Numeric_Schubert_Conditions.Expand(n,k,nq,ic,xp,vf);
   -- nf := nf*nt;
    h := new Poly_Sys'(Filter_Zero_Equations(sc));
  end One_Flag_Homotopy;

  procedure One_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    dim : constant natural32 := Degree_of_Freedom(locmap);
    st,sf : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),ic));
    sc : Poly_Sys(1..nq);

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := QuadDobl_Embed_Polynomials.Add_Variables(x,1);
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
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Standard_Complex_Numbers.Complex_Number := mf(n+1-a,f);
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

  procedure One_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
          -- := Column_Pattern(n,k,p,rows,cols); -- change from p
           := Column_Pattern(n,k,q,rows,cols);   -- to current q
    dim : constant integer32 := integer32(Degree_of_Freedom(locmap));
   -- nt : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),ic));
    sc : Poly_Sys(1..nq);

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Matrices;

  begin
   -- put(file,"dimension : "); put(file,dim,1); new_line(file);
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := DoblDobl_Embed_Polynomials.Add_Variables(x,1);
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

  procedure One_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
          -- := Column_Pattern(n,k,p,rows,cols); -- change from p
           := Column_Pattern(n,k,q,rows,cols);   -- to current q
    dim : constant integer32 := integer32(Degree_of_Freedom(locmap));
   -- nt : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),ic));
    sc : Poly_Sys(1..nq);

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Matrices;

  begin
   -- put(file,"dimension : "); put(file,dim,1); new_line(file);
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := QuadDobl_Embed_Polynomials.Add_Variables(x,1);
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

  procedure Flag_Conditions
             ( n,k : in integer32;
               x : in Standard_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

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
               x : in DoblDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

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
               x : in QuadDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

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

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               x : in Standard_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    s : Array_of_Poly_Sys(ic'range);

  begin
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               x : in DoblDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    s : Array_of_Poly_Sys(ic'range);

  begin
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               x : in QuadDobl_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    s : Array_of_Poly_Sys(ic'range);

  begin
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

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
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
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
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
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

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

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
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,x,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in Standard_Complex_Matrices.Matrix;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x,mfx : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    mfx := Moving_Flag(mf,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,mfx,vf(i).all);
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
               mf : in DoblDobl_Complex_Matrices.Matrix;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x,mfx : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    mfx := Moving_Flag(mf,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,mfx,vf(i).all);
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
               mf : in QuadDobl_Complex_Matrices.Matrix;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x,mfx : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    mfx := Moving_Flag(mf,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Expand(n,k,nq,c,mfx,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in Standard_Complex_Matrices.Matrix;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x,mfx : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    mfx := Moving_Flag(mf,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand
               (n,k,nq,c,mfx,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in DoblDobl_Complex_Matrices.Matrix;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x,mfx : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    mfx := Moving_Flag(mf,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand
               (n,k,nq,c,mfx,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               mf : in QuadDobl_Complex_Matrices.Matrix;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    x,mfx : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    mfx := Moving_Flag(mf,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand
               (n,k,nq,c,mfx,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

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
    st := mf; -- nf*nt;
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

  procedure Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

   -- fc : constant natural := Checker_Moves.Falling_Checker(p);
   -- ac : constant natural := Checker_Moves.Ascending_Checker(p,fc);
   -- t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
   --   := Transformation(n,q(fc));
   -- gamma : constant Complex_Number := Create(0.0); -- start solution !
    nt : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
      -- := Numeric_Transformation(t,gamma);
       := Identity(n);
    st : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
   -- dim : constant natural := Degree_of_Freedom(locmap);
    x,xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    st := mf; -- nf*nt;
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

  procedure Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

   -- fc : constant natural := Checker_Moves.Falling_Checker(p);
   -- ac : constant natural := Checker_Moves.Ascending_Checker(p,fc);
   -- t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
   --   := Transformation(n,q(fc));
   -- gamma : constant Complex_Number := Create(0.0); -- start solution !
    nt : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
      -- := Numeric_Transformation(t,gamma);
       := Identity(n);
    st : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
   -- dim : constant natural := Degree_of_Freedom(locmap);
    x,xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    st := mf; -- nf*nt;
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

  procedure Minimal_Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

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
    st := mf; -- nf*nt;
    put_line(file,"Moving flag in intersection conditions : ");
    put(file,st,2);
    xp := Moving_Flag(st,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,xp,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

   -- fc : constant natural := Checker_Moves.Falling_Checker(p);
   -- ac : constant natural := Checker_Moves.Ascending_Checker(p,fc);
   -- t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
   --   := Transformation(n,q(fc));
   -- gamma : constant Complex_Number := Create(0.0); -- start solution !
    nt : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
      -- := Numeric_Transformation(t,gamma);
       := Identity(n);
    st : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
   -- dim : constant natural := Degree_of_Freedom(locmap);
    x,xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    st := mf; -- nf*nt;
    put_line(file,"Moving flag in intersection conditions : ");
    put(file,st,2);
    xp := Moving_Flag(st,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,xp,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Minimal_Flag_Conditions
             ( file: in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

   -- fc : constant natural := Checker_Moves.Falling_Checker(p);
   -- ac : constant natural := Checker_Moves.Ascending_Checker(p,fc);
   -- t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
   --   := Transformation(n,q(fc));
   -- gamma : constant Complex_Number := Create(0.0); -- start solution !
    nt : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
      -- := Numeric_Transformation(t,gamma);
       := Identity(n);
    st : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
   -- dim : constant natural := Degree_of_Freedom(locmap);
    x,xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    st := mf; -- nf*nt;
    put_line(file,"Moving flag in intersection conditions : ");
    put(file,st,2);
    xp := Moving_Flag(st,x);
    for i in ic'range loop
      declare
        c : constant Brackets.Bracket(1..k) := Brackets.Bracket(ic(i).all);
        nq : constant integer32
           := integer32(Symbolic_Schubert_Conditions.Number_of_NotAbove
                          (natural32(n),c));
        sc : constant Poly_Sys(1..nq)
           := Numeric_Schubert_Conditions.Minimal_Expand(n,k,nq,c,xp,vf(i).all);
      begin
        s(i) := new Poly_Sys'(Filter_Zero_Equations(sc));
      end;
    end loop;
    f := Concatenate(s);
  end Minimal_Flag_Conditions;

  procedure Moving_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Standard_Complex_Numbers.Complex_Number := mf(n+1-a,f);
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
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    dim : constant natural32 := Degree_of_Freedom(locmap);
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := DoblDobl_Embed_Polynomials.Add_Variables(x,1);
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
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    dim : constant natural32 := Degree_of_Freedom(locmap);
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Matrices;

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := QuadDobl_Embed_Polynomials.Add_Variables(x,1);
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
               h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant Standard_Complex_Numbers.Complex_Number := mf(n+1-a,f);
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

  procedure Moving_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in DoblDobl_Complex_VecMats.VecMat;
               mf,nf : in DoblDobl_Complex_Matrices.Matrix;
               h : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols); -- change from p
          -- rows and cols must be with p
          -- := Column_Pattern(n,k,p,rows,cols); -- change from p
          -- := Column_Pattern(n,k,q,rows,cols);   -- to current q
    dim : constant natural32 := Degree_of_Freedom(locmap);
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Matrices;

  begin
    put(file,"dimension : "); put(file,dim,1); new_line(file);
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := DoblDobl_Embed_Polynomials.Add_Variables(x,1);
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

  procedure Moving_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in QuadDobl_Complex_VecMats.VecMat;
               mf,nf : in QuadDobl_Complex_Matrices.Matrix;
               h : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    f : constant integer32 := Checker_Moves.Falling_Checker(p);
    a : constant integer32 := Checker_Moves.Ascending_Checker(p,f);
    t : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Transformation(n,integer32(q(f)));
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number := mf(n+1-a,f);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols); -- change from p
          -- rows and cols must be with p
          -- := Column_Pattern(n,k,p,rows,cols); -- change from p
          -- := Column_Pattern(n,k,q,rows,cols);   -- to current q
    dim : constant natural32 := Degree_of_Freedom(locmap);
   -- nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    st,sf : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..n);
    x,xt,xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    s : Array_of_Poly_Sys(ic'range);

    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Matrices;

  begin
    put(file,"dimension : "); put(file,dim,1); new_line(file);
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := QuadDobl_Embed_Polynomials.Add_Variables(x,1);
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

  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;
 
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

  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;
 
    res : DoblDobl_Complex_Poly_Matrices.Matrix(start'range(1),start'range(2));
    z : Complex_Number;
    zero : constant double_double := create(0.0);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nv+1 => 0);
    for i in start'range(1) loop
      for j in start'range(2) loop
        z := target(i,j);
        if ((REAL_PART(z) = zero) and (IMAG_PART(z) = zero)) then
          res(i,j) := Null_Poly;
        else
          t.cf := z;
          t.dg(nv+1) := 1;
          res(i,j) := Create(t);
        end if;
        z := start(i,j);
        if ((REAL_PART(z) /= zero) or (IMAG_PART(z) /= zero)) then
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

  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;
 
    res : QuadDobl_Complex_Poly_Matrices.Matrix(start'range(1),start'range(2));
    z : Complex_Number;
    zero : constant quad_double := create(0.0);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nv+1 => 0);
    for i in start'range(1) loop
      for j in start'range(2) loop
        z := target(i,j);
        if ((REAL_PART(z) = zero) and (IMAG_PART(z) = zero)) then
          res(i,j) := Null_Poly;
        else
          t.cf := z;
          t.dg(nv+1) := 1;
          res(i,j) := Create(t);
        end if;
        z := start(i,j);
        if ((REAL_PART(z) /= zero) or (IMAG_PART(z) /= zero)) then
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
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),cond));
    sc : Poly_Sys(1..nq); -- Poly_Sys(1..n);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    nv : constant integer32 := integer32(Degree_of_Freedom(locmap));
    x,xt,mfxt : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    flag : constant Standard_Complex_Poly_Matrices.Matrix(1..n,1..n)
         := Cheater_Homotopy_Flag(nv,start,target);
    mf : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Setup_Flag_Homotopies.Moved_Flag(n);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := Standard_Embed_Polynomials.Add_Variables(x,1);
    mfxt := Setup_Flag_Homotopies.Moving_Flag(mf,xt);
    sc := Numeric_Schubert_Conditions.Expand(n,k,nq,cond,mfxt,flag);
    f := new Poly_Sys'(Filter_Zero_Equations(sc));
  end Flag_Conditions;

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               cond : in Bracket;
               start,target : in DoblDobl_Complex_Matrices.Matrix;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),cond));
    sc : Poly_Sys(1..nq); -- used to be n ??
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    nv : constant integer32 := integer32(Degree_of_Freedom(locmap));
    x,xt,mfxt : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    flag : constant DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..n)
         := Cheater_Homotopy_Flag(nv,start,target);
    mf : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
       := Setup_Flag_Homotopies.Moved_Flag(n);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := DoblDobl_Embed_Polynomials.Add_Variables(x,1);
    mfxt := Setup_Flag_Homotopies.Moving_Flag(mf,xt);
    sc := Numeric_Schubert_Conditions.Expand(n,k,nq,cond,mfxt,flag);
    f := new Poly_Sys'(Filter_Zero_Equations(sc));
  end Flag_Conditions;

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               cond : in Bracket;
               start,target : in QuadDobl_Complex_Matrices.Matrix;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    nq : constant integer32
       := integer32(Symbolic_Schubert_Conditions.Number_of_Equations
                      (natural32(n),cond));
    sc : Poly_Sys(1..nq);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    nv : constant integer32 := integer32(Degree_of_Freedom(locmap));
    x,xt,mfxt : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    flag : constant QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..n)
         := Cheater_Homotopy_Flag(nv,start,target);
    mf : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
       := Setup_Flag_Homotopies.Moved_Flag(n);

  begin
    x := Symbolic_Schubert_Conditions.Symbolic_Form_of_Plane(n,k,locmap);
    xt := QuadDobl_Embed_Polynomials.Add_Variables(x,1);
    mfxt := Setup_Flag_Homotopies.Moving_Flag(mf,xt);
    sc := Numeric_Schubert_Conditions.Expand(n,k,nq,cond,mfxt,flag);
    f := new Poly_Sys'(Filter_Zero_Equations(sc));
  end Flag_Conditions;

  procedure Many_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               conds : in Array_of_Brackets;
               start,target : in Standard_Complex_VecMats.VecMat;
               f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    eqs : Array_of_Poly_Sys(conds'range);
    cnt : integer32 := 0;
    ind : integer32;
    current : Link_to_Poly_Sys;

  begin
    for i in conds'range loop
      Flag_Conditions
        (n,k,p,rows,cols,conds(i).all,start(i).all,target(i).all,eqs(i));
      cnt := cnt + eqs(i)'last;
    end loop;
    f := new Poly_Sys(1..cnt);
    ind := 0;
    for i in eqs'range loop       -- copy the pointers to the systems
      current := eqs(i);
      for k in current'range loop -- flattening the double indexed array
        ind := ind + 1;
        f(ind) := current(k);
      end loop;
    end loop;
  end Many_Flag_Conditions;

  procedure Many_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               conds : in Array_of_Brackets;
               start,target : in DoblDobl_Complex_VecMats.VecMat;
               f : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    eqs : Array_of_Poly_Sys(conds'range);
    cnt : integer32 := 0;
    ind : integer32;
    current : Link_to_Poly_Sys;

  begin
    for i in conds'range loop
      Flag_Conditions
        (n,k,p,rows,cols,conds(i).all,start(i).all,target(i).all,eqs(i));
      cnt := cnt + eqs(i)'last;
    end loop;
    f := new Poly_Sys(1..cnt);
    ind := 0;
    for i in eqs'range loop       -- copy the pointers to the systems
      current := eqs(i);
      for k in current'range loop -- flattening the double indexed array
        ind := ind + 1;
        f(ind) := current(k);
      end loop;
    end loop;
  end Many_Flag_Conditions;

  procedure Many_Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               conds : in Array_of_Brackets;
               start,target : in QuadDobl_Complex_VecMats.VecMat;
               f : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    eqs : Array_of_Poly_Sys(conds'range);
    cnt : integer32 := 0;
    ind : integer32;
    current : Link_to_Poly_Sys;

  begin
    for i in conds'range loop
      Flag_Conditions
        (n,k,p,rows,cols,conds(i).all,start(i).all,target(i).all,eqs(i));
      cnt := cnt + eqs(i)'last;
    end loop;
    f := new Poly_Sys(1..cnt);
    ind := 0;
    for i in eqs'range loop       -- copy the pointers to the systems
      current := eqs(i);
      for k in current'range loop -- flattening the double indexed array
        ind := ind + 1;
        f(ind) := current(k);
      end loop;
    end loop;
  end Many_Flag_Conditions;

end Moving_Flag_Homotopies;
