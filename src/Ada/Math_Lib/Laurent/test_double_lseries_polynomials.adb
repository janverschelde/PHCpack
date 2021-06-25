with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_VecVecVecs;
with Symbol_Table_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;   use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Double_Laurent_Series;
with Random_Laurent_Series;             use Random_Laurent_Series;
with Double_Linear_Laurent_Solvers;

package body Test_Double_Lseries_Polynomials is

  procedure Test ( dim,nbr,deg,pwr,low,upp : in integer32 ) is

    plead : Standard_Integer_Vectors.Vector(1..nbr);
    pcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    pmons : Standard_Integer_VecVecs.VecVec(1..nbr);
    xlead : Standard_Integer_Vectors.Vector(1..dim);
    xcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    ye : integer32;
    yc : Standard_Complex_Vectors.Vector(0..deg);

  begin
    Make_Random_Polynomial(dim,nbr,deg,pwr,low,upp,plead,pcffs,pmons);
    put_line("A random polynomial with Laurent series coefficients :");
    Write(plead,pcffs,pmons);
    Random_Vector(dim,deg,low,upp,xlead,xcffs);
    put_line("A random vector of Laurent series :");
    Double_Linear_Laurent_Solvers.Write(xlead,xcffs,"x");
    Eval(deg,plead,pcffs,pmons,xlead,xcffs,ye,yc);
    put_line("The result of the evaluation :");
    Double_Laurent_Series.Write(ye,yc);
  end Test;

  procedure Set_Parameters
              ( p : in Link_to_Laur_Sys;
                neq,dim,tdx,deg : out integer32 ) is

    ans : character;

  begin
    neq := p'last;
    dim := integer32(Number_of_Unknowns(p(p'first)));
    put("Read "); put(neq,1); put(" polynomials in "); put(dim,1);
    put(" variables :"); Symbol_Table_io.Write; new_line;
    if neq = dim then
      tdx := 0;
    else
      tdx := tsymbol_Index;
      put("-> index of t : "); put(tdx,1); new_line;
    end if;
    if tdx = 0 then
      deg := 0;
      new_line;
      put("Give the degree of the series : "); get(deg);
    else
      deg := Minimum_Laurent_Series_Degree(p.all,tdx);
      put("=> the proposed degree of the series : "); put(deg,1); new_line;
      put("Raise the proposed degree ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put("Give the degree of the series : "); get(deg);
      end if;
    end if;
  end Set_Parameters;

  procedure Test_Input is

    p : Link_to_Laur_Sys;
    neq,dim,tdx,deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ..."); get(p);
    new_line;
    put_line("-> your system :"); put(p.all);
    new_line;
    Set_Parameters(p,neq,dim,tdx,deg);
    new_line;
    declare
      plead : Standard_Integer_VecVecs.VecVec(p'range);
      pcffs : Standard_Complex_VecVecs.Array_of_VecVecs(p'range);
      pmons : Standard_Integer_VecVecs.Array_of_VecVecs(p'range);
    begin
      if tdx = 0
       then Make_Series_System(p.all,dim,dim,0,deg,plead,pcffs,pmons);
       else Make_Series_System(p.all,dim,dim-1,tdx,deg,plead,pcffs,pmons);
      end if;
      for k in p'range loop
        put("Polynomial "); put(k,1);
        put_line(" with Laurent series coefficients :");
        Write(plead(k).all,pcffs(k),pmons(k).all);
      end loop;
    end;
  end Test_Input;

  procedure Chop ( tdx : in integer32;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ylead : out Standard_Integer_Vectors.Vector;
                   ycffs : in Standard_Complex_VecVecs.Link_to_VecVec ) is
  begin
    for k in 1..(tdx-1) loop
      ylead(k) := xlead(k);
      ycffs(k) := xcffs(k);
    end loop;
    for k in (tdx+1)..xlead'last loop
      ylead(k-1) := xlead(k);
      ycffs(k-1) := xcffs(k);
    end loop;
  end Chop;

  procedure Test_Jacobian_Evaluation
              ( dim,nvr,tdx,deg : integer32; jp : in Jaco_Mat;
                tva : in Table_Vector_Array ) is

    xlead : Standard_Integer_Vectors.Vector(1..dim);
    ylead : Standard_Integer_Vectors.Vector(1..nvr);
    xcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    ycff : Standard_Complex_VecVecs.VecVec(1..nvr);
    ycffs : constant Standard_Complex_VecVecs.Link_to_VecVec
          := new Standard_Complex_VecVecs.VecVec'(ycff);
    cff : Standard_Complex_Vectors.Link_to_Vector;
    x0cff : Standard_Complex_Vectors.Vector(1..dim);
    low,upp : constant integer32 := 0;
    evamat : Standard_Complex_Matrices.Matrix(jp'range(1),jp'range(2));
    Alead : Standard_Integer_Matrices.Matrix(jp'range(1),1..nvr);
    Acffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    neq : constant integer32 := jp'last(1);
    row : Standard_Complex_VecVecs.Link_to_VecVec;
    cffrow : Standard_Complex_Vectors.Link_to_Vector;

  begin
    new_line;
   -- put("Give the lower bound on the leading exponents : "); get(low);
   -- put("Give the upper bound on the leading exponents : "); get(upp);
    Random_Vector(dim,deg,low,upp,xlead,xcffs);
    put("Leading exponents of a random vector :"); put(xlead,1); new_line;
    put("A "); put(dim,1); put_line("-vector of Laurent series :");
    Double_Linear_Laurent_Solvers.Write(xlead,xcffs,"x");
    for k in 1..dim loop
      cff := xcffs(k);
      x0cff(k) := cff(0);
    end loop;
    if tdx /= 0 then -- set value for t at one
      x0cff(tdx) := Standard_Complex_Numbers.Create(1.0);
      for i in 1..dim loop
        cff := xcffs(i);        -- make the series constant
        for k in 1..deg loop
          cff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop; 
      end loop;
    end if;
    put_line("The leading coefficients :"); put_line(x0cff);
    evamat := Eval(jp,x0cff); 
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,neq,1,nvr,0,deg);
    if tdx = 0 then
      Eval(deg,tva,xlead,xcffs,Alead,Acffs);
    else
      Chop(tdx,xlead,xcffs,ylead,ycffs);
      Eval(deg,tva,ylead,ycffs,Alead,Acffs);
    end if;
    if tdx = 0 then
      for i in evamat'range(1) loop
        row := Acffs(i);
        for j in evamat'range(2) loop
          cffrow := row(j);
          put("J("); put(i,1); put(","); put(j,1); put(")    : ");
          put(evamat(i,j)); new_line;
          put("A("); put(i,1); put(","); put(j,1); put(")(0) : ");
          put(cffrow(0)); new_line;
        end loop;
      end loop;
    else
      for i in evamat'range(1) loop -- recall that t equals one 
        row := Acffs(i);
        for j in 1..(tdx-1) loop
          cffrow := row(j);
          put("J("); put(i,1); put(","); put(j,1); put(")    : ");
          put(evamat(i,j)); new_line;
          for k in 1..deg loop
            cffrow(0) := cffrow(0) + cffrow(k);
          end loop;
          put("A("); put(i,1); put(","); put(j,1); put(")(0) : ");
          put(cffrow(0)); new_line;
        end loop;
        for j in (tdx+1)..evamat'last(2) loop
          cffrow := row(j-1);
          put("J("); put(i,1); put(","); put(j,1); put(")    : ");
          put(evamat(i,j)); new_line;
          for k in 1..deg loop
            cffrow(0) := cffrow(0) + cffrow(k);
          end loop;
          put("A("); put(i,1); put(","); put(j,1); put(")(0) : ");
          put(cffrow(0)); new_line;
        end loop;
      end loop;
    end if;
  end Test_Jacobian_Evaluation;

  procedure Test_Jacobian is

    p : Link_to_Laur_Sys;
    neq,dim,tdx,deg,nvr : integer32 := 0;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ..."); get(p);
    new_line;
    put_line("-> your system :"); put(p.all);
    new_line;
    Set_Parameters(p,neq,dim,tdx,deg);
    if tdx = 0
     then nvr := dim;
     else nvr := dim-1;
    end if;
    new_line;
    declare
      tv : constant Table_Vector(neq)
         := Make_Table_Vector(p.all,dim,nvr,tdx,deg);
      jp : constant Jaco_Mat(1..neq,1..dim) := Create(p.all);
      tva : constant Table_Vector_Array(1..neq)
          := Make_Table_Vector_Array(jp,tdx,deg);
    begin
      put_line("The table representation :"); Write(tv);
      put_line("The symbolic Jacobian matrix :");
      for i in 1..neq loop
        for j in 1..dim loop
          put("J("); put(i,1); put(","); put(j,1); put_line(") :");
          put(jp(i,j)); new_line;
        end loop;
      end loop;
      put_line("The table representation of the Jacobian :"); Write(tva);
      Test_Jacobian_Evaluation(dim,nvr,tdx,deg,jp,tva);
    end;
  end Test_Jacobian;

  procedure Main is

    dim,nbr,deg,pwr,low,upp : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Test the Jacobian matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Jacobian;
    else
      new_line;
      put("Generate a random polynomial ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans /= 'y' then
        Test_Input;
      else
        new_line;
        put("Give the number of variables : "); get(dim);
        put("Give the number of monomials : "); get(nbr);
        put("Give the degree of the series : "); get(deg);
        put("Give the largest power of the variables : "); get(pwr);
        put("Give the lower bound on the leading exponents : "); get(low);
        put("Give the upper bound on the leading exponents : "); get(upp);
        Test(dim,nbr,deg,pwr,low,upp);
      end if;
    end if;
  end Main;

end Test_Double_Lseries_Polynomials;
