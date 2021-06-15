with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Integer_Linear_Solvers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Solutions;
with Supports_of_Polynomial_Systems;
with Standard_Binomial_Systems;
with Standard_Radial_Solvers;
with Standard_Binomial_Solvers;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;

procedure ts_cspsol is

-- DESCRIPTION :
--   Test on solving initial cell systems given by supports
--   (the points in a mixed cell) and coefficient vectors.

  procedure Make_Binomial_System
              ( s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                c : in Standard_Complex_VecVecs.VecVec;
                A : out Standard_Integer_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Sets up the exponent matrix A and right-hand-side vector b,
  --   given the supports in s and corresponding coefficients c
  --   of a fully mixed binomial system: x^A = b.

  -- REQUIRED : s'range = c'range = A'range(1) = A'range(2) = b'range,
  --   and all coefficients in c are nonzero.

    first,second : Standard_Floating_Vectors.Link_to_Vector;
    use Standard_Complex_Numbers;
    use Lists_of_Floating_Vectors;

  begin
    for i in s'range loop
      first := Head_Of(s(i));
      second := Head_Of(Tail_Of(s(i)));
      for j in A'range(2) loop
        A(i,j) := integer32(first(j)) - integer32(second(j));
      end loop;
      b(i) := (-c(i)(1))/c(i)(2);
    end loop;
  end Make_Binomial_System;

  procedure Solve_Binomial_System
               ( A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector; r : out integer32;
                 M,U : out Standard_Integer_Matrices.Matrix;
                 Usols : out Standard_Complex_Solutions.Solution_List;
                 Asols : out Standard_Complex_Solutions.Solution_List ) is

    nv : constant integer32 := A'last(1);
    nq : constant integer32 := A'last(2);
    rd : constant Standard_Floating_Vectors.Vector(c'range)
       := Standard_Radial_Solvers.Radii(c);
    ec : constant Standard_Complex_Vectors.Vector(c'range)
       := Standard_Radial_Solvers.Scale(c,rd);
    logrd : constant Standard_Floating_Vectors.Vector(c'range)
          := Standard_Radial_Solvers.Log10(rd);
    logx,e10x : Standard_Floating_Vectors.Vector(c'range);

  begin
    U := A;
    Standard_Integer_Linear_Solvers.Upper_Triangulate(M,U);
    r := Standard_Binomial_Solvers.Rank(U);
    if r = nq then
      if nv = nq then
        Usols := Standard_Binomial_Solvers.Solve_Upper_Square(U,ec);
        logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logrd);
        logx := Standard_Radial_Solvers.Multiply(M,logx);
        e10x := Standard_Radial_Solvers.Exp10(logx);
        Asols := Standard_Binomial_Systems.Eval(M,Usols);
        if nv = nq
         then Standard_Radial_Solvers.Multiply(Asols,e10x);
        end if;
      end if;
    end if;
  end Solve_Binomial_System;

  procedure Select_Subsystems
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mcc : in Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Selects all supported subsystems in the fully mixed case.

    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    dim : constant natural32
        := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    deg : Standard_Complex_Polynomials.Degrees
        := new Standard_Natural_Vectors.Vector(1..integer32(dim));
    cff : Standard_Complex_VecVecs.VecVec(p'range);
    A,M,U : Standard_Integer_Matrices.Matrix(p'range,p'range);
    b : Standard_Complex_Vectors.Vector(p'range);
    Usols,Asols : Standard_Complex_Solutions.Solution_List;
    mv,vol : natural32 := 0;
    r : integer32;
    sum : double_float := 0.0;

  begin
    for i in cff'range loop
      cff(i) := new Standard_Complex_Vectors.Vector(1..2);
    end loop;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Supports_of_Polynomial_Systems.Select_Coefficients
        (p,mic.pts.all,dim,deg,cff);
      Make_Binomial_System(mic.pts.all,cff,A,b);
      Standard_Binomial_Solvers.Solve(A,b,r,M,U,Usols,Asols);
      vol := Standard_Complex_Solutions.Length_Of(Asols);
      sum := sum + Standard_Binomial_Solvers.Sum_Residuals(A,b,Asols);
      mv := mv + vol;
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Complex_Polynomials.Clear(deg);
    put("the mixed volume : "); put(mv,1); new_line;
    put("sum of residuals : "); put(sum); new_line;
  end Select_Subsystems;

  procedure Main is

    rq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    infile : file_type;
    mcc : Mixed_Subdivision;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    n,r : natural32;

  begin
    new_line;
    put_line("Reading a random coefficient start system ...");
    get(rq);
    new_line;
    put_line("Reading mixed cells induced by floating-point lifting ...");
    Read_Name_and_Open_File(infile);
    get(infile,n,r,mix,mcc);
    new_line;
    put("Read "); put(Length_Of(mcc),1); put_line(" mixed cells.");
    Select_Subsystems(rq.all,mcc);
  end Main;

begin
  Main;
end ts_cspsol;
