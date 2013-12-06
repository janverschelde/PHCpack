with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;     use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;      use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;   use Standard_Complex_Numbers_Polar;
with Multprec_Integer_Numbers_io;      use Multprec_Integer_Numbers_io;
with Standard_Random_Vectors;          use Standard_Random_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Norms_Equals;
with Standard_Floating_Vectors_io;     use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors_io;      use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;    use Standard_Complex_Norms_Equals;
with Standard_Integer_Matrices_io;     use Standard_Integer_Matrices_io;
with Standard_Integer_Linear_Solvers;  use Standard_Integer_Linear_Solvers;
with Multprec_Integer_Linear_Solvers;  use Multprec_Integer_Linear_Solvers;
with Multprec_Integer_Matrices_io;     use Multprec_Integer_Matrices_io;
with Standard_Complex_Exponentiation;  use Standard_Complex_Exponentiation;
with Standard_Complex_Polynomials;     use Standard_Complex_Polynomials;
with Standard_Binomial_Systems;        use Standard_Binomial_Systems; 
with Standard_Radial_Solvers;          use Standard_Radial_Solvers;

package body Standard_Binomial_Solvers is

-- AUXILIARY QUICK TESTER :

  procedure Write_Residuals
             ( file : in file_type;
               A : in Multprec_Integer_Matrices.Matrix;
               c : in Standard_Complex_Vectors.Vector;
               sols : in Solution_List ) is

    tmp : Solution_List := sols;
    y : Standard_Complex_Vectors.Vector(A'range(2));
    AA : Standard_Integer_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        AA(i,j) := Multprec_Integer_Numbers.Create(A(i,j));
      end loop;
    end loop;
    while not Is_Null(tmp) loop
      y := Eval(AA,c,Head_Of(tmp).v);
      put(file,Max_Norm(y)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Residuals;

-- SOLVERS for UPPER TRIANGULAR systems :

  function Rank ( U : Standard_Integer_Matrices.Matrix ) return integer32 is
  begin
    for i in U'range(2) loop
      if i <= U'last(1) then
        if U(i,i) = 0
         then return i-1;
        end if;
      else
        return i;
      end if;
    end loop;
    return U'last(2);
  end Rank;

  function Rank ( U : Multprec_Integer_Matrices.Matrix ) return integer32 is
  begin
    for i in U'range(2) loop
      if i <= U'last(1) then
        if Equal(U(i,i),0)
         then return i-1;
        end if;
      else
        return i;
      end if;
    end loop;
    return U'last(2);
  end Rank;

  function Degree ( U : Standard_Integer_Matrices.Matrix ) return integer32 is

    res : integer32 := 1;

  begin
    for i in U'range(2) loop
      if i <= U'last(1)
       then res := res*U(i,i);
       else return res;
      end if;
    end loop;
    return res;
  end Degree;

  function Degree ( U : Multprec_Integer_Matrices.Matrix )
                  return Integer_Number is

    res : Integer_Number := create(integer(1));

  begin
    for i in U'range(2) loop
      if i <= U'last(1)
       then Mul(res,U(i,i));
       else return res;
      end if;
    end loop;
    return res;
  end Degree;

  function Create ( v : Standard_Complex_Vectors.Vector ) return Solution is

  -- DESCRIPTION :
  --   Creates a solution vector with values in v.

    res : Solution(v'last);

  begin
    res.t := Create(0.0);
    res.m := 1;
    res.v := v;
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Create;

  function Solve_Upper_Square
              ( U : Standard_Integer_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector ) return Solution_List is

    res,res_last : Solution_List;
    accu : Standard_Complex_Vectors.Vector(c'range);

    procedure Enumerate_Roots
                ( v : in out Standard_Complex_Vectors.Vector;
                  k : in integer32 ) is

      eva : Complex_Number;

    begin
      if k > c'last then
        Append(res,res_last,Create(v));
      else
        eva := Create(1.0);
        for j in 1..k-1 loop
          eva := eva*(v(j)**integer(U(j,k)));
        end loop;
        if U(k,k) > 0 then
          for j in 0..U(k,k)-1 loop
            v(k) := Root(c(k)/eva,natural32(U(k,k)),natural32(j));
            Enumerate_Roots(v,k+1);
          end loop;
        else
          for j in 0..(-U(k,k)-1) loop
            v(k) := Root(eva/c(k),natural32(-U(k,k)),natural32(j));
            Enumerate_Roots(v,k+1);
          end loop;
        end if;
      end if;
    end Enumerate_Roots;

  begin
    Enumerate_Roots(accu,1);
    return res;
  end Solve_Upper_Square;

  procedure Solve_Upper_Square
              ( U : in Standard_Integer_Matrices.Matrix;
                c : in Standard_Complex_Vectors.Vector;
                s : in Solution_List ) is

    accu : Link_to_Solution;
    tmp : Solution_List := s;

    procedure Enumerate_Roots ( k : in integer32 ) is

      eva : Complex_Number;

    begin
      if k > c'last then
        tmp := Tail_Of(tmp);
        if not Is_Null(tmp) then          -- move to next location
          for i in accu.v'range loop
            Head_Of(tmp).v(i) := accu.v(i); -- copy current values
          end loop;
          accu := Head_Of(tmp);
        end if;
      else
        eva := Create(1.0);
        for j in 1..k-1 loop
          eva := eva*(accu.v(j)**integer(U(j,k)));
        end loop;
        if U(k,k) > 0 then
          for j in 0..U(k,k)-1 loop
            accu.v(k) := Root(c(k)/eva,natural32(U(k,k)),natural32(j));
            Enumerate_Roots(k+1);
          end loop;
        else
          for j in 0..(-U(k,k)-1) loop
            accu.v(k) := Root(eva/c(k),natural32(-U(k,k)),natural32(j));
            Enumerate_Roots(k+1);
          end loop;
        end if;
      end if;
    end Enumerate_Roots;

  begin
    accu := Head_Of(tmp);
    Enumerate_Roots(1);
  end Solve_Upper_Square;

  function Solve_Upper_Square
              ( U : Multprec_Integer_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector ) return Solution_List is

    res,res_last : Solution_List;
    accu : Standard_Complex_Vectors.Vector(c'range);

    procedure Enumerate_Roots
                ( v : in out Standard_Complex_Vectors.Vector;
                  k : in integer32 ) is

      eva,pwr : Complex_Number;
      diagonal : integer32;

    begin
      if k > c'last then
        Append(res,res_last,Create(v));
      else
        eva := Create(1.0);
        for j in 1..k-1 loop
          pwr := Polar_Exponentiation_ModTwoPi_of_Unit(v(j),U(j,k));
          eva := eva*pwr; -- eva := eva*(v(j)**U(j,k));
        end loop;
        diagonal := Create(U(k,k));
        if diagonal > 0 then
          for j in 0..diagonal-1 loop
            v(k) := Root(c(k)/eva,natural32(diagonal),natural32(j)); 
            Enumerate_Roots(v,k+1);
          end loop;
        else
          for j in 0..(-diagonal-1) loop
            v(k) := Root(eva/c(k),natural32(-diagonal),natural32(j)); 
            Enumerate_Roots(v,k+1);
          end loop;
        end if;
      end if;
    end Enumerate_Roots;

  begin
    Enumerate_Roots(accu,1);
    return res;
  end Solve_Upper_Square;

  procedure Write_Vector_and_Modulus
              ( file : in file_type;
                x : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the vector x along with the modulus of all components.

  begin
    for i in x'range loop
      put(file,x(i));
      put(file,"abs : "); put(file,Radius(x(i))); new_line(file);
    end loop;
  end Write_Vector_and_Modulus;

  function Solve_Upper_Square
              ( file : file_type;
                U : Standard_Integer_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector ) return Solution_List is

    res,res_last : Solution_List;
    accu : Standard_Complex_Vectors.Vector(c'range);

    procedure Enumerate_Roots
                ( v : in out Standard_Complex_Vectors.Vector;
                  k : in integer32 ) is

      eva : Complex_Number;

    begin
      for i in 1..k loop
        put(file,"  "); 
      end loop;
      put(file,"At stage "); put(file,k,1);
      if k > c'last then
        put_line(file,"  storing new solution");
        put_line(file,v);
        Append(res,res_last,Create(v));
      else
        eva := Create(1.0);
        for j in 1..k-1 loop
          eva := eva*(v(j)**integer(U(j,k)));
          put(file,"U("); put(file,j,1); put(file,","); put(file,k,1);
          put(file,") = "); put(file,U(j,k),1); put(file,"  ");
          put(file,"  v("); put(file,j,1); put(file,") = ");
          put(file,v(j)); new_line(file);
          put(file,"eva = "); put(file,eva); new_line(file);
        end loop;
        put(file," exponent is "); put(file,U(k,k),1); new_line(file);
        if U(k,k) > 0 then
          for j in 0..U(k,k)-1 loop
            v(k) := Root(c(k)/eva,natural32(U(k,k)),natural32(j));            
            put(file," root #"); put(file,j,1); put(file," : ");
            put(file,v(k)); new_line(file);
            Enumerate_Roots(v,k+1);
          end loop;
        else
          for j in 0..(-U(k,k)-1) loop
            v(k) := Root(eva/c(k),natural32(-U(k,k)),natural32(j));            
            put(file," root #"); put(file,-j,1); put(file," : ");
            put(file,v(k)); new_line(file);
            Enumerate_Roots(v,k+1);
          end loop;
        end if;
      end if;
    exception
      when others => put_line(file,"Exception raised in Enumerate_Roots");
                     put(file,"k = "); put(file,k,1); new_line(file);
                     put(file,"eva = "); put(file,eva); new_line(file);
                     put(file,"v = "); put_line(file,v); raise;
    end Enumerate_Roots;

  begin
    put_line(file,"enumerating the roots...");
    Enumerate_Roots(accu,1);
    return res;
  exception
    when others => put_line(file,"Exception raised in Solve_Upper Square");
                   put_line(file,"The upper triangular matrix U :");
                   put(file,U);
                   put_line(file,"The coefficient vector c ");
                   Write_Vector_and_Modulus(file,c);
                  -- put_line(file,c);
                   raise;
  end Solve_Upper_Square;

  procedure Solve_Upper_Square
              ( file : in file_type;
                U : in Standard_Integer_Matrices.Matrix;
                c : in Standard_Complex_Vectors.Vector;
                s : in Solution_List ) is

    tmp : Solution_List := s;
    accu : Link_to_Solution;

    procedure Enumerate_Roots ( k : in integer32 ) is

      eva : Complex_Number;

    begin
      for i in 1..k loop
        put(file,"  "); 
      end loop;
      put(file,"At stage "); put(file,k,1);
      if k > c'last then
        put_line(file,"  storing new solution");
        put_line(file,accu.v);
        tmp := Tail_Of(tmp);
        if not Is_Null(tmp) then          -- move to next location
          for i in accu.v'range loop
            Head_Of(tmp).v(i) := accu.v(i); -- copy current values
          end loop;
          accu := Head_Of(tmp);
        end if;
      else
        eva := Create(1.0);
        for j in 1..k-1 loop
          eva := eva*(accu.v(j)**integer(U(j,k)));
          put(file,"U("); put(file,j,1); put(file,","); put(file,k,1);
          put(file,") = "); put(file,U(j,k),1); put(file,"  ");
          put(file,"  v("); put(file,j,1); put(file,") = ");
          put(file,accu.v(j)); new_line(file);
          put(file,"eva = "); put(file,eva); new_line(file);
        end loop;
        put(file," exponent is "); put(file,U(k,k),1); new_line(file);
        if U(k,k) > 0 then
          for j in 0..U(k,k)-1 loop
            accu.v(k) := Root(c(k)/eva,natural32(U(k,k)),natural32(j));
            put(file," root #"); put(file,j,1); put(file," : ");
            put(file,accu.v(k)); new_line(file);
            Enumerate_Roots(k+1);
          end loop;
        else
          for j in 0..(-U(k,k)-1) loop
            accu.v(k) := Root(eva/c(k),natural32(-U(k,k)),natural32(j)); 
            put(file," root #"); put(file,-j,1); put(file," : ");
            put(file,accu.v(k)); new_line(file);
            Enumerate_Roots(k+1);
          end loop;
        end if;
      end if;
    exception
      when others => put_line(file,"Exception raised in Enumerate_Roots");
                     put(file,"k = "); put(file,k,1); new_line(file);
                     put(file,"eva = "); put(file,eva); new_line(file);
                     put(file,"v = "); put_line(file,accu.v); raise;
    end Enumerate_Roots;

  begin
    put_line(file,"enumerating the roots...");
    accu := Head_Of(tmp);
    Enumerate_Roots(1);
  exception
    when others => put_line(file,"Exception raised in Solve_Upper Square");
                   put_line(file,"The upper triangular matrix U :");
                   put(file,U);
                   put_line(file,"The coefficient vector c ");
                   Write_Vector_and_Modulus(file,c);
                  -- put_line(file,c);
                   raise;
  end Solve_Upper_Square;

  function Solve_Upper_Square
              ( file : file_type;
                U : Multprec_Integer_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector ) return Solution_List is

    res,res_last : Solution_List;
    accu : Standard_Complex_Vectors.Vector(c'range);

    procedure Enumerate_Roots
                ( v : in out Standard_Complex_Vectors.Vector;
                  k : in integer32 ) is

      eva,pwr : Complex_Number;
      diagonal : integer32;

    begin
      for i in 1..k loop
        put(file,"  "); 
      end loop;
      put(file,"At stage "); put(file,k,1);
      if k > c'last then
        put_line(file,"  storing new solution");
        put_line(file,v);
        Append(res,res_last,Create(v));
      else
        eva := Create(1.0);
        for j in 1..k-1 loop
          put(file,"U("); put(file,j,1); put(file,","); put(file,k,1);
          put(file,") = "); put(file,U(j,k),1); put(file,"  ");
          put(file,"  v("); put(file,j,1); put(file,") = ");
          put(file,v(j)); new_line(file);
          put(file,"eva = "); put(file,eva); new_line(file);
          pwr := Polar_Exponentiation_ModTwoPi_of_Unit(v(j),U(j,k));
          eva := eva*pwr; -- eva := eva*(v(j)**U(j,k));
        end loop;
        diagonal := Create(U(k,k));
        if diagonal > 0 then
          for j in 0..diagonal-1 loop
            v(k) := Root(c(k)/eva,natural32(diagonal),natural32(j));            
            put(file," root #"); put(file,j,1); put(file," : ");
            put(file,v(k)); new_line(file);
            Enumerate_Roots(v,k+1);
          end loop;
        else
          for j in 0..(-diagonal-1) loop
            v(k) := Root(eva/c(k),natural32(-diagonal),natural32(j));
            put(file," root #"); put(file,-j,1); put(file," : ");
            put(file,v(k)); new_line(file);
            Enumerate_Roots(v,k+1);
          end loop;
        end if;
      end if;
    end Enumerate_Roots;

  begin
    put_line(file,"enumerating the roots...");
    Enumerate_Roots(accu,1);
    return res;
  end Solve_Upper_Square;

  function Extend_to_Square
             ( U : Standard_Integer_Matrices.Matrix )
             return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(U'range(1),U'range(1));

  begin
    for j in U'range(2) loop
       for i in U'range(1) loop   -- copy j-th column from U
         res(i,j) := U(i,j);
       end loop;
    end loop;
    for j in U'last(2)+1..U'last(1) loop
      for i in U'range(1) loop   -- add j-th unit vector to U
        res(i,j) := 0;
      end loop;
      res(j,j) := 1;
    end loop;
    return res;
  end Extend_to_Square;

  function Extend_with_Vector
             ( c,v : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'first..c'last+v'last+1-v'first);

  begin
    res(c'range) := c;
    for i in v'range loop
      res(c'last+i+1-v'first) := v(i); 
    end loop;
    return res;
  end Extend_with_Vector;

  procedure Solve_Upper_General
             ( U : in Standard_Integer_Matrices.Matrix;
               c : in Standard_Complex_Vectors.Vector;
               sols : out Solution_List;
               ext_U : out Standard_Integer_Matrices.Matrix; 
               f,ext_c : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := U'last(1) - U'last(2);

  begin
    f := Random_Vector(1,n);   
    ext_U := Extend_to_Square(U);
    ext_c := Extend_with_Vector(c,f);
    sols := Solve_Upper_Square(ext_U,ext_c);
  end Solve_Upper_General;

  procedure Solve_Upper_General
             ( file : in file_type;
               U : in Standard_Integer_Matrices.Matrix;
               c : in Standard_Complex_Vectors.Vector;
               sols : out Solution_List;
               ext_U : out Standard_Integer_Matrices.Matrix;
               f,ext_c : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := U'last(1) - U'last(2);

  begin
    f := Random_Vector(1,n);   
    ext_U := Extend_to_Square(U);
    ext_c := Extend_with_Vector(c,f);
    put_line(file,"The extended upper triangular matrix U :");
    put(file,ext_U);
    put_line(file,"The extended right hand side vector c :");
    put_line(file,ext_c);
    sols := Solve_Upper_Square(file,ext_U,ext_c);
  end Solve_Upper_General;

  procedure Solve_Upper_General
             ( U : in Standard_Integer_Matrices.Matrix;
               c,f : in Standard_Complex_Vectors.Vector;
               sols : out Solution_List;
               ext_U : out Standard_Integer_Matrices.Matrix;
               ext_c : out Standard_Complex_Vectors.Vector ) is
  begin
    ext_U := Extend_to_Square(U);
    ext_c := Extend_with_Vector(c,f);
    sols := Solve_Upper_Square(ext_U,ext_c);
  end Solve_Upper_General;

  procedure Solve_Upper_General
             ( file : in file_type;
               U : in Standard_Integer_Matrices.Matrix;
               c,f : in Standard_Complex_Vectors.Vector;
               sols : out Solution_List;
               ext_U : out Standard_Integer_Matrices.Matrix;
               ext_c : out Standard_Complex_Vectors.Vector ) is
  begin
    ext_U := Extend_to_Square(U);
    ext_c := Extend_with_Vector(c,f);
    put_line(file,"The extended upper triangular matrix U :");
    put(file,ext_U);
    put_line(file,"The extended right hand side vector c :");
    put_line(file,ext_c);
    sols := Solve_Upper_Square(file,ext_U,ext_c);
  end Solve_Upper_General;

-- GENERAL SOLVERS :

  procedure Solve ( A : in Standard_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Standard_Integer_Matrices.Matrix;
                    Asols : out Solution_List ) is

    Usols : Solution_List;

  begin
    Solve(A,c,r,M,U,Usols,Asols);
  end Solve;

  procedure Solve ( A : in Multprec_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Asols : out Solution_List ) is

    Usols : Solution_List;

  begin
    Solve(A,c,r,M,U,Usols,Asols);
  end Solve;

  procedure Solve ( A : in Standard_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Standard_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List ) is

    nv : constant integer32 := A'last(1);
    nq : constant integer32 := A'last(2);
    rd : constant Standard_Floating_Vectors.Vector(c'range) := Radii(c);
    ec : constant Standard_Complex_Vectors.Vector(c'range) := Scale(c,rd);
    logrd : constant Standard_Floating_Vectors.Vector(c'range) := Log10(rd);
    logx,e10x : Standard_Floating_Vectors.Vector(c'range);

  begin
    U := A;
    Upper_Triangulate(M,U);
    r := Rank(U);
    if r = nq then
      if nv = nq then
        Usols := Solve_Upper_Square(U,ec);
        logx := Radial_Upper_Solve(U,logrd);
        logx := Multiply(M,logx);
        e10x := Exp10(logx);
      else
        declare
          ext_U : Standard_Integer_Matrices.Matrix(U'range(1),U'range(1));
          n : constant integer32 := U'last(1)-U'last(2);
          f : Standard_Complex_Vectors.Vector(1..n);
          ext_c : Standard_Complex_Vectors.Vector(U'range(1));
        begin
          Solve_Upper_General(U,c,Usols,ext_U,f,ext_c);
        end;
      end if;
      Asols := Eval(M,Usols);
      if nv = nq
       then Multiply(Asols,e10x);
      end if;
    end if;
  end Solve;

  procedure Solve ( A : in Multprec_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List ) is

    nv : constant integer32 := A'last(1);
    nq : constant integer32 := A'last(2);
    rd : constant Standard_Floating_Vectors.Vector(c'range) := Radii(c);
    ec : constant Standard_Complex_Vectors.Vector(c'range) := Scale(c,rd);
    logrd : constant Standard_Floating_Vectors.Vector(c'range) := Log10(rd);
    maxlogrd : double_float;
    tolzero : constant double_float := 1.0E-12;
    logx,e10x : Standard_Floating_Vectors.Vector(c'range);

  begin
    put_line("this is what is called ...");
    U := A;
    Upper_Triangulate(M,U);
    r := Rank(U);
    if r = nq then
      if nv = nq then
       -- Usols := Solve_Upper_Square(standard_output,U,ec);
        Usols := Solve_Upper_Square(U,ec);
       -- put_line("logrd = "); put_line(logrd); new_line;
        maxlogrd := Standard_Floating_Norms_Equals.Max_Norm(logrd);
        put("max norm of logrd : "); put(maxlogrd); new_line;
       -- put_line("setting all entries in logrd to zero ...");
        if maxlogrd > tolzero then
          logx := Radial_Upper_Solve(U,logrd);
          logx := Multiply(M,logx);
          e10x := Exp10(logx);
        end if;
      end if;
      put_line("writing residuals once, for Usols ...");
      Write_Residuals(standard_output,A,ec,Usols);
      Asols := Eval(M,Usols);
      put_line("Writing residuals twice ...");
      Write_Residuals(standard_output,A,ec,Asols);
      if maxlogrd > tolzero and nv = nq
       then Multiply(Asols,e10x);
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    A : in Standard_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Standard_Integer_Matrices.Matrix;
                    Asols : out Solution_List ) is

    Usols : Solution_List;

  begin
    Solve(file,A,c,r,M,U,Usols,Asols);
  end Solve;

  procedure Solve ( file : in file_type;
                    A : in Multprec_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Asols : out Solution_List ) is

    Usols : Solution_List;

  begin
    Solve(file,A,c,r,M,U,Usols,Asols);
    Write_Residuals(standard_output,A,c,Asols);
  end Solve;

  procedure Solve ( file : in file_type;
                    A : in Standard_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Standard_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List ) is

    nv : constant integer32 := A'last(1);
    nq : constant integer32 := A'last(2);
    rd : constant Standard_Floating_Vectors.Vector(c'range) := Radii(c);
    ec : constant Standard_Complex_Vectors.Vector(c'range) := Scale(c,rd);
    logrd : constant Standard_Floating_Vectors.Vector(c'range) := Log10(rd);
    logx,e10x : Standard_Floating_Vectors.Vector(c'range);

    use Standard_Integer_Matrices;

  begin
    put_line(file,"The Hermite Normal Form of A: M*A = U.");
    put_line(file,"-> the matrix A : "); put(file,A);
    U := A;
    Upper_Triangulate(M,U);
    put_line(file,"-> the unimodular transformation M :"); put(file,M);
    put_line(file,"-> the upper triangular matrix U :"); put(file,U);
    put_line(file,"-> the product M*A : "); put(file,M*A);
    r := Rank(U);
    put(file,"the rank of A (and U) : "); put(file,r,1); new_line(file);
    if r = nq then
      put_line(file,"-> the system has full rank.");
      if nv = nq  then -- the test A'last(1) = A'last(2) did not work !
        Usols := Solve_Upper_Square(file,U,ec);
        logx := Radial_Upper_Solve(U,logrd);
        logx := Multiply(M,logx);
        e10x := Exp10(logx);
        put_line(file,"the magnitude of the roots : ");
        put_line(e10x); new_line(file);
      else
        declare
          ext_U : Standard_Integer_Matrices.Matrix(U'range(1),U'range(1));
          n : constant integer32 := U'last(1)-U'last(2);
          f : Standard_Complex_Vectors.Vector(1..n);
          ext_c : Standard_Complex_Vectors.Vector(U'range(1));
        begin
          Solve_Upper_General(file,U,c,Usols,ext_U,f,ext_c);
        end;
      end if;
      put(file,"Found "); put(file,Length_Of(Usols),1);
      put_line(file," solutions.");
      Asols := Eval(M,Usols);
      if nv = nq
       then Multiply(Asols,e10x);
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    A : in Multprec_Integer_Matrices.Matrix;
                    c : in Standard_Complex_Vectors.Vector; r : out integer32;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List ) is

    nv : constant integer32 := A'last(1);
    nq : constant integer32 := A'last(2);
    rd : constant Standard_Floating_Vectors.Vector(c'range) := Radii(c);
    ec : constant Standard_Complex_Vectors.Vector(c'range) := Scale(c,rd);
    logrd : constant Standard_Floating_Vectors.Vector(c'range) := Log10(rd);
    logx,e10x : Standard_Floating_Vectors.Vector(c'range);

    use Multprec_Integer_Matrices;

  begin
    put_line(file,"The Hermite Normal Form of A: M*A = U.");
    put_line(file,"-> the matrix A : "); put(file,A);
    U := A;
    Upper_Triangulate(M,U);
    put_line(file,"-> the unimodular transformation M :"); put(file,M);
    put_line(file,"-> the upper triangular matrix U :"); put(file,U);
    put_line(file,"-> the product M*A : "); put(file,M*A);
    r := Rank(U);
    put(file,"the rank of A (and U) : "); put(file,r,1); new_line(file);
    if r = nq then
      put_line(file,"-> the system has full rank.");
      if nv = nq  then -- the test A'last(1) = A'last(2) did not work !
        Usols := Solve_Upper_Square(file,U,ec);
        logx := Radial_Upper_Solve(U,logrd);
        logx := Multiply(M,logx);
        e10x := Exp10(logx);
        put_line(file,"the magnitude of the roots : ");
        put_line(e10x); new_line(file);
      else
        null;
       -- declare
       --   ext_U : Standard_Integer_Matrices.Matrix(U'range(1),U'range(1));
       --   n : constant natural := U'last(1)-U'last(2);
       --   f : Standard_Complex_Vectors.Vector(1..n);
       --   ext_c : Standard_Complex_Vectors.Vector(U'range(1));
       -- begin
       --   Solve_Upper_General(file,U,c,Usols,ext_U,f,ext_c);
       -- end;
      end if;
      put(file,"Found "); put(file,Length_Of(Usols),1);
      put_line(file," solutions.");
      Asols := Eval(M,Usols);
      if nv = nq
       then Multiply(Asols,e10x);
      end if;
    end if;
  end Solve;

-- RESIDUAL CALCULATION :

  function Sum_Residuals
             ( A : Standard_Integer_Matrices.Matrix;
               c : Standard_Complex_Vectors.Vector; sols : Solution_List )
             return double_float is

    res : double_float := 0.0;
    tmp : Solution_List := sols;
    y : Standard_Complex_Vectors.Vector(A'range(2));

  begin
    while not Is_Null(tmp) loop
      y := Eval(A,c,Head_Of(tmp).v);
      res := res + Max_Norm(y);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Sum_Residuals;

  procedure Write_Residuals
             ( file : in file_type;
               A : in Standard_Integer_Matrices.Matrix;
               c : in Standard_Complex_Vectors.Vector;
               sols : in Solution_List ) is

    tmp : Solution_List := sols;
    y : Standard_Complex_Vectors.Vector(A'range(2));

  begin
    while not Is_Null(tmp) loop
      y := Eval(A,c,Head_Of(tmp).v);
      put(file,Max_Norm(y)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Residuals;

end Standard_Binomial_Solvers;
