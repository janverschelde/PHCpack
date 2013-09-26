with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Integer64_Linear_Solvers;
with Symbol_Table;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
--with Standard_Smith_Normal_Form;
with Standard_Exponent_Transformations;  use Standard_Exponent_Transformations;

procedure ts_unimod is

-- DESCRIPTION :
--   Development of unimodular coordinate transformations.

  procedure Unimodular_Coordinate_Transformation
              ( A : in Standard_Integer64_Matrices.Matrix;
                M : out Standard_Integer64_Matrices.Matrix;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Computes the Smith form of the matrix A and returns in M
  --   a unimodular coordinate transformation, if fail is false.

  -- REQUIRED : 
  --   The rows of A contain the generators for the cone of tropisms
  --   and M is a square matrix of A'range(2).

    U : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(1));
   --   := Standard_Smith_Normal_Form.Identity(A'last(1));
    V : Standard_Integer64_Matrices.Matrix(A'range(2),A'range(2));
   --   := Standard_Smith_Normal_Form.Identity(A'last(2));
    D : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2)) := A;
   -- B : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2));
    C : constant Standard_Integer64_Matrices.Matrix(A'range(2),A'range(1))
      := Standard_Integer64_Matrices.Transpose(A);
    p,dM : integer64;
    uv : boolean;

    use Standard_Integer64_Matrices;

  begin
   -- Standard_Smith_Normal_Form.Diagonalize(U,D,V);
   -- put_line("The Smith form of A : "); put(D);
   -- put_line("The matrix U : "); put(U);
   -- put_line("The matrix V : "); put(V);
   -- B := U*A*V;
   -- put_line("Checking U*A*V : "); put(B); 
    Unimodular_Coordinate_Transformation(C,U,D,V,M,p,uv,fail);
    if fail then
      put_line("No rational parameterization possible.");
    else
      if uv 
       then put_line("UV-type of unimodular transformation");
      end if;
      put_line("The unimodular transformation :"); put(M);
      dM := Standard_Integer64_Linear_Solvers.Det(M);
      put("its determinant : "); Standard_Integer_Numbers_io.put(dM,1);
      new_line;
    end if;
  end Unimodular_Coordinate_Transformation;

  procedure Interactive_Input
              ( n,k : out integer32;
                A : out Standard_Integer64_Matrices.Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Prompts the user for n, k, and a k-by-n matrix A,
  --   where n is the ambient dimension and k the codimension,
  --   equal to the number of pretropisms.

  begin
    new_line;
    n := 0; put("Give the ambient dimension : "); get(n);
    k := 0; put("Give the number of rows : "); get(k);
    put("Reading a "); put(k,1); put("-by-"); put(n,1);
    put_line(" matrix...");
    A := new Standard_Integer64_Matrices.Matrix(1..k,1..n);
    get(A.all);
  end Interactive_Input;

  procedure Input_from_File
              ( n,k : out integer32;
                A : out Standard_Integer64_Matrices.Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and reads n, k, and
  --   a k-by-n matrix A from that file, where n is the ambient
  --   dimension and k the codimension.

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of an input file ...");
    Read_Name_and_Open_File(file);
    k := 0; get(file,k);
    n := 0; get(file,n);
    A := new Standard_Integer64_Matrices.Matrix(1..k,1..n);
    get(file,A.all);
    close(file);
  end Input_from_File;

  function Transform_Coordinates
             ( t : Term; M : Standard_Integer64_Matrices.Matrix )
             return Term is

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector'(t.dg'range => 0);
    for i in res.dg'range loop
      for j in M'range(2) loop
        res.dg(i) := res.dg(i) + integer32(M(i,j))*t.dg(j);
      end loop;
    end loop;
    return res;
  end Transform_Coordinates;

  function Transform_Coordinates
             ( p : Poly; M : Standard_Integer64_Matrices.Matrix )
             return Poly is

  -- DESCRIPTION :
  --   Transforms the coordinates of the polynomial p with M: x = y^M.

    res : Poly := Null_Poly;

    procedure Transform_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Transform_Coordinates(t,M);

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Transform_Term;
    procedure Transform_Terms is new Visiting_Iterator(Transform_Term);

  begin
    Transform_Terms(p);
    return res;
  end Transform_Coordinates;

  function Transform_Coordinates
             ( p : Laur_Sys; M : Standard_Integer64_Matrices.Matrix )
             return Laur_Sys is

  -- DESCRIPTION :
  --   Transforms the coordinates of the polynomials in p with M: x = y^M.

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Transform_Coordinates(p(i),M);
    end loop;
    return res;
  end Transform_Coordinates;

  function Check_Transform ( p : Laur_Sys; k : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the first k variables in every polynomial of p
  --   have the same exponents.  Returns false otherwise.

    res : boolean := false;
    min,max : integer32;

  begin
    for i in p'range loop
      for v in 1..k loop
        min := Minimal_Degree(p(i),v);
        max := Maximal_Degree(p(i),v);
        res := (min /= max);
        if res then
          new_line;
          put("At equation "); put(i,1);
          put(" and variable "); put(v,1); 
          put(" min = "); put(min,1);
          put(" max = "); put(max,1); new_line;
        end if;
        exit when res;
      end loop;
    end loop;
    return res;
  end Check_Transform;

  function Project ( t : Term; k : integer32 ) return Term is

  -- DESCRIPTION :
  --   Removes the first k variables of the term t.

  -- REQUIRED : k < t.dg'last.

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-k);
    for i in res.dg'range loop
      res.dg(i) := t.dg(i+k);
    end loop;
    return res;
  end Project;

  function Project ( p : Poly; k : integer32 ) return Poly is

  -- DESCRIPTION :
  --   Removes the first k coordinates of the polynomial p.

    res : Poly := Null_Poly;

    procedure Project_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Project(t,k);

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Project_Term;
    procedure Project_Terms is new Visiting_Iterator(Project_Term);
 
  begin
    Project_Terms(p);
    return res;
  end Project;

  function Project ( p : Laur_Sys; k : integer32 ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Project(p(i),k);
    end loop;
    return res;
  end Project;

  procedure Initialize_Symbol_Table ( n,k : in integer32 ) is

  -- DESCRIPTION :
  --   Initializes the symbol table with the symbols 'y'
  --   followed by the index k up to k+n.

    sb : Symbol_Table.symbol;

  begin
    sb := (sb'range => ' ');
    sb(sb'first) := 'y';
    Symbol_Table.Clear;
    Symbol_Table.Init(natural32(n));
    for i in k..(k+n-1) loop
      declare
        nb : constant string := Characters_and_Numbers.Convert(i);
      begin
        for j in nb'range loop
          sb(j+1) := nb(j);
        end loop;
        Symbol_Table.Add(sb);
      end;
    end loop;
  end Initialize_Symbol_Table;

  procedure Write_Header ( file : in file_type; n,k : in integer32 ) is

  -- DESCRIPTION :
  --   Writes the header polynomial to file, to initialize the symbol
  --   table properly, starting with k and ending at n+k-1.

  begin
    for i in k..(n+k-2) loop
      put(file,"y"); put(file,i,1); put(file," + ");
    end loop;
    put(file,"y"); put(file,n+k-1,1); put(file," - ");
    for i in k..(n+k-2) loop
      put(file,"y"); put(file,i,1); put(file," - ");
    end loop;
    put(file,"y"); put(file,n+k-1,1); put(file," + "); new_line(file);
  end Write_Header;

  procedure Write_Projection_to_File
              ( p : in Laur_Sys; k : in integer32;
                M : in Standard_Integer64_Matrices.Matrix ) is 

    file : file_type;
    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));

  begin
    new_line;
    put_line("Writing the projected system to file ...");
    Read_Name_and_Create_File(file);
    put(file,p'last,1);
    put(file," "); put(file,n,1); new_line(file);
    Write_Header(file,n,k);
    put(file,p);
    new_line(file);
    put_line(file,"The transformation matrix M :");
    put(file,M);
    close(file);
  end Write_Projection_to_File;

  procedure Transform_and_Project
              ( p : in Laur_Sys; k : in integer32;
                M : in Standard_Integer64_Matrices.Matrix;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Performs the coordinate transformation x = y^M to the system p,
  --   checks if the first k variables can be factored out.

    q : Laur_Sys(p'range) := Transform_Coordinates(p,M);
    r : Laur_Sys(p'range);
    ans : character;
    n : integer32;

  begin
    new_line;
    put("checking the transformation ...");
    fail := Check_Transform(q,k);
    if fail
     then put_line(" failed!");
     else put_line(" okay.");
    end if;
    new_line;
    put_line("The transformed system : "); put(q);
    new_line;
    put("Remove the first "); put(k,1); put(" variables ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      r := Project(q,k);
      n := integer32(Number_of_Unknowns(r(r'first)));
      Initialize_Symbol_Table(n,k);
      put_line("After projection : "); put(r);
      Write_Projection_to_File(r,k,M);
    end if;
    Clear(q); Clear(r);
  end Transform_and_Project;

  procedure Transform ( A : in Standard_Integer64_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in A the generators for a cone of tropisms,
  --   computes a unimodular coordinate transformation
  --   to transform a polynomial system.

    M : Standard_Integer64_Matrices.Matrix(A'range(2),A'range(2));
    fail : boolean;
    ans : character;
    lp : Link_to_Laur_Sys;

  begin
    Unimodular_Coordinate_Transformation(A,M,fail);
    if not fail then
      new_line;
      put("Transform a polynomial system ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading a Laurent polynomial system ...");
        get(lp);
        Transform_and_Project(lp.all,A'last(1),M,fail);
      end if;
    end if;
  end Transform;

  procedure Main is

    ans : character;
    n,k : integer32 := 0;
    A : Standard_Integer64_Matrices.Link_to_Matrix;

  begin
    new_line;
    put_line("MENU to enter the input matrix : ");
    put_line("  1. type in the entries of the matrix;");
    put_line("  2. give a file name to read the input matrix.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Interactive_Input(n,k,A);
     else Input_from_File(n,k,A);
    end if;
    put("The "); put(k,1); put("-by-"); put(n,1);
    put_line(" matrix : "); put(A.all);
    Transform(A.all);
  end Main;

begin
  Main;
end ts_unimod;
