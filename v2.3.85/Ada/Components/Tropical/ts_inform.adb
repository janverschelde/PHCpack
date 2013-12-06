with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Symbol_Table,Symbol_Table_io;
--with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
--with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
--with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Polynomial_Drops;
with Standard_Power_Transformations;     use Standard_Power_Transformations;
with Standard_Initial_Forms;             use Standard_Initial_Forms;

procedure ts_inform is

-- DESCRIPTION :
--   Computes the initial form of a polynomial system.

  function Number_of_Variables ( p : Laur_Sys ) return natural32 is
  begin
    return Number_of_Unknowns(p(p'first));
  end Number_of_Variables;

  procedure Read_Tropism_Coordinates ( v : in out Vector ) is

  -- DESCRIPTION :
  --   Uses the symbol table to prompt the user for the coordinates
  --   of a tropism returned in v.

  begin
    for i in v'range loop
      Symbol_Table_io.put(Symbol_Table.get(natural32(i)));
      put(" : "); get(v(i));
    end loop;
  end Read_Tropism_Coordinates;

  procedure Input_from_File
              ( n,k : out integer32;
                v : out Standard_Integer_VecVecs.Link_to_VecVec ) is

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
    declare
      A : Standard_Integer_Matrices.Matrix(1..k,1..n);
    begin
      get(file,A);
      v := new Standard_Integer_VecVecs.VecVec(1..k);
      for i in 1..k loop
        v(i) := new Standard_Integer_Vectors.Vector(1..n);
        for j in 1..n loop
          v(i)(j) := A(i,j);
        end loop;
      end loop;
    end;
    close(file);
  end Input_from_File;

  procedure Write_to_File ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name if the system p
  --   has to be written to file.

    ans : character;
    file : file_type;

  begin
    put("write to file ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("reading name of output file ...");
      Read_Name_and_Create_File(file);
      put_line(file,p);
      close(file);
    end if;
  end Write_to_File;

  procedure Read_Cone_of_Tropisms
              ( n : in integer32;
                v : in out Standard_Integer_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Uses the symbol table to prompt the user for the coordinates
  --   for as many n-vectors as the range of v.

  begin
    for i in v'range loop
      put("reading tropism "); put(i,1); put_line(" ...");
      declare
        w : Standard_Integer_Vectors.Vector(1..n);
      begin
        Read_Tropism_Coordinates(w);
        v(i) := new Standard_Integer_Vectors.Vector'(w);
      end;
    end loop;
  end Read_Cone_of_Tropisms;

  function Special_Eliminate ( v : Vector; p : integer32 ) return Matrix is

  -- DESCRIPTION :
  --   Returns a special elimination matrix, assuming v(p) = 1,
  --   where p is the pivot in v.

    res : Matrix(v'range,v'range) := Identity_Matrix(natural32(v'last));

  begin
    for i in v'range loop
      res(p,i) := v(i);
    end loop;
    return res;
  end Special_Eliminate;

  procedure Transform ( q : in Laur_Sys; v : in Vector ) is

  -- DESCRIPTION :
  --   Transforms an initial form system with respect to v.

    p : constant integer32 := Pivot(v);
    T : Matrix(v'range,v'range);
    tq : Laur_Sys(q'range);
    ans : character;
 
  begin
    put("special elimination matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then T := Special_Eliminate(v,p);
     else T := Eliminate(v,p);
    end if;
    put_line("An eliminating power transformation :"); put(T);
    tq := Transform(q,T);
    put_line("The transformed polynomial system :"); put(tq);
    Write_to_File(tq);
    new_line;
    put("Drop the first variable (y/n) ? ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        work : constant Laur_Sys(tq'range)
              := Polynomial_Drops.Remove_Variable(tq,1);
      begin
        put_line(work);
        Write_to_File(work);
      end;
    end if;
  end Transform;

  function Monomial_Count ( p : Laur_Sys ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there are at least two monomials in every 
  --   polynomial of p.

    res : boolean := true;
    nt : natural32;

  begin
    for i in p'range loop
      nt := Number_of_Terms(p(i));
      if nt < 2
       then res := false;
      end if;
      put(" "); put(nt,1);
    end loop;
    new_line;
    return res;
  end Monomial_Count;

  procedure Curve_Tropism ( p : in Laur_Sys ) is

    n : constant natural32 := Number_of_Variables(p);
    v : Vector(1..integer32(n));
    invp : Laur_Sys(p'range);

  begin
    put("Give "); put(n,1); put(" coordinates of tropism : ");
    new_line; Read_Tropism_Coordinates(v);
    invp := Initial(p,v);
    put_line("The initial form system :"); put(invp);
    Transform(invp,v);
  end Curve_Tropism;

  procedure Surface_Tropism ( p : in Laur_Sys ) is

    n,d : integer32 := 0;
    ans : character;
    v : Standard_Integer_VecVecs.Link_to_VecVec;
    work,invp : Laur_Sys(p'range);

  begin
    new_line;
    put("Input on file ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Input_from_File(n,d,v);
    else
      new_line;
      put("What is the dimension of the surface ? ");
      get(d);
      v := new Standard_Integer_VecVecs.VecVec(1..d);
      n := integer32(Number_of_Variables(p));
      Read_Cone_of_Tropisms(n,v.all);
    end if;
    Copy(p,work);
    for i in 1..d loop
      put("The initial form for "); put(v(i).all); put_line(" : ");
      invp := Initial(work,v(i).all);
      put_line(invp);
      Copy(invp,work); Clear(invp);
    end loop;
    Write_to_File(work);
  end Surface_Tropism;

  procedure Fish_for_Tropisms ( p : in Laur_Sys ) is

    n : constant natural32 := Number_of_Variables(p);
    v : Vector(1..integer32(n));
    work,invp : Laur_Sys(p'range);
    ans : character;
    ispretrop : boolean;

  begin
    Copy(p,work);
    loop
      put("Give "); put(n,1); put(" coordinates of tropism : ");
      new_line; Read_Tropism_Coordinates(v);
      invp := Initial(work,v);
      put_line("The initial form system :"); put(invp);
      ispretrop := Monomial_Count(invp);
      if ispretrop then
        put(v); put_line(" is pretropism, continue");
        Copy(invp,work);
      else
        put(v); put(" is not a pretropism, continue ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Clear(invp);
    end loop;
  end Fish_for_Tropisms;

  procedure Main is

    lp : Link_to_Laur_Sys;
    ans : character;

  begin
    new_line;
    put_line("MENU to compute initial forms of a Laurent system :");
    put_line("  1. give one vector as pretropism for a curve;");
    put_line("  2. give a matrix of pretropisms for a surface;");
    put_line("  3. interactively fishing for pretropisms.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    put_line("Reading a polynomial system..."); get(lp);
    if ans = '1' then
      Curve_Tropism(lp.all);
    elsif ans = '2' then
      Surface_Tropism(lp.all);
    else
      Fish_for_Tropisms(lp.all);
    end if;
  end Main;

begin
  Main;
end ts_inform;
