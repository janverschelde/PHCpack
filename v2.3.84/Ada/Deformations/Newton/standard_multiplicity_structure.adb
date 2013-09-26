with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Complex_QR_Least_Squares; use Standard_Complex_QR_Least_Squares;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Numerical_Rank;           use Standard_Numerical_Rank;
with Standard_Nullity_Matrices;         use Standard_Nullity_Matrices;
 
package body Standard_Multiplicity_Structure is

  procedure Orthogonal_Basis
              ( nr,nc : in natural32; a : in Matrix; b : out Matrix ) is

    qraux : Standard_Complex_Vectors.Vector(1..integer32(nc))
          := (1..integer32(nc) => Create(0.0));
    jpvt : Standard_Integer_Vectors.Vector(1..integer32(nc))
         := (1..integer32(nc) => 0);

  begin
    b := a;
    QRD(b,qraux,jpvt,false);
  end Orthogonal_Basis;

  procedure Multiplicity_Structure
             ( f : in Poly_Sys; z : in Standard_Complex_Vectors.Vector;
               tol : in double_float; max_h : in natural32;
               h : out Standard_Natural_Vectors.Vector; m : out natural32 ) is

    nq : constant natural32 := natural32(f'last);
    nv : constant natural32 := natural32(z'last);
    nr,nc : natural32;
    nullity : Standard_Natural_Vectors.Vector(0..integer32(max_h));

  begin
    h := (0..integer32(max_h) => 0);
    h(0) := 1;
    nullity(0) := 0;          -- makes a nice exit in loop below
    for i in 1..h'last loop
      Dimensions_of_Nullity_Matrix(nq,nv,natural32(i),nr,nc);
      declare
        eva : constant Matrix(1..integer32(nr),1..integer32(nc))
            := Evaluate_Nullity_Matrix(nq,nv,nr,nc,natural32(i),f,z);
        rnk : constant natural32 := Numerical_Rank(eva,tol);
      begin
        nullity(i) := nc - rnk;
      end;
      exit when (nullity(i) = nullity(i-1));
      if i = 1
       then h(i) := nullity(i) - 1;
       else h(i) := nullity(i) - nullity(i-1);
      end if;
    end loop;
    m := Standard_Natural_Vectors.Sum(h);
  end Multiplicity_Structure;

  procedure Multiplicity_Structure
             ( file : in file_type;
               f : in Poly_Sys; z : in Standard_Complex_Vectors.Vector;
               tol : in double_float; max_h : in natural32;
               h : out Standard_Natural_Vectors.Vector; m : out natural32 ) is

    nq : constant natural32 := natural32(f'last);
    nv : constant natural32 := natural32(z'last);
    nr,nc : natural32;
    nullity : Standard_Natural_Vectors.Vector(0..integer32(max_h));

  begin
    h := (0..integer32(max_h) => 0);
    h(0) := 1;
    nullity(0) := 0;          -- makes a nice exit in loop below
    for i in 1..h'last loop
      Dimensions_of_Nullity_Matrix(nq,nv,natural32(i),nr,nc);
      declare
        eva : constant Matrix(1..integer32(nr),1..integer32(nc))
            := Evaluate_Nullity_Matrix(file,nq,nv,nr,nc,natural32(i),f,z);
        rnk : constant natural32 := Numerical_Rank(eva,tol);
      begin
        nullity(i) := nc - rnk;
        put(file,"nullity("); put(file,i,1); put(file,") = ");
        put(file,nullity(i),1);
      end;
      exit when (nullity(i) = nullity(i-1));
      if i = 1
       then h(i) := nullity(i) - 1;
       else h(i) := nullity(i) - nullity(i-1);
      end if;
      put(file,"  h("); put(file,i,1); put(file,") = ");
      put(file,h(i),1); new_line(file);
    end loop;
    new_line(file);
    m := Standard_Natural_Vectors.Sum(h);
  end Multiplicity_Structure;

  procedure Multiplicity_Structure
             ( file : in file_type; output : in boolean;
               f : in Poly_Sys; sols : in out Solution_List;
               tol : in double_float; max_h : in natural32 ) is

    nv : constant natural32 := Number_of_Unknowns(f(f'first));
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    new_line(file);
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        z : constant Standard_Complex_Vectors.Vector(1..integer32(nv))
          := ls.v;
        h : Standard_Natural_Vectors.Vector(0..integer32(max_h));
        m : natural32;
      begin
        if output
         then Multiplicity_Structure(file,f,z,tol,max_h,h,m);
         else Multiplicity_Structure(f,z,tol,max_h,h,m);
        end if;
        put(file,"Multiplicity Structure of solution ");
        put(file,i,1); put_line(file," :");
        put(file,"  values of the Hilbert function : ");
        put(file,h); new_line(file);
        put(file,"  the multiplicity : "); put(file,m,1); new_line(file);
        ls.m := integer32(m);
      end;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Multiplicity_Structure;

  procedure Driver_to_Multiplicity_Structure
             ( file : in file_type;
               f : in Poly_Sys; sols : in out Solution_List ) is

    max_h : natural32 := 3;
    tol : double_float := 1.0E-8;
    ans : character;
    first,output : boolean;

  begin
    new_line;
    first := true;
    loop
      if first then
        put("Default value for maximal differential order is ");
        first := false;
      else
        put("Current value for maximal differential order is ");
      end if;
      put(max_h,1); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans = 'n';
      put("Give new value for the maximal differential order : ");
      get(max_h);
    end loop;
    first := true;
    loop
      if first then
        put("Default tolerance for numerical rank is ");
        first := false;
      else
        put("Current tolerance for numerical rank is ");
      end if;
      put(tol,3); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans = 'n';
      put("Give new value for tolerance : "); get(tol);
    end loop;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    Multiplicity_Structure(file,output,f,sols,tol,max_h);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Driver_to_Multiplicity_Structure;

end Standard_Multiplicity_Structure;
