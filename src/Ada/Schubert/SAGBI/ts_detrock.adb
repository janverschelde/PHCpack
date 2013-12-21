with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Osculating_Planes;                  use Osculating_Planes;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Partitions_of_Sets_Of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_Of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with m_Homogeneous_Bezout_Numbers;       use m_Homogeneous_Bezout_Numbers;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Standard_Dynamic32_Triangulations;  use Standard_Dynamic32_Triangulations;
with Matrix_Indeterminates;
with SAGBI_Homotopies;                   use SAGBI_Homotopies;

procedure ts_detrock is

-- DESCRIPTION :
--   Generates (m,p)-system and performs root counting.

  m,p : natural32 := 0;
  ans : character;
  file : file_type;

  function Determinant_System ( m,p : natural32 ) return Poly_Sys is

    res : Poly_Sys(1..integer32(m*p));
    lp : Poly;
    s : double_float := Random;
    inc : constant double_float := 2.0/double_float(m*p);
    mat : Matrix(1..integer32(m+p),1..integer32(m));

  begin
    Matrix_Indeterminates.Initialize_Symbols(m+p,p);
    lp := Lifted_Localized_Laplace_Expansion(m+p,p);
    for i in res'range loop
      s := s+inc;
      mat := Orthogonal_Basis(m+p,m,s);
      res(i) := Intersection_Condition(mat,lp);
      if s > 2.0
       then s := s - 2.0;
      end if;
    end loop;
    return res;
  end Determinant_System;

  procedure Count_Roots ( file : in file_type; h : in Poly_Sys;
                          m,p : in natural32; title_banner : in string ) is

    function Minimum ( a,b : natural32 ) return natural32 is
    begin
      if a <= b
       then return a;
       else return b;
      end if;
    end Minimum;

    function Construct_Partition ( m,p : natural32 ) return Partition is

      min_mp : constant natural32 := Minimum(m,p);
      z : Partition(1..min_mp);
      cnt : natural32 := 0;

    begin
      for i in z'range loop
        z(i) := Create(m*p);
      end loop;
      if m <= p then
        for i in 1..m loop
          for j in 1..p loop
            cnt := cnt+1;
            Add(z(i),cnt);
          end loop;
        end loop;
      else
        cnt := 1;
        for i in 1..p loop
          for j in 1..m loop
            Add(z(i),cnt);
            cnt := cnt+p;
            if cnt > m*p
             then cnt := cnt-m*p;
            end if;
          end loop;
        end loop;
      end if;
      return z;
    end Construct_Partition;

    procedure Multi_Homogeneous_Bound ( f : in Poly_Sys ) is

      nz : natural32;
      b : natural64;
     -- z : Partition(p'range);
      min_mp : constant natural32 := Minimum(m,p);
      z : Partition(1..min_mp) := Construct_Partition(m,p);

    begin
     -- PB(f,b,nz,z);
      nz := z'last;
      b := Bezout_Number(f,z);
      put(file,nz,1); put(file,"-homogeneous Bezout number : ");
      put(file,b,1); new_line(file);
      put(file,"  with partition "); put(file,z); new_line(file);
      Clear(z);
    end Multi_Homogeneous_Bound;

    procedure Apply_Root_Counts ( f : in Poly_Sys; cmpvol : in boolean ) is

      d : constant natural32 := Total_Degree(f);
      sup,lifted,lifted_last : List;
      t : Triangulation;
      vol : natural32;
    
    begin
      new_line(file);
      put_line(file,"ROOT COUNTS : ");
      new_line(file);
      put(file,"total degree : "); put(file,d,1); new_line(file);
      Multi_Homogeneous_Bound(f);
      if cmpvol then
        sup := Create(f(f'first));
        Dynamic_Lifting(sup,false,false,0,lifted,lifted_last,t);
        vol := Volume(t);
        put(file,"mixed volume : "); put(file,vol,1); new_line(file);
        new_line(file);
        put_line(file,"The lifted support : ");
        new_line(file);
        put(file,lifted);
       -- new_line(file);
       -- put_line(file,"The regular triangulation : ");
       -- new_line(file);
       -- put(file,p'length,t,vol);
        Clear(t); Clear(sup); Clear(lifted);
      end if;
    end Apply_Root_Counts;

    procedure Main is

      target,start : Poly_Sys(h'range);

    begin
      target := Eval(h,Create(1.0),integer32(m*p+1));
      put(file,target'length,target);
      new_line(file);
      put_line(file,title_banner);
      Apply_Root_Counts(target,false);
      start := Eval(h,Create(0.0),integer32(m*p+1));
      new_line(file);
      put(file,start'length,start);
      new_line(file);
      put_line(file,"TITLE : start system in SAGBI homotopy.");
      Apply_Root_Counts(start,true);
      Clear(target); Clear(start);
    end Main;

  begin
    Main;
  end Count_Roots;

begin
  new_line;
  put_line("Performing root counts on determinantal (m,p)-systems.");
  loop
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    put("Give m : "); get(m);
    put("Give p : "); get(p);
    declare
      title : constant string := "TITLE : all " & Convert(integer32(p))
        & "-planes that intersect " & Convert(integer32(m*p))
        & " random real osculating " & Convert(integer32(m))
        & "-planes.";
	begin
   -- put(file,"Determinantal ("); put(file,m,1); put(file,",");
   -- put(file,p,1); put_line(file,")-system :");
      Count_Roots(file,Determinant_System(m,p),m,p,title);
    end;
    Close(file);
    put("Do you want other systems to test ? (y/n) ");
    Ask_Yes_or_No(ans);
    exit when (ans /= 'y');
  end loop;
end ts_detrock;
