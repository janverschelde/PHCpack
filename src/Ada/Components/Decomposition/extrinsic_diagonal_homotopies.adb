-- io added for debugging
--with text_io,integer_io; use text_io,integer_io;

with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Witness_Sets;                      use Witness_Sets;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Standard_Diagonal_Polynomials;
with Standard_Diagonal_Solutions;
with DoblDobl_Diagonal_Polynomials;
with DoblDobl_Diagonal_Solutions;
with QuadDobl_Diagonal_Polynomials;
with QuadDobl_Diagonal_Solutions;

package body Extrinsic_Diagonal_Homotopies is

  function Cascade_Dimension ( n1,n2,a,b : natural32 ) return natural32 is

    k : constant natural32 := n1-a;
    res : natural32;

  begin
    if a + b >= k
     then res := 3*k - a;
     else res := 2*k + b;
    end if;
    return res;
  end Cascade_Dimension;

  function Cascade_Dimension
             ( p1e,p2e : Standard_Complex_Poly_Systems.Poly_Sys;
               a,b : natural32 ) return natural32 is

    use Standard_Complex_Polynomials;

    n1 : constant natural32 := Number_of_Unknowns(p1e(p1e'first));
    n2 : constant natural32 := Number_of_Unknowns(p2e(p2e'first));

  begin
    return Cascade_Dimension(n1,n2,a,b);
  end Cascade_Dimension;

  function Cascade_Dimension
             ( p1e,p2e : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               a,b : natural32 ) return natural32 is

    use DoblDobl_Complex_Polynomials;

    n1 : constant natural32 := Number_of_Unknowns(p1e(p1e'first));
    n2 : constant natural32 := Number_of_Unknowns(p2e(p2e'first));

  begin
    return Cascade_Dimension(n1,n2,a,b);
  end Cascade_Dimension;

  function Cascade_Dimension
             ( p1e,p2e : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               a,b : natural32 ) return natural32 is

    use QuadDobl_Complex_Polynomials;

    n1 : constant natural32 := Number_of_Unknowns(p1e(p1e'first));
    n2 : constant natural32 := Number_of_Unknowns(p2e(p2e'first));

  begin
    return Cascade_Dimension(n1,n2,a,b);
  end Cascade_Dimension;

  function Cascade_Dimension
             ( p1e,p2e : Standard_Complex_Laur_Systems.Laur_Sys;
               a,b : natural32 ) return natural32 is

    use Standard_Complex_Laurentials;

    n1 : constant natural32 := Number_of_Unknowns(p1e(p1e'first));
    n2 : constant natural32 := Number_of_Unknowns(p2e(p2e'first));

  begin
    return Cascade_Dimension(n1,n2,a,b);
  end Cascade_Dimension;

  function Cascade_Dimension
             ( p1e,p2e : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               a,b : natural32 ) return natural32 is

    use DoblDobl_Complex_Laurentials;

    n1 : constant natural32 := Number_of_Unknowns(p1e(p1e'first));
    n2 : constant natural32 := Number_of_Unknowns(p2e(p2e'first));

  begin
    return Cascade_Dimension(n1,n2,a,b);
  end Cascade_Dimension;

  function Cascade_Dimension
             ( p1e,p2e : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               a,b : natural32 ) return natural32 is

    use QuadDobl_Complex_Laurentials;

    n1 : constant natural32 := Number_of_Unknowns(p1e(p1e'first));
    n2 : constant natural32 := Number_of_Unknowns(p2e(p2e'first));

  begin
    return Cascade_Dimension(n1,n2,a,b);
  end Cascade_Dimension;

  procedure Cascade1
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_VecVecs;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Diagonal_Polynomials;

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    cdi : constant Poly_Sys
        := Complete(natural32(2*n1),natural32(2*n1)-a-b,dia);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(integer32(b),rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
    ind_target := ind_start;
    for i in cdi'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(cdi(i),b);
    end loop;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in 1..integer32(b) loop      -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-integer32(b)+i);
    end loop;
    for i in 1..integer32(b) loop      -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(i).all);
    end loop;
  end Cascade1;

  procedure Cascade1
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Diagonal_Polynomials;

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    cdi : constant Poly_Sys
        := Complete(natural32(2*n1),natural32(2*n1)-a-b,dia);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(integer32(b),rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
    ind_target := ind_start;
    for i in cdi'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(cdi(i),b);
    end loop;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in 1..integer32(b) loop      -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-integer32(b)+i);
    end loop;
    for i in 1..integer32(b) loop      -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(i).all);
    end loop;
  end Cascade1;

  procedure Cascade1
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Diagonal_Polynomials;

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    cdi : constant Poly_Sys
        := Complete(natural32(2*n1),natural32(2*n1)-a-b,dia);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(integer32(b),rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
    ind_target := ind_start;
    for i in cdi'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(cdi(i),b);
    end loop;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in 1..integer32(b) loop      -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-integer32(b)+i);
    end loop;
    for i in 1..integer32(b) loop      -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(i).all);
    end loop;
  end Cascade1;

  procedure Cascade1
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_VecVecs;
    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;
    use Standard_Diagonal_Polynomials;

    p1 : constant Laur_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Laur_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    rp1 : constant Laur_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Laur_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Laur_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Laur_Sys := Diagonal(n1);
    cdi : constant Laur_Sys
        := Complete(natural32(2*n1),natural32(2*n1)-a-b,dia);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(integer32(b),rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
    ind_target := ind_start;
    for i in cdi'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(cdi(i),b);
    end loop;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in 1..integer32(b) loop      -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-integer32(b)+i);
    end loop;
    for i in 1..integer32(b) loop      -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(i).all);
    end loop;
  end Cascade1;

  procedure Cascade1
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Diagonal_Polynomials;

    p1 : constant Laur_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Laur_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    rp1 : constant Laur_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Laur_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Laur_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Laur_Sys := Diagonal(n1);
    cdi : constant Laur_Sys
        := Complete(natural32(2*n1),natural32(2*n1)-a-b,dia);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(integer32(b),rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
    ind_target := ind_start;
    for i in cdi'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(cdi(i),b);
    end loop;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in 1..integer32(b) loop      -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-integer32(b)+i);
    end loop;
    for i in 1..integer32(b) loop      -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(i).all);
    end loop;
  end Cascade1;

  procedure Cascade1
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Diagonal_Polynomials;

    p1 : constant Laur_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Laur_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    rp1 : constant Laur_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Laur_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Laur_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Laur_Sys := Diagonal(n1);
    cdi : constant Laur_Sys
        := Complete(natural32(2*n1),natural32(2*n1)-a-b,dia);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(integer32(b),rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
    ind_target := ind_start;
    for i in cdi'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(cdi(i),b);
    end loop;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,b);
        Clear(hp); Clear(rhp);
      end;
    end loop;
    for i in 1..integer32(b) loop      -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-integer32(b)+i);
    end loop;
    for i in 1..integer32(b) loop      -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(i).all);
    end loop;
  end Cascade1;

  procedure Cascade2
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_VecVecs;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Diagonal_Polynomials;

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    codim1 : constant integer32 := n1 - integer32(a);
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
   -- put("in Cascade2, n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1); new_line;
   -- put("  a = "); put(a,1); put("  b = "); put(b,1);
   -- put("  nz1 = "); put(nz1,1); put("  nz2 = "); put(nz2,1); new_line;
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(codim1,rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
   -- put("ind_start = "); put(ind_start,1); new_line;
    ind_target := ind_start;
    for i in dia'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(dia(i),natural32(codim1));
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..codim1 loop            -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-codim1+i);
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..integer32(b)-codim1 loop -- hyperplanes without slack to target
      ind_target := ind_target + 1;
      declare
        tsl : Standard_Complex_Vectors.Vector(0..target'last);
      begin
        tsl := sli(i)(0..target'last);
        for j in 1..codim1 loop
          tsl(target'last-j+1) := Create(0.0);
        end loop;
        target(ind_target) := Hyperplane(tsl);
      end;
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in 1..codim1 loop            -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(integer32(b)-codim1+i).all);
    end loop;
   -- put("ind_target : "); put(ind_target,1);
   -- put("  target'last : "); put(target'last,1); new_line;
   -- put("ind_start : "); put(ind_start,1);
   -- put("  start'last : "); put(start'last,1); new_line;
  end Cascade2;

  procedure Cascade2
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Diagonal_Polynomials;

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    codim1 : constant integer32 := n1 - integer32(a);
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
   -- put("in Cascade2, n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1); new_line;
   -- put("  a = "); put(a,1); put("  b = "); put(b,1);
   -- put("  nz1 = "); put(nz1,1); put("  nz2 = "); put(nz2,1); new_line;
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(codim1,rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
   -- put("ind_start = "); put(ind_start,1); new_line;
    ind_target := ind_start;
    for i in dia'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(dia(i),natural32(codim1));
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..codim1 loop            -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-codim1+i);
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..integer32(b)-codim1 loop -- hyperplanes without slack to target
      ind_target := ind_target + 1;
      declare
        tsl : DoblDobl_Complex_Vectors.Vector(0..target'last);
      begin
        tsl := sli(i)(0..target'last);
        for j in 1..codim1 loop
          tsl(target'last-j+1) := Create(integer(0));
        end loop;
        target(ind_target) := Hyperplane(tsl);
      end;
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in 1..codim1 loop            -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(integer32(b)-codim1+i).all);
    end loop;
   -- put("ind_target : "); put(ind_target,1);
   -- put("  target'last : "); put(target'last,1); new_line;
   -- put("ind_start : "); put(ind_start,1);
   -- put("  start'last : "); put(start'last,1); new_line;
  end Cascade2;

  procedure Cascade2
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Diagonal_Polynomials;

    p1 : constant Poly_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Poly_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    codim1 : constant integer32 := n1 - integer32(a);
    rp1 : constant Poly_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Poly_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Poly_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Poly_Sys := Diagonal(n1);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
   -- put("in Cascade2, n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1); new_line;
   -- put("  a = "); put(a,1); put("  b = "); put(b,1);
   -- put("  nz1 = "); put(nz1,1); put("  nz2 = "); put(nz2,1); new_line;
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(codim1,rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
   -- put("ind_start = "); put(ind_start,1); new_line;
    ind_target := ind_start;
    for i in dia'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(dia(i),natural32(codim1));
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..codim1 loop            -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-codim1+i);
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..integer32(b)-codim1 loop -- hyperplanes without slack to target
      ind_target := ind_target + 1;
      declare
        tsl : QuadDobl_Complex_Vectors.Vector(0..target'last);
      begin
        tsl := sli(i)(0..target'last);
        for j in 1..codim1 loop
          tsl(target'last-j+1) := Create(integer(0));
        end loop;
        target(ind_target) := Hyperplane(tsl);
      end;
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in 1..codim1 loop            -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(integer32(b)-codim1+i).all);
    end loop;
   -- put("ind_target : "); put(ind_target,1);
   -- put("  target'last : "); put(target'last,1); new_line;
   -- put("ind_start : "); put(ind_start,1);
   -- put("  start'last : "); put(start'last,1); new_line;
  end Cascade2;

  procedure Cascade2
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_VecVecs;
    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;
    use Standard_Diagonal_Polynomials;

    p1 : constant Laur_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Laur_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    codim1 : constant integer32 := n1 - integer32(a);
    rp1 : constant Laur_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Laur_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Laur_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Laur_Sys := Diagonal(n1);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
   -- put("in Cascade2, n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1); new_line;
   -- put("  a = "); put(a,1); put("  b = "); put(b,1);
   -- put("  nz1 = "); put(nz1,1); put("  nz2 = "); put(nz2,1); new_line;
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(codim1,rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
   -- put("ind_start = "); put(ind_start,1); new_line;
    ind_target := ind_start;
    for i in dia'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(dia(i),natural32(codim1));
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..codim1 loop            -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-codim1+i);
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..integer32(b)-codim1 loop -- hyperplanes without slack to target
      ind_target := ind_target + 1;
      declare
        tsl : Standard_Complex_Vectors.Vector(0..target'last);
      begin
        tsl := sli(i)(0..target'last);
        for j in 1..codim1 loop
          tsl(target'last-j+1) := Create(0.0);
        end loop;
        target(ind_target) := Hyperplane(tsl);
      end;
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in 1..codim1 loop            -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(integer32(b)-codim1+i).all);
    end loop;
   -- put("ind_target : "); put(ind_target,1);
   -- put("  target'last : "); put(target'last,1); new_line;
   -- put("ind_start : "); put(ind_start,1);
   -- put("  start'last : "); put(start'last,1); new_line;
  end Cascade2;

  procedure Cascade2
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Diagonal_Polynomials;

    p1 : constant Laur_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Laur_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    codim1 : constant integer32 := n1 - integer32(a);
    rp1 : constant Laur_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Laur_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Laur_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Laur_Sys := Diagonal(n1);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
   -- put("in Cascade2, n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1); new_line;
   -- put("  a = "); put(a,1); put("  b = "); put(b,1);
   -- put("  nz1 = "); put(nz1,1); put("  nz2 = "); put(nz2,1); new_line;
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(codim1,rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
   -- put("ind_start = "); put(ind_start,1); new_line;
    ind_target := ind_start;
    for i in dia'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(dia(i),natural32(codim1));
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..codim1 loop            -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-codim1+i);
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..integer32(b)-codim1 loop -- hyperplanes without slack to target
      ind_target := ind_target + 1;
      declare
        tsl : DoblDobl_Complex_Vectors.Vector(0..target'last);
      begin
        tsl := sli(i)(0..target'last);
        for j in 1..codim1 loop
          tsl(target'last-j+1) := Create(integer(0));
        end loop;
        target(ind_target) := Hyperplane(tsl);
      end;
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in 1..codim1 loop            -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(integer32(b)-codim1+i).all);
    end loop;
   -- put("ind_target : "); put(ind_target,1);
   -- put("  target'last : "); put(target'last,1); new_line;
   -- put("ind_start : "); put(ind_start,1);
   -- put("  start'last : "); put(start'last,1); new_line;
  end Cascade2;

  procedure Cascade2
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Diagonal_Polynomials;

    p1 : constant Laur_Sys := Remove_Embedding1(p1e,a);
    p2 : constant Laur_Sys := Remove_Embedding1(p2e,b);
    n1 : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := integer32(Number_of_Unknowns(p2(p2'first)));
    nz1 : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nz2 : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    codim1 : constant integer32 := n1 - integer32(a);
    rp1 : constant Laur_Sys := Complete(natural32(n1),a,p1(1..p1'last-nz1));
    rp2 : constant Laur_Sys := Complete(natural32(n2),b,p2(1..p2'last-nz2));
    rp : constant Laur_Sys := Product(n1,n2,rp1,rp2);
    dia : constant Laur_Sys := Diagonal(n1);
    s1 : constant VecVec := Slices(p1e,a);
    s2 : constant VecVec := Slices(p2e,b);
    sli : constant VecVec := Random_Hyperplanes(b,natural32(target'last));
    ind_start,ind_target : integer32 := 0;

  begin
   -- put("in Cascade2, n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1); new_line;
   -- put("  a = "); put(a,1); put("  b = "); put(b,1);
   -- put("  nz1 = "); put(nz1,1); put("  nz2 = "); put(nz2,1); new_line;
    for i in rp'range loop             -- product of two systems
      ind_start := ind_start + 1;
      start(ind_start) := Append_Variables(codim1,rp(i));
      Copy(start(ind_start),target(ind_start));
    end loop;
   -- put("ind_start = "); put(ind_start,1); new_line;
    ind_target := ind_start;
    for i in dia'range loop            -- fill in the diagonal to target
      ind_target := ind_target + 1 ;
      target(ind_target) := Add_Embedding(dia(i),natural32(codim1));
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in s1'range loop             -- add hyperplanes of s1 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s1(i)(0..n1));
        rhp : Poly := Append_Variables(n2,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in s2'range loop             -- add hyperplanes of s2 to start
      ind_start := ind_start + 1;
      declare
        hp : Poly := Hyperplane(s2(i)(0..n2));
        rhp : Poly := Insert_Variables(n1,hp);
      begin
        start(ind_start) := Add_Embedding(rhp,natural32(codim1));
        Clear(hp); Clear(rhp);
      end;
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..codim1 loop            -- add dummy slacks to start
      ind_start := ind_start + 1;
      start(ind_start) := Create(start'last,start'last-codim1+i);
    end loop;
   -- put("ind_start : "); put(ind_start,1); new_line;
    for i in 1..integer32(b)-codim1 loop -- hyperplanes without slack to target
      ind_target := ind_target + 1;
      declare
        tsl : QuadDobl_Complex_Vectors.Vector(0..target'last);
      begin
        tsl := sli(i)(0..target'last);
        for j in 1..codim1 loop
          tsl(target'last-j+1) := Create(integer(0));
        end loop;
        target(ind_target) := Hyperplane(tsl);
      end;
    end loop;
   -- put("ind_target : "); put(ind_target,1); new_line;
    for i in 1..codim1 loop            -- add random hyperplanes to target
      ind_target := ind_target + 1;
      target(ind_target) := Hyperplane(sli(integer32(b)-codim1+i).all);
    end loop;
   -- put("ind_target : "); put(ind_target,1);
   -- put("  target'last : "); put(target'last,1); new_line;
   -- put("ind_start : "); put(ind_start,1);
   -- put("  start'last : "); put(start'last,1); new_line;
  end Cascade2;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Polynomials;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
 
  begin
    if a+b < k
     then Cascade1(p1e,p2e,a,b,start,target);
     else Cascade2(p1e,p2e,a,b,start,target);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Polynomials;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
 
  begin
    if a+b < k
     then Cascade1(p1e,p2e,a,b,start,target);
     else Cascade2(p1e,p2e,a,b,start,target);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Polynomials;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
 
  begin
    if a+b < k
     then Cascade1(p1e,p2e,a,b,start,target);
     else Cascade2(p1e,p2e,a,b,start,target);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_Laurentials;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
 
  begin
    if a+b < k
     then Cascade1(p1e,p2e,a,b,start,target);
     else Cascade2(p1e,p2e,a,b,start,target);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    use DoblDobl_Complex_Laurentials;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
 
  begin
    if a+b < k
     then Cascade1(p1e,p2e,a,b,start,target);
     else Cascade2(p1e,p2e,a,b,start,target);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    use QuadDobl_Complex_Laurentials;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
 
  begin
    if a+b < k
     then Cascade1(p1e,p2e,a,b,start,target);
     else Cascade2(p1e,p2e,a,b,start,target);
    end if;
  end Extrinsic_Cascade_Homotopy;

  function Extrinsic_Product
               ( a,b : natural32;
                 s1,s2 : Standard_Complex_Solutions.Solution )
               return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Solutions;
    use Standard_Diagonal_Solutions;

    k : constant natural32 := natural32(s1.n);
    s : constant Solution(s1.n+s2.n) := Product(s1,s2);

  begin
    if a+b < k
     then return Add_Embedding(s,b);
     else return Add_Embedding(s,k-a);
    end if;
  end Extrinsic_Product;

  function Extrinsic_Product
               ( a,b : natural32;
                 s1,s2 : DoblDobl_Complex_Solutions.Solution )
               return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Diagonal_Solutions;

    k : constant natural32 := natural32(s1.n);
    s : constant Solution(s1.n+s2.n) := Product(s1,s2);

  begin
    if a+b < k
     then return Add_Embedding(s,b);
     else return Add_Embedding(s,k-a);
    end if;
  end Extrinsic_Product;

  function Extrinsic_Product
               ( a,b : natural32;
                 s1,s2 : QuadDobl_Complex_Solutions.Solution )
               return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Diagonal_Solutions;

    k : constant natural32 := natural32(s1.n);
    s : constant Solution(s1.n+s2.n) := Product(s1,s2);

  begin
    if a+b < k
     then return Add_Embedding(s,b);
     else return Add_Embedding(s,k-a);
    end if;
  end Extrinsic_Product;

  function Extrinsic_Product
               ( a,b,k : natural32;
                 sols1,sols2 : Standard_Complex_Solutions.Solution_List )
               return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    use Standard_Diagonal_Solutions;

    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k
     then return Add_Embedding(sols,b);
     else return Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Product;

  function Extrinsic_Product
               ( a,b,k : natural32;
                 sols1,sols2 : DoblDobl_Complex_Solutions.Solution_List )
               return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Diagonal_Solutions;

    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k
     then return Add_Embedding(sols,b);
     else return Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Product;

  function Extrinsic_Product
               ( a,b,k : natural32;
                 sols1,sols2 : QuadDobl_Complex_Solutions.Solution_List )
               return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Diagonal_Solutions;

    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k
     then return Add_Embedding(sols,b);
     else return Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Product;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                sols1,sols2 : in Standard_Complex_Solutions.Solution_List;
                start,target : out Standard_Complex_Poly_Systems.Poly_Sys;
                esols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;
    use Standard_Diagonal_Solutions;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k then
      Cascade1(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,b);
    else
      Cascade2(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                sols1,sols2 : in DoblDobl_Complex_Solutions.Solution_List;
                start,target : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                esols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Diagonal_Solutions;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k then
      Cascade1(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,b);
    else
      Cascade2(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                a,b : in natural32;
                sols1,sols2 : in QuadDobl_Complex_Solutions.Solution_List;
                start,target : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                esols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Diagonal_Solutions;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k then
      Cascade1(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,b);
    else
      Cascade2(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in Standard_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                sols1,sols2 : in Standard_Complex_Solutions.Solution_List;
                start,target : out Standard_Complex_Laur_Systems.Laur_Sys;
                esols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Solutions;
    use Standard_Diagonal_Solutions;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k then
      Cascade1(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,b);
    else
      Cascade2(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                sols1,sols2 : in DoblDobl_Complex_Solutions.Solution_List;
                start,target : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                esols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Diagonal_Solutions;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k then
      Cascade1(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,b);
    else
      Cascade2(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Cascade_Homotopy;

  procedure Extrinsic_Cascade_Homotopy
              ( p1e,p2e : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                a,b : in natural32;
                sols1,sols2 : in QuadDobl_Complex_Solutions.Solution_List;
                start,target : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                esols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Diagonal_Solutions;

    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-a;
    sols : constant Solution_List := Product(sols1,sols2);
 
  begin
    if a+b < k then
      Cascade1(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,b);
    else
      Cascade2(p1e,p2e,a,b,start,target);
      esols := Add_Embedding(sols,k-a);
    end if;
  end Extrinsic_Cascade_Homotopy;

end Extrinsic_Diagonal_Homotopies;
