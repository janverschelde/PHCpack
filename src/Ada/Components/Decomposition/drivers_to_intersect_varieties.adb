with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Witness_Sets;                      use Witness_Sets;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Intrinsic_Diagonal_Continuation;   use Intrinsic_Diagonal_Continuation;
with Extrinsic_Diagonal_Continuation;   use Extrinsic_Diagonal_Continuation;

package body Drivers_to_Intersect_Varieties is

  function Square_of_the_Original
             ( n,d : natural32; ep : Poly_Sys ) return Poly_Sys is

    p : constant Poly_Sys(1..ep'last-integer32(d)) := Remove_Embedding1(ep,d);
    nz : constant integer32 := integer32(Number_of_Zero_Equations(p));

  begin
    return Complete(n,d,p(1..p'last-nz));
  end Square_of_the_Original;

  procedure Intrinsic_Diagonal_Homotopy
              ( file : in file_type; report : in boolean;
                ep1,ep2 : in Poly_Sys; esols1,esols2 : in Solution_List;
                a,b : in natural32 ) is

  -- DESCRIPTION :
  --   Diagonal homotopy, assuming a >= b. 

    n : constant integer32
      := integer32(Number_of_Unknowns(ep1(ep1'first)) - a);
    k1 : constant integer32 := n-integer32(a);
    k2 : constant integer32 := n-integer32(b);
    s1 : constant VecVec := Slices(ep1,a);
    s2 : constant VecVec := Slices(ep2,b);
    eqs1 : constant Matrix(1..integer32(a),0..n) := Equations_to_Matrix(s1,n);
    eqs2 : constant Matrix(1..integer32(b),0..n) := Equations_to_Matrix(s2,n);
    gen1 : constant Matrix(1..n,0..k1) := Generators(eqs1);
    gen2 : constant Matrix(1..n,0..k2) := Generators(eqs2);
    pla1 : constant Matrix(1..n,0..k1) := Orthogonalize(gen1);
    pla2 : constant Matrix(1..n,0..k2) := Orthogonalize(gen2);
    isols1 : constant Solution_List := Project(esols1,pla1);
    isols2 : constant Solution_List := Project(esols2,pla2);
    p1 : constant Poly_Sys := Square_of_the_Original(natural32(n),a,ep1);
    p2 : constant Poly_Sys := Square_of_the_Original(natural32(n),b,ep2);
    sols,esols : Solution_List;
    m : constant integer32 := k1+k2;
    plane : Matrix(1..n,0..m);

  begin
    Intrinsic_Diagonal_Homotopy
      (file,report,a,b,p1,p2,isols1,isols2,pla1,pla2,sols,plane);
    new_line(file);
    put(file,p1'last+p2'last,1); put(file,"  ");
    put(file,n,1); new_line(file);
    for i in p1'range loop
      put(file," "); put(file,p1(i)); new_line(file);
    end loop;
    for i in p2'range loop
      put(file," "); put(file,p2(i)); new_line(file);
    end loop;
    new_line(file);
    if not Is_Null(sols) then
      put_line(file,"THE SOLUTIONS :");
      esols := Expand(sols,plane);
      put(file,Length_Of(esols),natural32(Head_Of(esols).n),esols);
    end if;
  end Intrinsic_Diagonal_Homotopy;

  procedure Intrinsic_Diagonal_Homotopy
              ( file : in file_type; report : in boolean;
                ep1,ep2 : in Poly_Sys; esols1,esols2 : in Solution_List;
                a,b : in natural32; f : out Link_to_Poly_Sys;
                p : out Link_to_Matrix; s : out Solution_List ) is

  -- DESCRIPTION :
  --   Diagonal homotopy, assuming a >= b. 

    n : constant integer32
      := integer32(Number_of_Unknowns(ep1(ep1'first)) - a);
    k1 : constant integer32 := n - integer32(a);
    k2 : constant integer32 := n - integer32(b);
    s1 : constant VecVec := Slices(ep1,a);
    s2 : constant VecVec := Slices(ep2,b);
    eqs1 : constant Matrix(1..integer32(a),0..n) := Equations_to_Matrix(s1,n);
    eqs2 : constant Matrix(1..integer32(b),0..n) := Equations_to_Matrix(s2,n);
    gen1 : constant Matrix(1..n,0..k1) := Generators(eqs1);
    gen2 : constant Matrix(1..n,0..k2) := Generators(eqs2);
    pla1 : constant Matrix(1..n,0..k1) := Orthogonalize(gen1);
    pla2 : constant Matrix(1..n,0..k2) := Orthogonalize(gen2);
    isols1 : constant Solution_List := Project(esols1,pla1);
    isols2 : constant Solution_List := Project(esols2,pla2);
    p1 : constant Poly_Sys := Square_of_the_Original(natural32(n),a,ep1);
    p2 : constant Poly_Sys := Square_of_the_Original(natural32(n),b,ep2);
    sols : Solution_List;
    m : constant integer32 := k1+k2;
    plane : Matrix(1..n,0..m);

  begin
    Intrinsic_Diagonal_Homotopy
      (file,report,a,b,p1,p2,isols1,isols2,pla1,pla2,sols,plane);
    f := new Poly_Sys(1..p1'last+p2'last);
    for i in p1'range loop
      Copy(p1(i),f(i));
    end loop;
    for i in p2'range loop
      Copy(p2(i),f(p1'last+i));
    end loop;
    p := new Matrix'(plane);
    s := sols;
  end Intrinsic_Diagonal_Homotopy;

  procedure Generic_Diagonal_Homotopy
              ( file : in file_type; nefA,nefB,n,a,b : in natural32;
                esols1,esols2 : in Solution_List; p1,p2 : in VecVec ) is

    k1 : constant integer32 := integer32(n-a);
    k2 : constant integer32 := integer32(n-b);
    eqs1 : constant Matrix(1..integer32(a),0..integer32(n))
         := Equations_to_Matrix(p1,integer32(n));
    eqs2 : constant Matrix(1..integer32(b),0..integer32(n))
         := Equations_to_Matrix(p2,integer32(n));
    gen1 : constant Matrix(1..integer32(n),0..k1) := Generators(eqs1);
    gen2 : constant Matrix(1..integer32(n),0..k2) := Generators(eqs2);
    pla1 : constant Matrix(1..integer32(n),0..k1) := Orthogonalize(gen1);
    pla2 : constant Matrix(1..integer32(n),0..k2) := Orthogonalize(gen2);
    isols1 : constant Solution_List := Project(esols1,pla1);
    isols2 : constant Solution_List := Project(esols2,pla2);
    sols : Solution_List;
    m : constant integer32 := k1+k2;
    plane : Matrix(1..integer32(n),0..m);

    procedure Diagonal_Homotopy is new
      Intrinsic_Diagonal_Continuation.Generic_Diagonal_Homotopy(fA,jfA,fB,jfB);

  begin
    Diagonal_Homotopy(file,
       integer32(nefA),integer32(nefB),integer32(n),integer32(a),integer32(b),
       isols1,isols2,pla1,pla2,sols,plane);
  end Generic_Diagonal_Homotopy;

  function Complete
             ( r,c : integer32; x : Vector; v : Matrix ) return Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 1..r, copying from x(1..r) and adding
  --   components of x with index higher than r, multiplied with v entries.

    res : Vector(1..r);

  begin
    for i in 1..r loop
      res(i) := x(i);
      for j in 1..c loop
        res(i) := res(i) + v(i,j)*x(r+j);
      end loop;
    end loop;
    return res;
  end Complete;

  function Complete ( r,c,n : integer32; x,v : Matrix ) return Matrix is

  -- DESCRIPTION :
  --   Returns a matrix with r rows and n columns, copying the first r
  --   rows of x and adding to it entries of x with row indices > r,
  --   multiplied with elements of v.

    res : Matrix(1..r,1..n);

  begin
    for i in 1..r loop
      for k in 1..n loop
        res(i,k) := x(i,k);
        for j in 1..c loop
          res(i,k) := res(i,k) + v(i,j)*x(r+j,k);
        end loop;
      end loop;
    end loop;
    return res;
  end Complete;

  procedure B_Call_Generic_Diagonal_Homotopy
              ( file : in file_type; nef1,n,a,b : in natural32;
                ep2 : in Poly_Sys; esols1,esols2 : in Solution_List;
                s1 : in VecVec ) is

    s2 : constant VecVec := Slices(ep2,b);
    p2 : constant Poly_Sys := Remove_Embedding1(ep2,b);
    nz : constant integer32 := integer32(Number_of_Zero_Equations(p2));
    nef2 : constant integer32 := p2'last-nz;
    f2 : constant Eval_Poly_Sys(1..nef2) := Create(p2(1..nef2));
    jmf2 : constant Jaco_Mat(p2'range,1..integer32(n)) := Create(p2);
    ejf2 : constant Eval_Jaco_Mat(p2'range,1..integer32(n)) := Create(jmf2);

  begin
    if n-b = natural32(nef2) then
      declare
        function fB ( x : Vector ) return Vector is
        begin
         return Eval(f2,x);
        end fB;
        function jfB ( x : Vector ) return Matrix is
        begin
         return Eval(ejf2,x);
        end jfB;
        procedure Diagonal_Homotopy is
          new Generic_Diagonal_Homotopy(fA,jfA,fB,jfB);
      begin
        Diagonal_Homotopy(file,nef1,natural32(nef2),n,a,b,esols1,esols2,s1,s2);
      end;
    else
      declare
        r : constant integer32 := integer32(n-b);
        c : constant integer32 := nef2-r;
        v : constant Matrix(1..r,1..c)
          := Standard_Random_Matrices.Random_Matrix(natural32(r),natural32(c));
        function fB ( x : Vector ) return Vector is
        begin
         return Complete(r,c,Eval(f2,x),v);
        end fB;
        function jfB ( x : Vector ) return Matrix is
        begin
         return Complete(r,c,integer32(n),Eval(ejf2,x),v);
        end jfB;
        procedure Diagonal_Homotopy is
          new Generic_Diagonal_Homotopy(fA,jfA,fB,jfB);
      begin
        Diagonal_Homotopy(file,nef1,natural32(r),n,a,b,esols1,esols2,s1,s2);
      end;
    end if;
  end B_Call_Generic_Diagonal_Homotopy;

  procedure A_Call_Generic_Diagonal_Homotopy
              ( file : in file_type; ep1,ep2 : in Poly_Sys;
                esols1,esols2 : in Solution_List; a,b : in natural32 ) is

    n : constant integer32
      := integer32(Number_of_Unknowns(ep1(ep1'first)) - a);
    s1 : constant VecVec := Slices(ep1,a);
    p1 : constant Poly_Sys := Remove_Embedding1(ep1,a);
    nz : constant integer32 := integer32(Number_of_Zero_Equations(p1));
    nef1 : constant integer32 := p1'last - nz;
    f1 : constant Eval_Poly_Sys(1..nef1) := Create(p1(1..nef1));
    jmf1 : constant Jaco_Mat(p1'range,1..n) := Create(p1);
    ejf1 : constant Eval_Jaco_Mat(p1'range,1..n) := Create(jmf1);

  begin
    if n-integer32(a) = nef1 then
      declare
        function fA ( x : Vector ) return Vector is
        begin
          return Eval(f1,x);
        end fA;
        function jfA ( x : Vector ) return Matrix is
        begin
          return Eval(ejf1,x);
        end jfA;
        procedure A_Diagonal_Homotopy is
          new B_Call_Generic_Diagonal_Homotopy(fA,jfA);
      begin
        A_Diagonal_Homotopy(file,natural32(nef1),natural32(n),a,b,
                            ep2,esols1,esols2,s1);
      end;
    else
      declare
        r : constant integer32 := n - integer32(a);
        c : constant integer32 := nef1-r;
        v : constant Matrix(1..r,1..c)
          := Standard_Random_Matrices.Random_Matrix(natural32(r),natural32(c));
        function fA ( x : Vector ) return Vector is
        begin
          return Complete(r,c,Eval(f1,x),v);
        end fA;
        function jfA ( x : Vector ) return Matrix is
        begin
          return Complete(r,c,n,Eval(ejf1,x),v);
        end jfA;
        procedure A_Diagonal_Homotopy is
          new B_Call_Generic_Diagonal_Homotopy(fA,jfA);
      begin
        A_Diagonal_Homotopy(file,natural32(r),natural32(n),a,b,
                            ep2,esols1,esols2,s1);
      end;
    end if;
  end A_Call_Generic_Diagonal_Homotopy;

end Drivers_to_Intersect_Varieties;
