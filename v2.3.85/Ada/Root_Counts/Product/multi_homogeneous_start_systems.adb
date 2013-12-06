with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Degree_Structure;                   use Degree_Structure;
with Standard_Linear_Product_System;
with Standard_Complex_Prod_Planes;

package body Multi_Homogeneous_Start_Systems is

  procedure Add_Hyperplanes ( i : in natural32; p : in Poly;
                              z : in partition; m : in natural32;
                              dg : in Standard_Natural_Vectors.Vector) is
  -- DESCRIPTION : 
  --   This routine adds hyperplanes to the Random Product System 
  --   according to the partition and the degrees of the polynomial 

  -- ON ENTRY :
  --   i        the number of the polynomial in the system;
  --   p        the i-the polynomial of the system;
  --   z        a partition of the set of unknowns of p;
  --   m        the number of sets in z;
  --   dg       the degrees of the polynomial in the sets of z,
  --            for j in 1..m : dg(j) = Degree(p,z(j)).

    n : constant natural32 := Number_Of_Unknowns(p);
    h : Standard_Complex_Vectors.Vector(0..integer32(n));

  begin
    for k in 1..m loop
      for l in 1..dg(integer32(k)) loop
        h(0) := random1;
        for j in 1..n loop
          if Is_In(z(k),j)
           then h(integer32(j)) := Random1;
           else h(integer32(j)) := Create(0.0);
          end if;
        end loop;
        Standard_Linear_Product_System.Add_Hyperplane(i,h);
      end loop;
    end loop;
  end Add_Hyperplanes;

  procedure Construct_Random_Product_System ( p : in Poly_Sys ) is

  -- DESCRIPTION : 
  --   This procedure constructs a random product system by
  --   finding a good partition for each equation of the system p.

    m : natural32;

  begin
    if Degree_Structure.Empty
     then Degree_Structure.Create(p);
    end if;
    for i in p'range loop
      m := Degree_Structure.Get(natural32(i));
      declare
        z : Partition(1..m);
        dg : Standard_Natural_Vectors.Vector(1..integer32(m));
      begin
        Degree_Structure.Get(natural32(i),z,dg);
        Add_Hyperplanes(natural32(i),p(i),z,m,dg);
        Clear(z);
      end;
    end loop;
  end Construct_Random_Product_System;

  procedure RPQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List; nl : out natural32 ) is

    n : constant natural32 := natural32(p'length);

  begin
    Standard_Linear_Product_System.Init(n);
    Construct_Random_Product_System(p);
    q := Standard_Linear_Product_System.Polynomial_System;
    Standard_Linear_Product_System.Solve(sols,nl);
    Degree_Structure.Clear;
    Standard_Linear_Product_System.Clear;
  end RPQ;

  procedure Solve ( i,n : in integer32; first,last : in out Solution_List;
                    a : in out Matrix;
                    b : in out Standard_Complex_Vectors.Vector;
                    acc : in out partition ) is

  begin
    if i > n then
      declare
        aa : Matrix(a'range(1),a'range(2));
        bb : Standard_Complex_Vectors.Vector(b'range);
        rcond : double_float;
        ipvt : Standard_Integer_Vectors.Vector(b'range);
      begin
        for k in a'range(1) loop
          for l in a'range(2) loop
            aa(k,l) := a(k,l);
          end loop;
          bb(k) := b(k);
        end loop;
        lufco(aa,n,ipvt,rcond);
       -- put("rcond : "); put(rcond); new_line;
        if rcond + 1.0 /= 1.0 then
          lusolve(aa,n,ipvt,bb);
          declare
            s : Solution(n);
          begin
            s.t := Create(0.0);
            s.m := 1;
            s.v := bb;
            Append(first,last,s);
          end;
        end if;
      exception
        when constraint_error => return;
      end;
    else
      declare
        h : Standard_Complex_Vectors.Vector(0..n);
        count : natural32 := 0;
        z : Partition(1..natural32(n));
        m : natural32;
        d : Standard_Natural_Vectors.Vector(1..n);
      begin
        Degree_Structure.Get(natural32(i),z,d);
        m := Degree_Structure.Get(natural32(i));
        for j1 in 1..m loop
          if Degree_Structure.Admissible(acc,natural32(i)-1,z(j1)) then
            acc(natural32(i)) := Create(z(j1));
            for j2 in 1..d(integer32(j1)) loop
              count := count + 1;
              h := Standard_Linear_Product_System.Get_Hyperplane
                     (natural32(i),count);
              b(i) := -h(0);
              for k in 1..n loop
                a(i,k) := h(k);
              end loop;
              Solve(i+1,n,first,last,a,b,acc);
            end loop;
            Clear(acc(natural32(i)));
          else 
            count := count + d(integer32(j1));
          end if;
        end loop;
        Clear(z);
      end;
    end if;
  end Solve;

  procedure Solve_Start_System 
                ( n : in natural32; sols : in out Solution_List ) is

    m : Matrix(1..integer32(n),1..integer32(n));
    v : Standard_Complex_Vectors.Vector(1..integer32(n));
    acc : Partition(1..n);
    last : Solution_List;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        m(i,j) := Create(0.0);
      end loop;
      v(i) := Create(0.0);
    end loop;
    Solve(1,integer32(n),sols,last,m,v,acc);
  end Solve_Start_System;

  procedure GBQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List ) is

    n : constant natural32 := natural32(p'length);

  begin
    Standard_Linear_Product_System.Init(n);
    Construct_Random_Product_System(p);
    q := Standard_Linear_Product_System.Polynomial_System;
    Solve_Start_System(n,sols);
    Degree_Structure.Clear;
    Standard_Linear_Product_System.Clear;
  end GBQ;

  procedure GBQ ( p : in Poly_Sys; q : out Poly_Sys; rq : out Prod_Sys;
                  sols : in out Solution_List ) is

    n : constant natural32 := natural32(p'length);

  begin
    Standard_Linear_Product_System.Init(n);
    Construct_Random_Product_System(p);
    q := Standard_Linear_Product_System.Polynomial_System;
    rq := Standard_Complex_Prod_Planes.Create;
    Solve_Start_System(n,sols);
    Degree_Structure.Clear;
    Standard_Linear_Product_System.Clear;
  end GBQ;

end Multi_Homogeneous_Start_Systems;
