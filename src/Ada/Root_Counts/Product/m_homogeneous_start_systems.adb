with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Integer_Matrices;
with Lists_of_Integer_Vectors;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;        use Degrees_in_Sets_of_Unknowns;
with Standard_Linear_Product_System;
with Standard_Complex_Prod_Planes;
with m_Homogeneous_Permanent_Factors;    use m_Homogeneous_Permanent_Factors;

package body m_Homogeneous_Start_Systems is

  procedure Create_Random_Hyperplanes
               ( index,n,d : in integer32; s : in Set ) is
  begin
    for i in 1..d loop
      declare
        h : Standard_Complex_Vectors.Vector(0..n);
      begin
        h(0) := Random1;
        for j in 1..Dimension(s) loop
          if Is_In(s,j)
           then h(integer32(j)) := Random1;
           else h(integer32(j)) := Create(0.0);
          end if;
        end loop;
        Standard_Linear_Product_System.Add_Hyperplane(natural32(index),h);
      end;
    end loop;
  end Create_Random_Hyperplanes;

  procedure Create_Random_System 
              ( n,m : natural32; z : partition;
                d : Standard_Integer_Matrices.Matrix ) is

  begin
    for j in d'range(2) loop
      for i in d'range(1) loop
        Create_Random_Hyperplanes(i,integer32(n),d(i,j),z(natural32(j)));
      end loop;
    end loop;
  end Create_Random_System;

  procedure m_Homogeneous_Start_System
              ( p : in Poly_Sys; z : in partition ) is

    n : constant natural32 := natural32(p'length);
    m : constant natural32 := z'last;
    d : constant Standard_Integer_Matrices.Matrix := Degree_Table(p,z);

  begin
    Standard_Linear_Product_System.Init(n);
    Create_Random_System(n,m,z,d);
  end m_Homogeneous_Start_System;

  procedure m_Homogeneous_Start_System
                 ( p : in Poly_Sys; z : in partition;
                   q : out Poly_Sys; qsols : in out Solution_List ) is

    n : constant natural32 := natural32(p'length);
    m : constant natural32 := z'last;
    d : constant Standard_Integer_Matrices.Matrix := Degree_Table(p,z);
   -- nl : natural32 := 0;
    ind_sols : Lists_of_Integer_Vectors.List;

  begin
    Standard_Linear_Product_System.Init(n);
    Create_Random_System(n,m,z,d);
   -- using Solve of Standard_Linear_Product_System is too inefficient
   -- q := Standard_Linear_Product_System.Polynomial_System;
   -- Standard_Linear_Product_System.Solve(qsols,nl);
    Permanent_Factors(p,z,ind_sols);
    Solve_m_Homogeneous_Start_System(ind_sols,q,qsols);
    Lists_of_Integer_Vectors.Clear(ind_sols);
    Standard_Linear_Product_System.Clear;
  end m_Homogeneous_Start_System;

  procedure m_Homogeneous_Start_System
                 ( p : in Poly_Sys; z : in partition;
                   q : out Poly_Sys; rq : out Prod_Sys;
                   qsols : in out Solution_List ) is

    n : constant natural32 := natural32(p'length);
    m : constant natural32 := z'last;
    d : constant Standard_Integer_Matrices.Matrix := Degree_Table(p,z);
   -- nl : natural32 := 0;
    ind_sols : Lists_of_Integer_Vectors.List;

  begin
    Standard_Linear_Product_System.Init(n);
    Create_Random_System(n,m,z,d);
    Permanent_Factors(p,z,ind_sols);
    Solve_m_Homogeneous_Start_System(ind_sols,q,qsols);
    rq := Standard_Complex_Prod_Planes.Create;
   -- using Solve of Standard_Linear_Product_System is too inefficient
   -- q := Standard_Linear_Product_System.Polynomial_System;
   -- Standard_Linear_Product_System.Solve(qsols,nl);
    Lists_of_Integer_Vectors.Clear(ind_sols);
    Standard_Linear_Product_System.Clear;
  end m_Homogeneous_Start_System;

end m_Homogeneous_Start_Systems;
