with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

function Solve_Pieri_Leaves
           ( file : file_type; b1,b2 : Bracket; m : Matrix ) return Matrix is

  res : Matrix(m'range(1),b1'range);
  n : constant integer32 := m'last(1);
  dim : constant integer32 := m'last(2)-1;  -- dimension of linear system
  equ : Bracket(1..dim);                    -- equations : res(equ(i),*) = 0
  sys : Matrix(1..dim,1..dim);
  rhs : Standard_Complex_Vectors.Vector(1..dim);
  gen : Standard_Complex_Vectors.Vector(1..n);

  procedure Initialize_Solution_Plane is

  -- DESCRIPTION :
  --   Fills in the patterns of the resulting p-plane and sets up the
  --   equations that the resulting plane has to satisfy.
  --   These equations define the space C.

    cnt : integer32 := 0;
    zero_row : boolean;

  begin
    for i in 1..n loop
      zero_row := true;
      for j in b1'range loop
        if ((i < integer32(b2(j))) or (i > n+1-integer32(b1(b1'last+1-j))))
         then null;
         else zero_row := false;
        end if;
      end loop;
      if zero_row then
        cnt := cnt+1;
        equ(cnt) := natural32(i);
      end if;
    end loop;
  end Initialize_Solution_Plane;

  procedure Solve_Linear_System is

  -- DESCRIPTION :
  --   Solves the linear system that determines the solution plane.
  --   The solution to the linear system gives the coefficients in the
  --   linear combination of the columns of the given m-plane that
  --   determines the generator for the line of intersection with the
  --   space C and the given m-plane.

    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
   -- orgsys : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
   -- orgrhs,resi : Standard_Complex_Vectors.Vector(1..dim);
    use Standard_Complex_Vectors;

  begin
    for j in 1..dim loop
      for i in equ'range loop
        sys(i,j) := m(integer32(equ(i)),j);
      end loop;
    end loop;
   -- orgsys := sys;
    for i in 1..dim loop
      rhs(i) := m(integer32(equ(i)),m'last(2));
    end loop;
   -- orgrhs := rhs;
    lufac(sys,dim,ipvt,info);
    lusolve(sys,dim,ipvt,rhs);
    for i in 1..n loop       -- linear combination of the columns of m-plane
      gen(i) := -m(i,m'last(2));
      for j in rhs'range loop
        gen(i) := gen(i) + rhs(j)*m(i,j);
      end loop;
    end loop;
   -- put_line(file,"The solution to the linear system : ");
   -- put_line(file,rhs); new_line(file);
   -- resi := orgrhs - orgsys*rhs;
   -- put_line(file,"The residual vector of the linear system : ");
   -- put(file,resi,2); new_line(file);
   -- put_line(file,"The generator of the intersection : ");
   -- put_line(file,gen); new_line(file);
  end Solve_Linear_System;

  procedure Determine_Solution_Plane is

  -- DESCRIPTION :
  --   Determines the solution plane from the generator of the intersection
  --   of the space C with the given m-plane.

  begin
    for i in 1..n loop
      for j in b1'range loop
        if ((i < integer32(b2(j))) or (i > n+1-integer32(b1(b1'last+1-j))))
         then res(i,j) := Create(0.0);
         else res(i,j) := gen(i);
        end if;
      end loop;
    end loop;
  end Determine_Solution_Plane;

begin
  Initialize_Solution_Plane;
  Solve_Linear_System;
  Determine_Solution_Plane;
  return res;
end Solve_Pieri_Leaves;
