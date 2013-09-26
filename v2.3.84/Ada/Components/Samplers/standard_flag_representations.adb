with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Plane_Representations;    use Standard_Plane_Representations;

package body Standard_Flag_Representations is

  function Create_Intrinsic ( f : Matrix ) return Matrix is

    n : constant integer32 := f'last(2);
    k : constant integer32 := f'last(1);
    res : Matrix(1..n,0..n);
    smallest : constant Matrix(1..n,0..n-k) := Generators(f);
    c : Complex_Number;

  begin
    for i in 1..n loop               -- initialize with smallest space
      for j in 0..n-k loop
        res(i,j) := smallest(i,j);
      end loop;
    end loop;
    for d in reverse 1..k-1 loop     -- add spaces of co-dimension d
      declare
        p : Matrix(1..d,0..n);       -- extrinsic plane
        g : Matrix(1..n,0..n-d);     -- intrinsic representation
      begin
        for i in 1..d loop           -- copy d hyperplanes from f
          for j in 0..n loop
            p(i,j) := f(i,j);
          end loop;
        end loop;
        g := Generators(p);          -- compute intrinsic representation
        for i in 1..n loop 
          res(i,n-d) := Create(0.0); -- initialize new direction
        end loop;
        for j in 1..n-d loop
          c := Random1;              -- random direction
          for i in 1..n loop
            res(i,n-d) := res(i,n-d) + c*g(i,j);
          end loop;
        end loop;
      end;
    end loop;
    for i in 1..n loop               -- last direction is arbitrary
      res(i,n) := Random1;
    end loop;
    return Standard_Plane_Representations.Orthogonalize(res);
  end Create_Intrinsic;

  procedure Test_Flag ( file : in file_type;
                        extrinsic_flag,intrinsic_flag : in Matrix ) is

    n : constant integer32 := extrinsic_flag'last(2);
    k : constant integer32 := extrinsic_flag'last(1);
    x : Vector(1..n);
    c : Complex_Number;
    y : Vector(1..k);

  begin
    new_line(file);
    put_line(file,
             "Evaluation of flag at random points from intrinsic flag...");
    for d in n-k..n loop                -- space of dimension d
      for i in 1..n loop 
        x(i) := intrinsic_flag(i,0);    -- initialize with offset vector
      end loop;
      for j in 1..d loop                -- add in random directions
        c := Random1;
        for i in 1..n loop
          x(i) := x(i) + c*intrinsic_flag(i,j);
        end loop;
      end loop;                         -- x is random point in d-space
      for i in 1..k loop                -- evaluate in hyperplanes
        y(i) := extrinsic_flag(i,0);
        for j in 1..n loop
          y(i) := y(i) + extrinsic_flag(i,j)*x(j);
        end loop;
      end loop;
      put(file,"Evaluation of point in ");
      put(file,d,1); put_line(file,"-space :");
      put_line(file,y);
    end loop;
  end Test_Flag;

end Standard_Flag_Representations;
