with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;    use Standard_complex_Numbers_Polar;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Sets_of_Unknowns;                  use Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;       use Degrees_in_Sets_of_Unknowns;
with Set_Structure;
with Standard_Lined_Hypersurfaces;      use Standard_Lined_Hypersurfaces;

package body Hypersurface_Samplers is

  function Random_Multihomogeneous_Directions 
               ( n : natural32; z : Partition )
               return Standard_Complex_VecVecs.VecVec is

    use Standard_Complex_Vectors;
    use Standard_Complex_VecVecs;

    res : VecVec(integer32(z'first)..integer32(z'last));

  begin
    for i in z'range loop
      declare
        v : Vector(1..integer32(n)) := (1..integer32(n) => Create(0.0));
      begin
        for j in 1..n loop
          if Is_In(z(i),j)
           then v(integer32(j)) := Random1;
          end if;
        end loop;
        res(integer32(i)) := new Vector'(v);
      end;
    end loop;
    return res;
  end Random_Multihomogeneous_Directions;

  function Random_Set_Structure_Directions
              ( n,i : natural32 ) return Standard_Complex_VecVecs.VecVec is

    use Standard_Complex_Vectors;
    use Standard_Complex_VecVecs;

    res : VecVec(1..integer32(Set_Structure.Number_of_Sets(i)));

  begin
    for j in res'range loop
      declare
        v : Vector(1..integer32(n)) := (1..integer32(n) => Create(0.0));
      begin
        for k in 1..n loop
          if Set_Structure.Is_In(i,natural32(j),k)
           then v(integer32(k)) := Random1;
          end if;
        end loop;
        res(j) := new Vector'(v);
      end;
    end loop;
    return res;
  end Random_Set_Structure_Directions;

  procedure Generic_Points
              ( p : in Standard_Complex_Polynomials.Poly;
                ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                z : in Partition;
                b : in Standard_Complex_Vectors.Vector;
                v : in Standard_Complex_VecVecs.VecVec;
                eps : in double_float; maxit : in natural32;
                t : out Standard_Complex_VecVecs.VecVec;
                fail : out boolean ) is

    use Standard_Complex_Vectors;

  begin
    fail := false;
    for i in z'range loop
      declare
        d : constant integer32 := Degree(p,z(i));
        r : Vector(1..d);
      begin
        if d > 0 then
          Generic_Points(p,ep,natural32(d),b,v(integer32(i)).all,
                         eps,maxit,r,fail);
          if not fail
           then t(integer32(i)) := new Vector'(r);
          end if;
        end if;
      end;
      exit when fail;
    end loop;
  end Generic_Points;

  procedure Generic_Points
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                z : in Partition;
                b : in Standard_Complex_Vectors.Vector;
                v : in Standard_Complex_VecVecs.VecVec;
                eps : in double_float; maxit : in natural32;
                t : out Standard_Complex_VecVecs.VecVec;
                fail : out boolean ) is

    use Standard_Complex_Vectors;

  begin
    fail := false;
    for i in z'range loop
      declare
        d : constant integer32 := Degree(p,z(i));
        r : Vector(1..d);
      begin
        if d = 0 then
          put(file,"There are no roots to compute for set ");
          put(file,i,1);
          put_line(file," with the special direction : ");
          put_line(file,v(integer32(i)).all);
        else
          put(file,"Computing "); put(file,d,1);
          put(file," roots for set "); put(file,i,1);
          put_line(file," with special direction :");
          put_line(file,v(integer32(i)).all);
          put_line(file,"Computing generic points"
                      & " on a special random line ...");
          Generic_Points(file,p,ep,natural32(d),b,v(integer32(i)).all,
                         eps,maxit,r,fail);
          put(file,"Calculation of generic points for set ");
          put(file,i,1); 
          if fail then
            put_line(file," failed.");
          else
            put_line(file," succeeded.");
            t(integer32(i)) := new Vector'(r);
          end if;
        end if;
      end;
      exit when fail;
    end loop;
  end Generic_Points;

end Hypersurface_Samplers;
