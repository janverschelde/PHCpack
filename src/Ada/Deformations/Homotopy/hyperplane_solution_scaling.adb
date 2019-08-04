with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vector_Norms;

package body Hyperplane_Solution_Scaling is

  procedure Sub ( p : in out Standard_Complex_Polynomials.Poly;
                  c : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Polynomials;

    t : Term;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    d : constant Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);

  begin
    t.cf := c;
    t.dg := new Standard_Natural_Vectors.Vector'(d);
    Sub(p,t);
    Clear(t);
  end Sub;

  procedure Scale ( v : in out Standard_Complex_Vectors.Vector ) is

    mxv : constant double_float := Standard_Complex_Vector_Norms.Max_Norm(v);

    use Standard_Complex_Numbers;

  begin
    for i in v'range loop
      v(i) := v(i)/mxv;
    end loop;
  end Scale;

  procedure Adjust ( c : in Standard_Complex_Vectors.Link_to_Vector;
                     v : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    val : Complex_Number := c(c'last);

  begin
    for i in v'range loop
      val := val + c(i)*v(i);
    end loop;
    c(c'last) := c(c'last) - val;
  end Adjust;

end Hyperplane_Solution_Scaling;
