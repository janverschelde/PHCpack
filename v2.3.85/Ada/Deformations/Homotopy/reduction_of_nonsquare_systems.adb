with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers; 
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Reduction_of_Polynomials;           use Reduction_of_Polynomials;

package body Reduction_of_Nonsquare_Systems is

  function Random_Square ( p : in Poly_Sys ) return Poly_Sys is

    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    res : Poly_Sys(1..m);

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in m+1..p'last loop
        declare
          a : constant Complex_Number := Random1;
          tmp : Poly := a*p(j);
        begin
          Add(res(i),tmp);
          Clear(tmp);
        end;
      end loop;
    end loop;
    return res;
  end Random_Square;

  function Reduced_Square ( p : in Poly_Sys ) return Poly_Sys is

    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    res : Poly_Sys(1..m);

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in m+1..p'last loop
        declare
          tmp : Poly := Rpoly(res(i),p(j));
        begin
          Copy(tmp,res(i)); Clear(tmp);
        end;
      end loop;
    end loop;
    return res;
  end Reduced_Square;

end Reduction_of_Nonsquare_Systems;
