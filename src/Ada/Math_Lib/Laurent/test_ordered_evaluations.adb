with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_IO;       use Standard_Integer_Vectors_IO;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;
with Random_Laurent_Homotopy;
with Double_Ordered_Evaluations;

package body Test_Ordered_Evaluations is

  procedure Test ( dim,nbr,ord : in integer32 ) is

    cff : Standard_Complex_VecVecs.VecVec(1..dim);
    pwr : Standard_Floating_VecVecs.VecVec(1..dim);
    nbt : constant Standard_Integer_Vectors.Vector(1..dim)
        := (1..dim => ord);
    pdg : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Laurent_Homotopy.Random_Polynomial_Support(nbr,dim,-9,9);
    pcf : constant Standard_Complex_Vectors.Vector(1..nbr)
        := Standard_Random_Vectors.Random_Vector(1,nbr);
    pct : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Random_Laurent_Homotopy.Random_Leading_Powers(nbr);
    ycf : Standard_Complex_Vectors.Vector(1..(dim+1)*nbr);
    ydg : Standard_Floating_Vectors.Vector(1..(dim+1)*nbr);

    use Double_Ordered_Evaluations;

  begin
    Random_Laurent_Homotopy.Random_Power_Series(dim,nbt,cff,pwr);
    for i in cff'range loop
      put("-> a random power series "); put(i,1); put_line(" : ");
      put(cff(i)(0)); new_line;
      for j in 1..nbt(i) loop
        put(cff(i)(j)); put("  t^"); put(pwr(i)(j)); new_line;
      end loop;
    end loop;
    put_line("A random polynomial :");
    for i in 1..nbr loop
      put(pcf(i)); 
      put(" t^"); put(pct(i));
      put("  "); put(pdg(i)); new_line;
    end loop;
    First_Order_Evaluation(pcf,pct,pdg,cff,pwr,ycf,ydg,1);
    put_line("the first order evaluation :");
    for i in ycf'range loop
       put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
    end loop;
  end Test;

  procedure Main is

    dim,nbr,ord : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the number of monomials : "); get(nbr);
    put("Give the order of the series : "); get(ord);
    Test(dim,nbr,ord);
  end Main;

end Test_Ordered_Evaluations;
