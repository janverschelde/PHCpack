with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Double_Real_Powered_Series;
with Random_Laurent_Homotopy;
with Double_Newton_Puiseux;

package body Test_Newton_Puiseux is

  procedure Define_Homotopy
              ( dim : in integer32;
                nbm,nbt : in Standard_Integer_Vectors.Vector;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                hct : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    pdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);
    pcf : Standard_Complex_VecVecs.VecVec(1..dim);
    pct : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    Random_Laurent_Homotopy.Random_Laurent_System
      (dim,dim,-9,9,nbm,pdg,pcf,pct);
    if vrblvl > 0 then
      for i in 1..dim loop
        put("-> coefficients and degrees of polynomial ");
        put(i,1); put_line(" :");
        for j in 1..nbm(i) loop
          put(pcf(i)(j)); 
          put(" t^"); put(pct(i)(j));
          put("  "); put(pdg(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    Double_Real_Powered_Series.Random_Power_Series(dim,nbt,cff,pwr);
    if vrblvl > 0 then
      for i in 1..dim loop
        put("-> a random power series "); put(i,1); put_line(" : ");
        put(cff(i)(0)); new_line;
        for j in 1..nbt(i) loop
          put(cff(i)(j)); put("  t^"); put(pwr(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    Random_Laurent_Homotopy.Random_Homotopy(pdg,pcf,pct,cff,pwr,hdg,hcf,hct,1);
    Random_Laurent_Homotopy.Scale_Homotopy_Powers(hct);
    if vrblvl > 0 then
      for i in 1..dim loop
        put("-> coefficients and degrees of homotopy ");
        put(i,1); put_line(" :");
        for j in hcf(i)'range loop
          put(hcf(i)(j)); 
          put(" t^"); put(hct(i)(j));
          put("  "); put(hdg(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    Random_Laurent_Homotopy.Test_Random_Homotopy(hdg,hcf,hct,cff,pwr);
  end Define_Homotopy;

  procedure Test ( dim : in integer32 ) is

    nbm : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    nbt : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    cff : Standard_Complex_VecVecs.VecVec(1..dim);
    pwr : Standard_Floating_VecVecs.VecVec(1..dim);
    hdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);
    hcf : Standard_Complex_VecVecs.VecVec(1..dim);
    hct : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    put_line("Reading the number of monomials for every polynomial ...");
    for i in 1..dim loop
      put("  Give the number of monomials in polynomial "); put(i,1);
      put(" : "); get(nbm(i));
    end loop;
    put_line("Reading the number of terms in every series ...");
    for i in 1..dim loop
      put("  Give the number of terms in series "); put(i,1);
      put(" : "); get(nbt(i));
    end loop;
    Define_Homotopy(dim,nbm,nbt,cff,pwr,hdg,hcf,hct,1);
    Double_Newton_Puiseux.Run_Newton_Step(hcf,hct,hdg,cff,pwr,vrblvl=>4);
  end Test;

  procedure Main is

    dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables : "); get(dim);
    Test(dim);
  end Main;

end Test_Newton_Puiseux;
