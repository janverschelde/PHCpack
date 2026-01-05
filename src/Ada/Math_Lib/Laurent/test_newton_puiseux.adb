with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
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

  procedure Run_Diagonal_Newton_Steps
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float := 1.0E-12 ) is

  -- DESCRIPTION :
  --   Calls the non-interactive version of the diagonal Newton steps
  --   and checks the errors.

  -- ON ENTRY :
  --   hdg      supports of the Laurent homotopy;
  --   hcf      coefficients of the polynomials in the homotopy;
  --   hct      powers of t in the homotopy for each monomial;
  --   cff      coefficients of a power series solution;
  --   pwr      exponents of a power series solution;
  --   tol      tolerance to decide if a number is zero.

    dim : constant integer32 := cff'last;
    cf0,cf1,cf2,cf3,cf4 : Standard_Complex_Vectors.Vector(1..dim);
    pw1,pw2,pw3,pw4 : Standard_Floating_Vectors.Vector(1..dim);
    ans : character;
    nbr,vrb : integer32 := 0;
    zero : constant character := '0';
    err,sumerr : double_float;

  begin
    for i in cf0'range loop
      cf0(i) := cff(i)(cff(i)'first);
    end loop;
    put("Give the number of Newton steps, 1, 2, 3, or 4 : ");
    Communications_with_User.Ask_Alternative(ans,"1234");
    nbr := integer32(character'pos(ans)) - integer32(character'pos(zero));
    put("Give the verbose level (0 is silent) : "); get(vrb);
    put("running "); put(nbr,1); put_line(" steps ...");
    Double_Newton_Puiseux.Diagonal_Newton_Steps
      (hcf,hct,hdg,cf0,nbr,cf1,cf2,cf3,cf4,pw1,pw2,pw3,pw4,tol,vrb);
    sumerr := 0.0;
    for i in 1..nbr loop
      put("coefficients and powers at "); put(i,1); put_line(" :");
      for k in 1..dim loop
        put(cff(k)(i)); put("  t^"); put(pwr(k)(i)); new_line;
        case i is
          when 1 => put(cf1(k)); put("  t^"); put(pw1(k)); new_line;
                    err := AbsVal(cff(k)(i) - cf1(k));
                    err := err + abs(pwr(k)(i) - pw1(k));
          when 2 => put(cf2(k)); put("  t^"); put(pw2(k)); new_line;
                    err := AbsVal(cff(k)(i) - cf2(k));
                    err := err + abs(pwr(k)(i) - pw2(k));
          when 3 => put(cf3(k)); put("  t^"); put(pw3(k)); new_line;
                    err := AbsVal(cff(k)(i) - cf3(k));
                    err := err + abs(pwr(k)(i) - pw3(k));
          when 4 => put(cf4(k)); put("  t^"); put(pw4(k)); new_line;
                    err := AbsVal(cff(k)(i) - cf4(k));
                    err := err + abs(pwr(k)(i) - pw4(k));
          when others => null;
        end case;
        put("error : "); put(err,3); new_line;
        sumerr := sumerr + err;
      end loop;
    end loop;
    put("sum of errors : "); put(sumerr,3);
    if sumerr < tol
     then put_line(" okay!");
     else put_line(" bug?!");
    end if;
  end Run_Diagonal_Newton_Steps;

  procedure Test ( dim : in integer32 ) is

    nbm : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    nbt : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    cff : Standard_Complex_VecVecs.VecVec(1..dim);
    pwr : Standard_Floating_VecVecs.VecVec(1..dim);
    hdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);
    hcf : Standard_Complex_VecVecs.VecVec(1..dim);
    hct : Standard_Floating_VecVecs.VecVec(1..dim);
    ans : character;

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
    new_line;
    put("Run interactive version ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y'
     then Double_Newton_Puiseux.Run_Newton_Step(hcf,hct,hdg,cff,pwr,vrblvl=>4);
     else Run_Diagonal_Newton_Steps(hcf,hct,hdg,cff,pwr);
    end if;
  end Test;

  procedure Main is

    dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables : "); get(dim);
    Test(dim);
  end Main;

end Test_Newton_Puiseux;
