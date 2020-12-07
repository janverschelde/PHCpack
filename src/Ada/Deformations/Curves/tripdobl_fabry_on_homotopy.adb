with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with TripDobl_Complex_Numbers;
with TripDobl_Random_Numbers;
with Standard_Integer_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_VecVecs_io;        use TripDobl_Complex_VecVecs_io;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with TripDobl_System_and_Solutions_io;
with TripDobl_Homotopy;
with Standard_Parameter_Systems;
with Solution_Drops;
with TripDobl_Homotopy_Convolutions_io;
with TripDobl_Newton_Convolutions;
with TripDobl_Newton_Convolution_Steps;
with Convergence_Radius_Estimates;
with Fabry_on_Homotopy_Helpers;

package body TripDobl_Fabry_on_Homotopy is

  procedure TripDobl_Newton_Fabry
              ( cfs : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in TripDobl_Complex_Vectors.Vector ) is

    dim : constant integer32 := sol'last;
    deg : constant integer32 := cfs.deg;
    scf : constant TripDobl_Complex_VecVecs.VecVec(1..dim)
        := TripDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
    maxit : integer32 := deg/2;
    nbrit : integer32 := 0;
    tol : double_float := 1.0E-20;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : TripDobl_Complex_Vectors.Link_to_Vector
        := new TripDobl_Complex_Vectors.Vector(1..dim); -- dim = #equations
    fail,verbose : boolean;
    absdx,rcond,rad,err : triple_double;
    scale : constant boolean := false;
    zpt : TripDobl_Complex_Numbers.Complex_Number;

  begin
    Fabry_on_Homotopy_Helpers.Prompt_for_Parameters(maxit,tol,verbose);
    if verbose then
      TripDobl_Newton_Convolution_Steps.LU_Newton_Steps
        (standard_output,cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,
         ipvt,wrk,scale);
    else
      TripDobl_Newton_Convolution_Steps.LU_Newton_Steps
        (cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale,false);
    end if;
    put_line("The coefficients of the series : "); put_line(scf);
    Convergence_Radius_Estimates.Fabry
      (standard_output,scf,zpt,rad,err,fail,1,true);
    Fabry_on_Homotopy_Helpers.Write_Report(standard_output,rad,err,zpt,fail);
    TripDobl_Complex_Vectors.Clear(wrk);
  end TripDobl_Newton_Fabry;

  procedure TripDobl_Run
              ( nbequ,idxpar,deg : in integer32;
                sols : in out TripDobl_Complex_Solutions.Solution_List ) is

    cvh : TripDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : TripDobl_Complex_Solutions.Solution_List := sols;
    ls : TripDobl_Complex_Solutions.Link_to_Solution;
    ans : character;

  begin
    cvh := TripDobl_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    while not TripDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := TripDobl_Complex_Solutions.Head_Of(tmp);
      TripDobl_Newton_Fabry(cvh,ls.v);
      put("Continue with the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := TripDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end TripDobl_Run;

  procedure TripDobl_Artificial_Setup is

    target,start : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : TripDobl_Complex_Solutions.Solution_List;
    gamma : TripDobl_Complex_Numbers.Complex_Number;
    nbequ,nbvar,nbsols,deg : integer32 := 0;
    ans : character;

    use TripDobl_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for the target system ...");
    get(target);
    nbequ := target'last;
    nbvar := integer32(Number_of_Unknowns(target(target'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(nbvar,1); put_line(" variables.");
    new_line;
    put_line
      ("Reading the name of a file for the start system and solutions ...");
    TripDobl_System_and_Solutions_io.get(start,sols);
    nbsols := integer32(TripDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      nbvar := TripDobl_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(nbvar,1); put_line(".");
      new_line;
      put("Random gamma ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then gamma := TripDobl_Random_Numbers.Random1;
       else gamma := TripDobl_Complex_Numbers.Create(integer(1));
      end if;
      TripDobl_Homotopy.Create(target.all,start.all,1,gamma);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      TripDobl_Run(nbequ,nbvar+1,deg,sols);
    end if;
  end TripDobl_Artificial_Setup;

  procedure TripDobl_Natural_Setup is

    hom : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : TripDobl_Complex_Solutions.Solution_List;
    nbequ,sysnbvar,solnbvar,nbsols,idxpar,deg : integer32 := 0;
    par : Standard_Integer_Vectors.Vector(1..1);

    use TripDobl_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for a homotopy ...");
    TripDobl_System_and_Solutions_io.get(hom,sols);
    nbequ := hom'last;
    sysnbvar := integer32(Number_of_Unknowns(hom(hom'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(sysnbvar,1); put_line(" variables.");
    nbsols := integer32(TripDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      solnbvar := TripDobl_Complex_Solutions.Head_Of(sols).n;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(solnbvar,1); put_line(".");
      new_line;
      par := Standard_Parameter_Systems.Define_Parameters(nbequ,sysnbvar,1);
      idxpar := par(1);
      put("The index to the continuation parameter : ");
      put(idxpar,1); new_line;
      if solnbvar = nbequ then
        put_line("Solution dimension is okay.");
      else
        put_line("Need to drop one coordinate of each solution.");
        dropsols := Solution_Drops.Drop(sols,natural32(idxpar));
      end if;
      TripDobl_Homotopy.Create(hom.all,idxpar);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      if solnbvar = nbequ
       then TripDobl_Run(nbequ,idxpar,deg,sols);
       else TripDobl_Run(nbequ,idxpar,deg,dropsols);
      end if;
    end if;
  end TripDobl_Natural_Setup;

  procedure Main is

    ans : character;

  begin
    new_line;
    put("Artificial-parameter homotopy ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then TripDobl_Artificial_Setup;
     else TripDobl_Natural_Setup;
    end if;
  end Main;

end TripDobl_Fabry_on_Homotopy;
