with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_System_and_Solutions_io;
with QuadDobl_Homotopy;
with Standard_Parameter_Systems;
with Solution_Drops;
with QuadDobl_Homotopy_Convolutions_io;
with QuadDobl_Newton_Convolutions;
with QuadDobl_Newton_Convolution_Steps;
with Convergence_Radius_Estimates;
with Fabry_on_Homotopy_Helpers;

package body QuadDobl_Fabry_on_Homotopy is

  procedure QuadDobl_Newton_Fabry
              ( cfs : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in QuadDobl_Complex_Vectors.Vector ) is

    dim : constant integer32 := sol'last;
    deg : constant integer32 := cfs.deg;
    scf : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
        := QuadDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
    maxit : integer32 := deg/2;
    nbrit : integer32 := 0;
    tol : double_float := 1.0E-32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector(1..dim); -- dim = #equations
    fail,verbose : boolean;
    absdx,rcond,rad,err : quad_double;
    scale : constant boolean := false;
    zpt : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    Fabry_on_Homotopy_Helpers.Prompt_for_Parameters(maxit,tol,verbose);
    if verbose then
      QuadDobl_Newton_Convolution_Steps.LU_Newton_Steps
        (standard_output,cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,
         ipvt,wrk,scale);
    else
      QuadDobl_Newton_Convolution_Steps.LU_Newton_Steps
        (cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale,false);
    end if;
    put_line("The coefficients of the series : "); put_line(scf);
    Convergence_Radius_Estimates.Fabry
      (standard_output,scf,zpt,rad,err,fail,1,true);
    Fabry_on_Homotopy_Helpers.Write_Report(standard_output,rad,err,zpt,fail);
    QuadDobl_Complex_Vectors.Clear(wrk);
  end QuadDobl_Newton_Fabry;

  procedure QuadDobl_Run
              ( nbequ,idxpar,deg : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    cvh : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    ans : character;

  begin
    cvh := QuadDobl_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Newton_Fabry(cvh,ls.v);
      put("Continue with the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end QuadDobl_Run;

  procedure QuadDobl_Artificial_Setup is

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;
    nbequ,nbvar,nbsols,deg : integer32 := 0;
    ans : character;

    use QuadDobl_Complex_Polynomials;

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
    QuadDobl_System_and_Solutions_io.get(start,sols);
    nbsols := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      nbvar := QuadDobl_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(nbvar,1); put_line(".");
      new_line;
      put("Random gamma ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then gamma := QuadDobl_Random_Numbers.Random1;
       else gamma := QuadDobl_Complex_Numbers.Create(integer(1));
      end if;
      QuadDobl_Homotopy.Create(target.all,start.all,1,gamma);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      QuadDobl_Run(nbequ,nbvar+1,deg,sols);
    end if;
  end QuadDobl_Artificial_Setup;

  procedure QuadDobl_Natural_Setup is

    hom : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : QuadDobl_Complex_Solutions.Solution_List;
    nbequ,sysnbvar,solnbvar,nbsols,idxpar,deg : integer32 := 0;
    par : Standard_Integer_Vectors.Vector(1..1);

    use QuadDobl_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for a homotopy ...");
    QuadDobl_System_and_Solutions_io.get(hom,sols);
    nbequ := hom'last;
    sysnbvar := integer32(Number_of_Unknowns(hom(hom'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(sysnbvar,1); put_line(" variables.");
    nbsols := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      solnbvar := QuadDobl_Complex_Solutions.Head_Of(sols).n;
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
      QuadDobl_Homotopy.Create(hom.all,idxpar);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      if solnbvar = nbequ
       then QuadDobl_Run(nbequ,idxpar,deg,sols);
       else QuadDobl_Run(nbequ,idxpar,deg,dropsols);
      end if;
    end if;
  end QuadDobl_Natural_Setup;

  procedure Main is

    ans : character;

  begin
    new_line;
    put("Artificial-parameter homotopy ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Artificial_Setup;
     else QuadDobl_Natural_Setup;
    end if;
  end Main;

end QuadDobl_Fabry_on_Homotopy;
