with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Series_Norms;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_Systems_io;   use Standard_CSeries_Poly_Systems_io;
-- with Complex_Series_and_Polynomials_io;
with Floating_Lifting_Utilities;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Supports_of_Polynomial_Systems;
with Random_Coefficient_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Series_and_Solutions;
with Standard_Simpomial_Solvers;
with Power_Series_Methods;
with Double_Taylor_Homotopies;           use Double_Taylor_Homotopies;
with Double_Taylor_Homotopies_io;        use Double_Taylor_Homotopies_io;
with Taylor_Homotopy_Series;
with Taylor_Polyhedral_Homotopies;       use Taylor_Polyhedral_Homotopies;

procedure ts_rptmccsol is

-- DESCRIPTION :
--   Tests the application of the robust path tracking to solve a
--   random coefficient system, given a mixed cell configuration.

  procedure Confirm_Input
              ( file : in file_type; p : in Link_to_Laur_Sys;
                n : in natural32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in out Mixed_Subdivision; mv : out natural32 ) is

  -- DESCRIPTION :
  --   Confirms the input, writing to output file,
  --   and computes the volumes of all cells.
 
  -- ON ENTRY :
  --   file     must be opened for input;
  --   p        a polynomial system;
  --   n        number of variables;
  --   mix      type of mixture;
  --   mcc      a regular mixed cell configuration.

  -- ON RETURN :
  --   mv       the mixed volume.

  begin
    put(file,n,p.all);
    new_line(file);
    put(file,n,mix.all,mcc,mv);
  end Confirm_Input;

  function Coefficient_Vectors
             ( q : Laur_Sys ) return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficient vectors of the polynomials in q.

    res : Standard_Complex_VecVecs.VecVec(q'range);

    use Standard_Complex_Laur_Functions;

  begin
    for i in q'range loop
      declare
        cff : constant Standard_Complex_Vectors.Vector := Coeff(q(i));
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Coefficient_Vectors;

  procedure Track
              ( file : in file_type; deg : in integer32;
                pnt : in double_float; q : in Laur_Sys;
                cfq : in Standard_Complex_VecVecs.VecVec;
                qsols,qsols_last : in out Solution_List;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell ) is

  -- DESCRIPTION :
  --   Tracks the paths defined by the mixed cell.
 
  -- ON ENTRY :
  --   file     must be opened for input;
  --   deg      truncation degree of the series;
  --   pnt      expansion point of the Taylor series;
  --   q        a random coefficient system;
  --   cfq      coefficient vectors of q;
  --   mv       the mixed volume.
  --   mix      type of mixture;
  --   lif      lifted supports;
  --   mic      cell in a regular mixed cell configuration.

  -- ON RETURN :
  --   q        a random coefficient system;
  --   qsols    solutions of q.

    sq : Laur_Sys(q'range)
       := Supports_of_Polynomial_Systems.Select_Terms(q,mix.all,mic.pts.all);
    sqsols : Solution_List;
    tol_zero : constant double_float := 1.0e-12;
    fail,zero_y : boolean;
    thm : Taylor_Homotopy(q'range);
    hom : Standard_CSeries_Poly_Systems.Poly_Sys(q'range);
    verbose : constant boolean := false;
    nbrit : constant integer32 := deg;

  begin
    put(file,natural32(sq'last),sq);
    Standard_Simpomial_Solvers.Solve
      (sq,tol_zero,sqsols,fail,zero_y,false);
    if fail then
      put_line(file,"Solving initial form system failed!");
    else
      put(file,"Computed "); put(file,Length_Of(sqsols),1);
      put_line(file," start solutions.");
      put(file,Length_Of(sqsols),natural32(sq'last),sqsols);
      Make_Homotopy(file,deg,pnt,cfq,mix,lif,mic,thm);
      put_line(file,"The Taylor monomial homotopy :"); put(file,thm);
      hom := Taylor_Homotopy_Series.Make(thm);
      put_line(file,"The Taylor homotopy as series system :"); put(file,hom);
      declare
        len : constant integer32 := integer32(Length_Of(sqsols));
        srv : constant Standard_Complex_Series_VecVecs.VecVec(1..len)
            := Series_and_Solutions.Create(sqsols,0);
       -- idx = 0 => copy all coordinates of the solution
        eva : Standard_Complex_Series_Vectors.Vector(q'range);
        rco,nrm : double_float;
      begin
        for i in 1..len loop
          put(file,"Running on solution "); put(file,i,1);
          put_line(file," ...");
          Power_Series_Methods.Run_SVD_Newton
            (file,deg,nbrit,hom,srv(i).all,rco,eva,verbose);
          put(file,"rcond :"); put(file,rco,3); new_line(file);
         -- put_line(file,"The evaluated solution :");
         -- Complex_Series_and_Polynomials_io.put(file,eva);
          put(file,"norms :");
          for i in eva'range loop
            nrm := Standard_Complex_Series_Norms.Max_Norm(eva(i).all);
            put(file," "); put(file,nrm,3);
          end loop;
          new_line(file);
          Standard_Complex_Series_Vectors.Clear(eva);
        end loop;
      end;
      Concat(qsols,qsols_last,sqsols);
    end if;
    Clear(sq);
    Clear(thm);
    Standard_CSeries_Poly_Systems.Clear(hom);
  end Track;

  procedure Polyhedral_Continuation
              ( file : in file_type; deg : in integer32;
                pnt : in double_float;
                q : out Link_to_Laur_Sys;
                qsols,qsols_last : out Solution_List;
                n,mv : in natural32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Tracks the paths defined by the polyhedral homotopies.
 
  -- ON ENTRY :
  --   file     must be opened for input;
  --   deg      truncation degree for the series;
  --   pnt      expansion point of the Taylor series;
  --   n        number of variables;
  --   mv       the mixed volume.
  --   mix      type of mixture;
  --   mcc      a regular mixed cell configuration.

  -- ON RETURN :
  --   q        a random coefficient system;
  --   qsols    solutions of q;
  --   qsols_last points to the last solution of qsols.

    cfq : Standard_Complex_VecVecs.VecVec(1..integer32(n));
    lif : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
        := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mcc);
    tmp : Mixed_Subdivision := mcc;
    cnt : natural32 := 0;
    mic : Mixed_Cell;

  begin
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    put(file,lif);
    q := new Laur_Sys'(Random_Coefficient_Systems.Create(n,mix.all,lif));
    cfq := Coefficient_Vectors(q.all);
    put_line(file,"A RANDOM COEFFICIENT SYSTEM :");
    put(file,n,q.all);
    new_line(file);
    put(file,"Tracking "); put(file,mv,1);
    put_line(file," solution paths ...");
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      new_line(file);
      put(file,"PROCESSING CELL "); put(file,cnt,1);
      put_line(file," :");
      mic := Head_Of(tmp);
      Track(file,deg,pnt,q.all,cfq,qsols,qsols_last,mix,lif,mic);
      tmp := Tail_of(tmp);
    end loop;
    Standard_Complex_VecVecs.Clear(cfq);
  end Polyhedral_Continuation;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for a polynomial system and
  --   for a mixed cell configuration on file.

    lp,lq : Link_to_Laur_Sys;
    qsols,qsols_last : Solution_List;
    n,m,mv : natural32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    infile,outfile : file_type;
    deg : integer32 := 0;
    pnt : double_float := 1.0e-2;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put_line("-> your polynomial system :");
    put(lp.all);
    new_line;
    put_line("Reading a file name for a mixed cell configuration ...");
    Read_Name_and_Open_File(infile);
    get(infile,n,m,mix,mcc);
    Close(infile);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    Confirm_Input(outfile,lp,n,mix,mcc,mv);
    new_line;
    put("-> the mixed volume : "); put(mv,1); new_line;
    new_line;
    put("Give the truncation degree of the series : "); get(deg);
    new_line;
    put("Give the expansion point : "); get(pnt);
    Polyhedral_Continuation(outfile,deg,pnt,lq,qsols,qsols_last,n,mv,mix,mcc);
  end Main;

begin
  Main;
end ts_rptmccsol;
