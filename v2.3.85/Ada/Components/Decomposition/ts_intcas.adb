with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Continuation_Parameters;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Flag_Representations;     use Standard_Flag_Representations;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Standard_Intrinsic_Continuation;   use Standard_Intrinsic_Continuation;

procedure ts_intcas is

-- DESCRIPTION :
--   Interactive development of a cascade of homotopies
--   using intrinsic coordinates.

  function Embed ( f : in Matrix; d : in natural ) return Matrix is

  -- DESCRIPTION :
  --   Adds d rows to the matrix f, embedding the flag to deal with
  --   components of dimension d and less.

    res : Matrix(1..f'last(1)+d,f'range(2));

  begin
    for i in f'range(1) loop
      for j in f'range(2) loop
        res(i,j) := f(i,j);
      end loop;
    end loop;
    for i in f'last(1)+1..f'last(1)+d loop
      for j in f'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for i in 1..d loop
      res(f'last(1)+i,f'last(2)-d+i) := Create(1.0);
    end loop;
    return res;
  end Embed;

  function Intrinsic_Plane ( s : VecVec; n,d : natural ) return Matrix is

  -- DESCRIPTION :
  --   Returns generators for the plane in the ambient n-space
  --   of the embedding of a solution set of dimension d.
   
    res : Matrix(1..n,0..n-d);

  begin
    for i in 1..n-d loop
      res(i,0) := Create(0.0);
      for j in 1..n-d loop
        if i = j
         then res(i,j) := Create(1.0);
         else res(i,j) := Create(0.0);
        end if;
      end loop;
    end loop;
    for i in 1..d loop
      for j in 0..n-d loop
        res(n-d+i,j) := -s(i)(j);
      end loop;
    end loop;
    return res;
  end Intrinsic_Plane;

  procedure Down_Continuation
              ( file : in file_type; ep : in Poly_Sys; d : in natural;
                f : in Matrix; sols : in Solution_List ) is

    n : constant natural := ep'last;
    k : constant natural := n - d;
    target : Poly_Sys(ep'range) := Remove_Slice(ep);
    rsols : Solution_List := Remove_Embedding(sols,d);
    esols : Solution_List;
   -- iflag : Matrix(1..f'last(1)+d,f'range(2)) := Embed(f,d);
    start_sli : VecVec(1..d) := Slices(ep,d);
    start_pla : Matrix(1..n,0..k) := Intrinsic_Plane(start_sli,n,d);
    target_sli : VecVec(1..d) := Slices(target,d);
    target_pla : Matrix(1..n,0..k) := Intrinsic_Plane(target_sli,n,d);
    epv : Eval_Poly_Sys(ep'range) := Create(ep);
    sjm : Jaco_Mat(ep'range,ep'range) := Create(ep);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    sa : Solu_Info_Array(1..Length_Of(sols));
    pp : Continuation_Parameters.Pred_Pars;
    cp : Continuation_Parameters.Corr_Pars;
    oc : natural;
    gamma : Complex_Number := Random1;

    function Path ( t : Complex_Number ) return Matrix is

     -- res : Matrix(1..n,0..k);
     -- t1 : Complex_Number := Create(1.0) - t;

    begin
     --  for i in 1..n loop
     --    for j in 0..k loop
     --      res(i,j) := t1*start_pla(i,j) + t*target_pla(i,j);
     --    end loop;
     --  end loop;
     -- return res;
      return Moving_Plane(start_pla,target_pla,gamma,t);
    end Path;
    procedure Cont is new Reporting_QR_Continue(Path);

  begin
    new_line(file);
    put_line(file,"The target system in the cascade :");
    put_line(file,target);
    new_line(file);
    put_line(file,"Generators of start plane :"); put(file,start_pla,2);
    put_line(file,"Generators of target plane :"); put(file,target_pla,2);
    Set_Continuation_Parameter(rsols,Create(0.0));
    esols := Expand(rsols,start_pla);
    new_line(file);
    put_line(file,"THE START SOLUTIONS in extrinsic coordinates : ");
    put(file,Length_Of(esols),Head_Of(esols).n,esols);
    sa := Deep_Create(rsols);
    Continuation_Parameters.Tune(2);
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    Driver_for_Process_io(file,oc);
    Cont(file,epv,sjf,sa,pp,cp);
  end Down_Continuation;

  procedure Run_Cascade
              ( file : in file_type; ep : in Poly_Sys;
                esols : in Solution_List; d : in natural ) is

  -- DESCRIPTION :
  --   Runs a cascade of homotopies, starting at the solutions
  --   of an embedded system at the top dimension.

  -- ON ENTRY :
  --   file     file for intermediate output and results;
  --   ep       embedded polynomial system in extrinsic coordinates;
  --   esols    solutions of ep in extrinsic coordinates;
  --   d        dimension of the top dimensional component.

    m : constant natural := ep'last;   -- dimension with embedding
    k : constant natural := m-d;       -- co-dimension
    n : constant natural := k;         -- assume complete intersection
   -- p : constant Poly_Sys := Remove_Embedding(ep,d);
    s : VecVec(1..d) := Slices(ep,d);
    eflag : Matrix(1..d,0..n) := Equations_to_Matrix(s,n);
    iflag : Matrix(1..n,0..n) := Create_Intrinsic(eflag);
   -- isols : Solution_List := Project(esols,iflag);

  begin
   -- put_line(file,"The original system : "); put(file,p);
   -- Test_Flag(file,eflag,iflag);
   -- new_line(file);
   -- put_line(file,"THE SOLUTIONS in intrinsic coordinates : ");
   -- put(file,Length_Of(isols),Head_Of(isols).n,isols);
   -- Down_Continuation(file,ep,d,iflag,isols);
    Down_Continuation(file,ep,d,iflag,esols);
  end Run_Cascade;

  procedure Main is

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural;
    file : file_type;

  begin
    new_line;
    put_line("Executing a cascade of homotopies intrinsically...");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("Running the cascade, see the output file for results...");
    new_line;
    Run_Cascade(file,lp.all,sols,dim);
  end Main;

begin
  Main;
end ts_intcas;
