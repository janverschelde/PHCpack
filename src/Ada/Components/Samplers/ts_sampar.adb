with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Symbol_Table;                      use Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Projective_Transformations;        use Projective_Transformations;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Continuation_Parameters;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
--with Standard_Affine_Planes;            use Standard_Affine_Planes;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Affine_Solutions;         use Standard_Affine_Solutions;
with Intrinsic_Sampling_Machine;        use Intrinsic_Sampling_Machine;
--with Affine_Sampling_Machine;           use Affine_Sampling_Machine;

procedure ts_sampar is

-- DESCRIPTION :
--   Interactive development of sampling a positive dimensional component
--   using a parametric description of the linear slicing planes.

  function Square ( p : Poly_Sys; nk : integer32 ) return Poly_Sys is

  -- DESCRIPTION :
  --   The system on return has nk equations: the original first nk
  --   polynomials in p, plus random combinations of the other ones.

  begin
    if p'length <= nk then
      return p;
    else
      declare
        res : Poly_Sys(1..nk);
        r : Complex_Number;
      begin
        for i in res'range loop
          res(i) := p(i);
        end loop;
        for i in nk+1..p'last loop
          for j in res'range loop
            r := Random1;
            declare
              rp : Poly := r*p(i);
            begin
              Add(res(j),rp);
              Clear(rp);
            end;
          end loop;
        end loop;
        return res;
      end;
    end if;
  end Square;

  procedure Test_Eval
              ( n,k : in integer32; ep : in Poly_Sys; sols : in Solution_List;
                b : in Vector; v : in VecVec ) is

  -- DESCRIPTION :
  --   Tests the evaluation of the solutions, using the parametric
  --   representation of a plane of codimension k.
  --   The affine plane is represented by an offset vector b and 
  --   a basis of orthonormal vectors in v.

  -- NOTICE :
  --   Be careful with remove embedding: it does not entirely remove
  --   the empty equations resulting from the dummy slack variables.

    p : constant Poly_Sys := Remove_Embedding(ep,k);
    sp : constant Poly_Sys := Square(p,n-k);
    esp : constant Eval_Poly_Sys -- := Create(p(p'first..p'last-k));
        := Create(sp);           -- := Create(p);
   -- jm : Jaco_Mat(p'first..p'last-k,1..n) := Create(p(p'first..p'last-k));
   -- jm : Jaco_Mat(p'range,1..n) := Create(p);
    jm : Jaco_Mat(sp'range,1..n) := Create(sp);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put_line("The original system with the embedding removed : ");
    put(p);
    put_line("The system made square : ");
    put(sp);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Affine_Solution_Evaluation(esp,ejm,ls.v(1..n),b,v);
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Eval;

  procedure Test_Affine_Refiner
              ( ep : in Poly_Sys; sols : in Solution_List;
                k : in integer32 ) is

    sli : VecVec(1..k) := Slices(ep,natural32(k));
    ans : character;
    n : constant integer32 := ep'last-k;
    v,w : VecVec(1..n-k);
    b : Vector(1..n);

  begin
    put("The dimension of the ambient space is "); put(n,1);
    put_line(".");
    put("The dimension of the linear slicing space is "); put(k,1);
    put_line(".");
    Affine_Orthonormal_Basis(n,k,sli,b,v,w);
    put("Do you wish to see the original slices ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line(sli);
    end if;
    put("Do you wish to see the computed basis point ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("A basis point : "); put_line(b);
    end if;
    put("Do you wish to see the computed directions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("The directions : "); put_line(v);
    end if;
    put("Do you wish to see an orthonormal basis ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("The computed orthonormal basis :"); put_line(w);
    end if;
    put("Do you wish to continue with Newton refinement ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Eval(n,k,ep,sols,b,w);
    end if;
  end Test_Affine_Refiner;

  procedure Test_Projective_Refiner
              ( ep : in Poly_Sys; sols : in Solution_List;
                k : in integer32 ) is

    sli : VecVec(1..k) := Slices(ep,k);
    n : constant integer32 := ep'last-k;
    p : constant Poly_Sys := Remove_Embedding(ep,k);
    pp : Poly_Sys(ep'first..n);
    v,w : VecVec(1..n-k);
    b : Vector(1..n);

  begin
    put("The dimension of the ambient space is "); put(n,1);
    put_line(".");
    for i in ep'first..n loop
      pp(i) := Projective_Transformation(p(i));
    end loop;
    Symbol_Table.Replace(natural32(n)+1,Create("z0"));
    put_line("The system in projective coordinates : ");
    put(pp);
    Affine_Orthonormal_Basis(n,k,sli,b,v,w);
  end Test_Projective_Refiner;

  procedure Validate_Results
              ( file : in file_type; p : in Poly_Sys; k,n : in integer32;
                b : in Vector; v : in VecVec; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Writes the polynomial system p, the linear slices spanned by (b,v),
  --   and the solution list.  A standard root validation is applied.

    sys : Poly_Sys(1..n);
    eqs : constant VecVec := Equations(b,v);
    ind : integer32 := p'last;
    epsxa : constant double_float := 1.0E-13;
    epsfa : constant double_float := 1.0E-13;
    tolsing : constant double_float := 1.0E-8;
    nit : natural32 := 0;
    deflate : boolean := false;

  begin
    put_line(file,"The hyperplanes : "); put_line(file,eqs);
    sys(p'range) := p;
    for i in eqs'range loop
      ind := ind + 1;
      sys(ind) := Hyperplane(eqs(i).all);
    end loop;
    put_line(file,sys);
    new_line(file);
   -- put_line(file,"THE SOLUTIONS :");
   -- put(file,Length_Of(sols),Head_Of(sols).n,sols);
    Reporting_Root_Refiner
      (file,sys,sols,epsxa,epsfa,tolsing,nit,4,deflate,false);
  end Validate_Results;

  procedure Test_Affine_Sampler
              ( ep : in Poly_Sys; sols : in Solution_List;
                k : in integer32 ) is

    file : file_type;
    sli : VecVec(1..k) := Slices(ep,natural32(k));
    n : constant integer32 := ep'last-k;
    start_v,start_w,target_v,target_w : VecVec(1..n-k);
    start_b,target_b : Vector(1..n);
    p : constant Poly_Sys := Remove_Embedding(ep,k);
    sp : constant Poly_Sys := Square(p,n-k);
    eps : constant Eval_Poly_Sys := Create(sp); 
        -- := Create(p); --Create(p(p'first..p'last-k));
   -- jm : Jaco_Mat(p'first..p'last-k,1..n) := Create(p(p'first..p'last-k));
   -- jm : Jaco_Mat(p'range,1..n) := Create(p);
    jm : Jaco_Mat(sp'range,1..n) := Create(sp);
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ans : character;
    output : boolean := false;
    oc : natural32;
    gstartsols,gnewsols,newsols : Solution_List;

  begin
    new_line;
    put("Do you wish intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      output := true;
      put_line("Reading the name of the output file...");
      Read_Name_and_Create_File(file);
      new_line;
      Driver_for_Process_io(file,oc);
    end if;
    Affine_Orthonormal_Basis(n,k,sli,start_b,start_v,start_w);
    gstartsols := Decompose(sols,start_b,start_w);
    Set_Continuation_Parameter(gstartsols,Create(0.0));
    loop
      Random_Affine_Plane(n,k,target_b,target_v);
      target_w := Orthogonalize(target_v);
      Continuation_Parameters.Tune(0);
      Copy(gstartsols,gnewsols);
      if output then
        Reporting_Affine_Sampler
          (file,eps,ejm,start_b,target_b,start_w,target_w,gnewsols);
        newsols := Combine(gnewsols,target_b,target_w);
        Validate_Results(file,sp,k,n,target_b,target_w,newsols);
        --Validate_Results(file,p,k,n,target_b,target_w,newsols);
        -- (file,p(p'first..p'last-k),k,n,target_b,target_w,newsols);
      else
        Silent_Affine_Sampler
          (eps,ejm,start_b,target_b,start_w,target_w,gnewsols);
        newsols := Combine(gnewsols,target_b,target_w);
        Validate_Results(Standard_Output,sp,k,n,target_b,target_w,newsols);
        --Validate_Results(Standard_Output,p,k,n,target_b,target_w,newsols);
        -- p(p'first..p'last-k),k,n,target_b,target_w,newsols);
      end if;
      put("Do you wish more samples ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Affine_Sampler;

  procedure Main is

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : integer32;
    ans : character;

  begin
    new_line;
    put_line("Sampling with a parametric description of the linear slices.");
    new_line;
    put_line("Choose one of the following :");
    put_line("  1. affine Newton's method on a list of witness points;");
    put_line("  2. projective Newton's method on a list of witness points;");
    put_line("  3. compute new samples on a random affine plane;");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    Standard_Read_Embedding(lp,sols,natural32(dim));
    case ans is
      when '1' => Test_Affine_Refiner(lp.all,sols,dim);
      when '2' => Test_Projective_Refiner(lp.all,sols,dim);
      when '3' => Test_Affine_Sampler(lp.all,sols,dim);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sampar;
