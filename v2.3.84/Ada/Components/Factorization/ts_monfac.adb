with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Continuation_Parameters;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
--with Standard_Affine_Planes;            use Standard_Affine_Planes;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Plane_Operations;         use Standard_Plane_Operations;
with Standard_Affine_Solutions;         use Standard_Affine_Solutions;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Affine_Sampling_Machine;           use Affine_Sampling_Machine;
with Monodromy_Partitions;              use Monodromy_Partitions;

procedure ts_monfac is

-- DESCRIPTION :
--   Interactive development of the factorization of a pure dimensional
--   solution set using monodromy.

  function Product ( d : integer32; sols : Solution_List; v : Vector )
                   return Vector is

  -- DESCRIPTION :
  --   Returns the vector of all componentwise products of the solution
  --   vectors in sols with the given vector v.

    res : Vector(1..d);
    tmp : Solution_List := sols;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp).v*v;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Product;

  procedure Monodromy_Loop
              ( file : in file_type; n,k,d : in integer32; h,s : in Vector;
                eps : in Eval_Poly_Sys; ejm : Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in Solution_List;
                nb : in out natural32;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Computes one loop of the monodromy breakup algorithm.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space;
  --   k        dimension of the solution set;
  --   d        degree of the solution set;
  --   h        random vector to separate the solutions;
  --   s        products of the random vector h with the solution vectors.
  --   eps      defining equations in evaluable form;
  --   ejm      evaluable form of the Jacobi matrix;
  --   b        basis point or offset of an affine plane;
  --   v        orthonormal basis of n-k directions of an affine plane;
  --   sols     d witness points of the solution set, represented by
  --            coefficients of the linear combination of the directions,
  --            all solutions satisfy eps and the plane defined by (b,v);
  --   nb       current number of irreducible components;
  --   deco     current decomposition of the witness points.

  -- ON RETURN :
  --   nb       updated number of irreducible components;
  --   deco     updated decomposition of the witness points.

    tol : constant double_float := 1.0E-8;
    b1,b2 : Vector(b'range);
    v1,v2,w1,w2 : VecVec(v'range);
    s1 : Vector(s'range);
    gnewsols : Solution_List;
    perm : Standard_Natural_Vectors.Vector(s'range);

  begin
    Random_Affine_Plane(n,k,b1,v1); w1 := Orthogonalize(v1);
    Random_Affine_Plane(n,k,b2,v2); w2 := Orthogonalize(v2);
    Copy(sols,gnewsols);
    Set_Continuation_Parameter(gnewsols,Create(0.0));
    Silent_Affine_Sampler(eps,ejm,b,b1,v,w1,gnewsols);
    Set_Continuation_Parameter(gnewsols,Create(0.0));
    Silent_Affine_Sampler(eps,ejm,b1,b2,w1,w2,gnewsols);
    Set_Continuation_Parameter(gnewsols,Create(0.0));
    Silent_Affine_Sampler(eps,ejm,b2,b,w2,v,gnewsols);
    s1 := Product(d,gnewsols,h);
    perm := Map(s,s1,tol);
    put(file,"A permutation : "); put(file,perm); new_line(file);
    Add_Map(deco,nb,perm);
  end Monodromy_Loop;

  procedure Monodromy_Breakup
              ( file : in file_type; n,k,d : in integer32; p : in Poly_Sys;
                eps : in Eval_Poly_Sys; ejm : Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Application of monodromy to decompose a solution set in n-space
  --   of dimension k and of degree d.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        dimension of the ambient space;
  --   k        dimension of the solution set;
  --   d        degree of the solution set;
  --   p        original polynomial system, only needed for validation;
  --   eps      defining equations in evaluable form;
  --   ejm      evaluable form of the Jacobi matrix;
  --   b        basis point or offset of an affine plane;
  --   v        orthonormal basis of n-k directions of an affine plane;
  --   sols     d witness points of the solution set, represented by
  --            coefficients of the linear combination of the directions,
  --            all solutions satisfy eps and the plane defined by (b,v).

    separator : Vector(1..n-k) := Random_Vector(1,n-k);
    ips : Vector(1..d) := Product(d,sols,separator);
    nb : natural32 := natural32(d);
    deco : Standard_Natural_VecVecs.Link_to_VecVec := Init_Factors(nb);

  begin
    Continuation_Parameters.Tune(0);
    for i in 1..10 loop
      Monodromy_Loop(file,n,k,d,separator,ips,eps,ejm,b,v,sols,nb,deco);
      put_line(file,"The current decomposition :");
      Write_Factors(file,deco.all);
    end loop;
  end Monodromy_Breakup;

  procedure Test_Monodromy_Breakup
              ( file : in file_type;
                ep : in Poly_Sys; sols : in Solution_List;
                k : in integer32 ) is

    sli : VecVec(1..k) := Slices(ep,natural32(k));
    n : constant integer32 := ep'last-k;
    d : constant integer32 := integer32(Length_Of(sols));
    start_v,start_w : VecVec(1..n-k);
    start_b : Vector(1..n);
    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(k));
    eps : constant Eval_Poly_Sys := Create(p(p'first..p'last-k));
    jm : Jaco_Mat(p'first..p'last-k,1..n) := Create(p(p'first..p'last-k));
    ejm : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    gstartsols : Solution_List;

  begin
    Affine_Orthonormal_Basis(n,k,sli,start_b,start_v,start_w);
   -- gstartsols := Decompose(sols,start_b,start_w);
    gstartsols := Project(sols,start_b,start_w);
    Monodromy_Breakup(file,n,k,d,p,eps,ejm,start_b,start_w,gstartsols);
  end Test_Monodromy_Breakup;

  procedure Main is

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;

  begin
    new_line;
    put_line("Factoring a pure dimensional solution set using monodromy.");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Test_Monodromy_Breakup(file,lp.all,sols,integer32(dim));
  end Main;

begin
  Main;
end ts_monfac;
