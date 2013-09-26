with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Matrices_io;       use Multprec_Complex_Matrices_io;
with Multprec_to_Standard_Convertors;    use Multprec_to_Standard_Convertors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;    use Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Embed_Polynomials;         use Standard_Embed_Polynomials;
with Multprec_Linear_Spaces;             use Multprec_Linear_Spaces;
with Planes_and_Polynomials;             use Planes_and_Polynomials; 

package body Multprec_Subspace_Restrictions is

  procedure Container_Dimension
              ( file : in file_type; k,n,size : in natural32;
                samples : in Multprec_Complex_VecVecs.VecVec;
                tol : in double_float; mat,vec : in out Matrix;
                dim : out natural32 ) is

   -- tmp : Floating_Number;

  begin
    for i in 1..integer32(k)+1 loop
      for j in 1..integer32(n) loop
        Copy(samples(i)(j),mat(i,j));
      end loop;
    end loop;
   -- for i in 1..k+1 loop
   --   put(file,"Sample no "); put(file,i,1); put_line(file," :");
   --   for j in 1..n loop
   --     put(file,j,1); put(file," : ");
   --     put(file,samples(i)(j)); new_line(file);
   --   end loop;
   -- end loop;
   -- put_line(file,"The matrix : "); put(file,mat);
    Rank(integer32(k),integer32(n),size,mat,tol,vec,dim);
   -- put_line(file,"The triangulated matrix : ");
   -- put(file,vec,3);
   -- put_line(file,"Structure of the matrix : ");
   -- for i in 1..k loop
   --   for j in 1..n loop
   --     tmp := AbsVal(vec(i,j));
   --     if tmp > tol
   --      then put(file,"*");
   --      else put(file,"0");
   --     end if;
   --     Clear(tmp);
   --   end loop;
   --   new_line(file);
   -- end loop;
   -- put(file,"The dimension of the space : ");
   -- put(file,dim,1); new_line(file);
  end Container_Dimension;

  procedure Container_Subspace
               ( file : in file_type; k,n,level,dim : in natural32;
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 embp : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mat,vec : in Matrix; tol : in double_float;
                 samples : in Multprec_Complex_VecVecs.VecVec;
                 restsamp : out Multprec_Complex_VecVecs.VecVec;
                 lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                 kerpols,restp
                   : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                 restembp
                   : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    piv : constant Standard_Integer_Vectors.Vector
        := Pivots(integer32(k),integer32(n),vec,dim,tol);
    ker : Multprec_Complex_VecVecs.VecVec(1..integer32(n-dim))
        := Kernel(mat,vec,dim,piv,tol);
    kerpolsys : Multprec_Complex_Poly_Systems.Poly_Sys(ker'range);
   -- subpols : Multprec_Complex_Poly_Systems.Poly_Sys(ker'range);
    nvar : constant integer32 := integer32(dim+level);
    acc,eva : Complex_Number;
    restsys : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);
    mpembp,restembsys : Multprec_Complex_Poly_Systems.Poly_Sys(embp'range);

  begin
   -- put_line(file,"The hyperplanes : ");
    for i in kerpolsys'range loop
      kerpolsys(i) := Hyperplane(ker(i).all,tol);
    end loop;
   -- Multprec_Complex_VecVecs.Clear(ker); --> raises storage error!
   -- put_line(file,kerpolsys);
   -- for i in ker'range loop
   --   put(file,"Evaluate points "); put(file,i,1);
   --   put_line(file,"-th hyperplane : ");
   --   for j in 1..k+1 loop
   --     Copy(ker(i)(0),eva);
   --     for k in mat'range(2) loop
   --       acc := ker(i)(k)*mat(j,k);
   --       Add(eva,acc);
   --       Clear(acc);
   --     end loop;
   --     put(file,"eval : "); put(file,eva); new_line(file);
   --     Clear(eva);
   --   end loop;
   -- end loop;
   -- subpols := Substituting_Polynomials(ker'last,nvar,piv,ker,tol);
   -- put_line(file,"The substituting polynomials : ");
   -- for i in subpols'range loop
   --   put(file,"sp("); put(file,i,1); put(file,") : ");
   --   put_line(file,subpols(i));
   -- end loop;
   -- restsys := Restrict_to_Linear_Space(p,level,piv,ker,tol);
   -- WARNING : level = 0 for original system
    restsys := Restrict_to_Linear_Space(p,0,piv,ker,tol);
   -- put_line(file,"The original system restricted to the subspace : ");
   -- put_line(file,restsys);
    kerpols := new Multprec_Complex_Poly_Systems.Poly_Sys'(kerpolsys);
    restp := new Multprec_Complex_Poly_Systems.Poly_Sys'(restsys);
   -- put_line(file,"Evaluating the solutions restricted to the subspace : ");
    for i in samples'range loop
      exit when (samples(i) = null);      -- this is a patch...
      restsamp(i) := new Multprec_Complex_Vectors.Vector'
        (Remove_Variables(samples(i).all,integer32(level),nvar,piv));
   --   put(file,"restricted solution "); put(file,i,1); put_line(file," : ");
   --   put_line(file,restsamp(i).all);
   --   put(file,"evaluation "); put(file,i,1); put_line(file," : ");
   --   for j in restsys'range loop
   --     eva := Eval(restsys(j),restsamp(i).all);
   --     put(file,eva); new_line(file);
   --     Clear(eva);
   --   end loop;
    end loop;
    mpembp := Convert(embp);
    restembsys := Restrict_to_Linear_Space
                    (mpembp,integer32(level),piv,ker,tol);
   -- put_line(file,"The embedded system restricted to the subspace : ");
   -- put_line(file,restembsys);
    restembp := new Standard_Complex_Poly_Systems.Poly_Sys'
      (Convert(restembsys));
   -- put_line(file,"The converted restricted embedded system : ");
   -- put_line(file,restembp.all);
    for i in restsamp'range loop
      exit when (restsamp(i) = null);
     -- put(file,"evaluation of restricted sample ");
     -- put(file,i,1); put_line(file," : ");
     -- for j in restembsys'first..restembsys'last-level loop
     --   eva := Eval(restembsys(j),restsamp(i).all);
     --   put(file,eva); new_line(file);
     --   Clear(eva);
     -- end loop;
    end loop;
    Multprec_Complex_Poly_Systems.Clear(mpembp);
    Multprec_Complex_Poly_Systems.Clear(restembsys);
    lpiv := new Standard_Integer_Vectors.Vector'(piv);
  end Container_Subspace;

  procedure Subspace_Restriction
                ( file : in file_type;
                  orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                  ind,k,n,level,size : in natural32;
                  samples : in Multprec_Complex_VecVecs.Array_of_VecVecs;
                  restsamp : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
                  lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                  kerpols,restorgsys
                    : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                  restembsys
                    : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                  dim : out natural32 ) is

    tol : constant double_float := 10.0**(-8);
    mat : Matrix(1..integer32(n)+1,1..integer32(n));
    vec : Matrix(1..integer32(n),1..integer32(n));

  begin
    Container_Dimension(file,k,n,size,samples(integer32(ind)).all,
                        tol,mat,vec,dim);
    if dim < k then
      Container_Subspace
        (file,k,n,level,dim,orgsys,embsys,mat,vec,tol,
         samples(integer32(ind)).all,
         restsamp(integer32(ind)).all,lpiv,kerpols,restorgsys,restembsys);
    end if;
    Clear(mat); Clear(vec);
  end Subspace_Restriction;

  function Collapse_Equations
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys; dim : natural32 )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    use Multprec_Complex_Polynomials;

    res : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(dim));
    ran : Complex_Number;
    acc : Poly;

  begin
    for i in 1..integer32(dim) loop
      if not (p(p'first) = Null_Poly)
       then Copy(p(p'first),res(i));
      end if;
      for j in p'first+1..p'last loop
        ran := Create(Random1);
        acc := ran*p(j);
        Add(res(integer32(i)),acc);
        Clear(acc); 
        Clear(ran);
      end loop;
    end loop;
    return res;
  end Collapse_Equations;

  procedure Add_Random_Terms
               ( p : in out Standard_Complex_Polynomials.Poly;
                 dim,level : in natural32;
                 ep : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Adds random terms of ep that have added z-variables to p.

    use Standard_Complex_Polynomials;

    procedure Add_Term ( t : in Term; continue : out boolean ) is

      to_add : boolean := true;

    begin
      for i in 1..dim loop
        if t.dg(integer32(i)) /= 0
         then to_add := false;               -- term has an original variable
        end if;
        exit when not to_add;
      end loop;
      if to_add then
        to_add := false;
        for i in dim+1..dim+level loop
          if t.dg(integer32(i)) /= 0
           then to_add := true;          -- term has a z-variable
          end if;
          exit when to_add;
        end loop;
        if to_add
         then Add(p,t);
        end if;
      end if;
      continue := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    Add_Terms(ep);
  end Add_Random_Terms;

  function Embed_Collapsed_Equations
             ( restorgsys : Multprec_Complex_Poly_Systems.Poly_Sys;
               restembsys : Standard_Complex_Poly_Systems.Poly_Sys;
               dim,level : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;

    n : constant integer32 := integer32(dim+level);
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..n);
    acc : Poly;

  begin
    for i in restorgsys'range loop
      acc := Convert(restorgsys(i));
      res(i) := Add_Variables(acc,level);
      Add_Random_Terms(res(i),dim,level,restembsys(i));
      Clear(acc);
    end loop;
    for i in 1..integer32(level) loop                 -- copy last equations
      Copy(restembsys(restembsys'last+1-i),res(n+1-i));
    end loop;
    return res;
  end Embed_Collapsed_Equations;

end Multprec_Subspace_Restrictions;
