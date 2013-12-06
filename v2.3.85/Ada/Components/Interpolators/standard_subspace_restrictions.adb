with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Planes_and_Polynomials;             use Planes_and_Polynomials;
with Standard_Linear_Spaces;             use Standard_Linear_Spaces;

package body Standard_Subspace_Restrictions is

  procedure Container_Dimension
              ( file : in file_type; k,n : in natural32;
                samples : in Standard_Complex_VecVecs.VecVec;
                tol : in double_float; mat,vec : out Matrix;
                dim : out natural32 ) is
  begin
    for i in 1..integer32(k)+1 loop
      for j in 1..integer32(n) loop
        mat(i,j) := samples(i)(j);
      end loop;
    end loop;
    Rank(integer32(k),integer32(n),mat,tol,vec,dim);
   -- put_line(file,"Structure of the matrix : ");
   -- for i in 1..k loop
   --   for j in 1..n loop
   --     if AbsVal(vec(i,j)) > tol
   --      then put(file,"*");
   --      else put(file,"0");
   --     end if;
   --   end loop;
   --   new_line(file);
   -- end loop;
   -- put(file,"The dimension of the space : ");
   -- put(file,dim,1); new_line(file);
  end Container_Dimension;

  procedure Container_Subspace
               ( file : in file_type; k,n,level,dim : in natural32;
                 p : in Poly_Sys; mat,vec : in Matrix; tol : in double_float;
                 samples : in Standard_Complex_VecVecs.VecVec;
                 restsamp : out Standard_Complex_VecVecs.VecVec;
                 lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                 kerpols,restp : out Link_to_Poly_Sys ) is

    piv : constant Standard_Integer_Vectors.Vector
        := Pivots(integer32(k),integer32(n),vec,dim,tol);
    ker : constant Standard_Complex_VecVecs.VecVec(1..integer32(n-dim))
        := Kernel(mat,vec,dim,piv,tol);
    kerpolsys : Poly_Sys(ker'range);
   -- subpols : Poly_Sys(ker'range);
    nvar : constant integer32 := integer32(dim+level);
   -- eva : Complex_Number;
    restsys : Poly_Sys(p'range);

  begin
   -- put_line(file,"The hyperplanes : ");
    for i in kerpolsys'range loop
      kerpolsys(i) := Hyperplane(ker(i).all,tol);
    end loop;
   -- put_line(file,kerpolsys);
   -- for i in ker'range loop
   --   put(file,"Evaluate points "); put(file,i,1);
   --   put_line(file,"-th hyperplane : ");
   --   for j in 1..k+1 loop
   --     eva := ker(i)(0);
   --     for k in mat'range(2) loop
   --       eva := eva + ker(i)(k)*mat(j,k);
   --     end loop;
   --     put(file,"eval : "); put(file,eva); new_line(file);
   --   end loop;
   -- end loop;
   -- subpols := Substituting_Polynomials(ker'last,nvar,piv,ker,tol);
   -- put_line(file,"The substituting polynomials : ");
   -- put_line(file,subpols);
    restsys := Restrict_to_Linear_Space(p,integer32(level),piv,ker,tol);
   -- put_line(file,"The system restricted to the subspace : ");
   -- put_line(file,restsys);
   -- put(file,"The number of variables of every polynomial : ");
   -- for i in restsys'range loop
   --   put(file,Number_of_Unknowns(restsys(i)),1); put(file," ");
   -- end loop;
   -- new_line(file);
    kerpols := new Poly_Sys'(kerpolsys);
    restp := new Poly_Sys'(restsys);
   -- put_line(file,"Evaluating the solutions restricted to the space : ");
    for i in samples'range loop
      exit when (samples(i) = null);  -- this is again a patch ...
      restsamp(i) := new Vector'
        (Remove_Variables(samples(i).all,integer32(level),nvar,piv));
     -- put(file,"restricted solution "); put(file,i,1); put_line(file," : ");
     -- put_line(file,restsamp(i).all);
     -- put(file,"evaluation "); put(file,i,1); put_line(file," : ");
     -- for j in restsys'first..restsys'last-level loop  -- other slices !
     --   put(file,Eval(restsys(j),restsamp(i).all)); new_line(file);
     -- end loop;
    end loop;
    lpiv := new Standard_Integer_Vectors.Vector'(piv);
  end Container_Subspace;

  procedure Subspace_Restriction
                ( file : in file_type; embsys : in Poly_Sys;
                  ind,k,n,level : in natural32;
                  samples : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  restsamp : in out Standard_Complex_VecVecs.Array_of_VecVecs;
                  lpiv : out Standard_Integer_Vectors.Link_to_Vector;
                  kerpols,restembsys : out Link_to_Poly_Sys;
                  dim : out natural32 ) is

    tol : constant double_float := 10.0**(-8);
    mat : Matrix(1..integer32(n)+1,1..integer32(n));
    vec : Matrix(1..integer32(n),1..integer32(n));

  begin
    Container_Dimension(file,k,n,samples(integer32(ind)).all,tol,mat,vec,dim);
    if dim < k then
      Container_Subspace
        (file,k,n,level,dim,embsys,mat,vec,tol,samples(integer32(ind)).all,
         restsamp(integer32(ind)).all,lpiv,kerpols,restembsys);
    end if;
  end Subspace_Restriction;

  function Collapse_Equations
             ( p : Poly_Sys; dim,level : natural32 ) return Poly_Sys is

    n : constant integer32 := integer32(dim+level);
    res : Poly_Sys(1..n);
    acc : Poly;

  begin
    for i in 1..integer32(dim) loop               -- make random combinations
      if not (p(p'first) = Null_Poly)
       then Copy(p(p'first),res(i));
      end if;
      for j in p'first+1..p'last-integer32(level) loop
        if not (p(j) = Null_Poly) then
          acc := Random1*p(j);
          Add(res(i),acc);
          Clear(acc);
        end if;
      end loop;
    end loop;
    for i in 1..integer32(level) loop             -- copy last equations
      Copy(p(p'last+1-i),res(n+1-i));
    end loop;
    return res;
  end Collapse_Equations;

end Standard_Subspace_Restrictions;
