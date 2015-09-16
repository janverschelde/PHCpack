with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Quad_Double_numbers;               use Quad_Double_numbers;
with QuadDobl_Random_Numbers;           use QuadDobl_Random_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Random_Vectors;           use QuadDobl_Random_Vectors;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;      use QuadDobl_Complex_Matrices_io;
with QuadDobl_Random_Matrices;          use QuadDobl_Random_Matrices;
with QuadDobl_Complex_Linear_Solvers;   use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Poly_Functions;   use QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_SysFun;      use QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Plane_Representations;    use QuadDobl_Plane_Representations;
with QuadDobl_Lined_Hypersurfaces;      use QuadDobl_Lined_Hypersurfaces;

package body QuadDobl_Gridded_Hypersurfaces is

-- INTERNAL DATA :

  ep : Link_to_Eval_Poly_Sys;

-- AUXILIARIES FOR STACKED GRIDS :

  function Offset ( hyp : VecVec; b,v : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns a new offset vector for the line b+t*v, given in its 
  --   explicit representations by the hyperplane equations in hyp.

    res : Vector(b'range);
    n1 : constant integer32 := b'length - 1;
    mat : Matrix(1..n1,1..n1);
    rhs : Vector(1..n1);
    ipvt : Standard_Integer_Vectors.Vector(1..n1);
    info : integer32;
    zero : constant quad_double := create(0.0);

  begin
    for i in 1..n1 loop
      for j in 1..n1 loop
        mat(i,j) := hyp(i)(j);
      end loop;
      rhs(i) := -hyp(i)(0);
    end loop;
    lufac(mat,n1,ipvt,info);
    lusolve(mat,n1,ipvt,rhs);
    res(rhs'range) := rhs;
    res(res'last) := Create(zero);
    return res;
  end Offset;

  function Deep_Copy ( v : VecVec ) return VecVec is

  -- DESCRIPTION :
  --   Returns a deep copy of the content of v.

    res : VecVec(v'range);

  begin
    for i in res'range loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'(v(i).all);
    end loop;
    return res;
  end Deep_Copy;

-- FORMAT CONVERSIONS :

  function Create ( b,v : Vector ) return VecVec is
  begin
    return Equations1(b,v);
  end Create;

  function Create ( b,v : Vector; t : Complex_Number ) return Solution is

    res : Solution(b'length);
    one : constant quad_double := create(1.0);

  begin
    res.t := Create(one);
    res.m := 1;
    for i in b'range loop
      res.v(i) := b(i) + t*v(i);
    end loop;
    res.err := create(0.0);
    res.rco := create(1.0);
    res.res := create(0.0);
    return res;
  end Create;

  function Create ( b,v : Vector; t : Complex_Number )
                  return QuadDobl_Sample is

    sol : constant Solution(b'length) := Create(b,v,t);
    hyp : constant VecVec := Create(b,v);

  begin
    return Create(sol,hyp);
  end Create;

  function Create ( b,v,t : Vector ) return QuadDobl_Sample_List is

    res,res_last : QuadDobl_Sample_List;

  begin
    for i in t'range loop
      Append(res,res_last,Create(b,v,t(i)));
    end loop;
    return res;
  end Create;

-- THE SAMPLING MACHINE :

  procedure Initialize ( p : in Poly ) is

    n : integer32;
    dp : Poly;

  begin
    if p /= Null_Poly then
      n := integer32(Number_of_Unknowns(p));
      ep := new Eval_Poly_Sys(0..n);
      ep(0) := Create(p);
      for i in 1..n loop
        dp := Diff(p,i);
        ep(i) := Create(dp);
        Clear(dp);
      end loop;
    end if;
  end Initialize;

  function Sample ( b,v,t : Vector ) return QuadDobl_Sample_List is

    b1 : constant Vector(b'range) := Random_Vector(b'first,b'last);
    v1 : constant Vector(v'range) := Random_Vector(v'first,v'last);
    t1 : Vector(t'range) := t;

  begin
    Silent_Hypersurface_Sampler(ep.all,b,v,b1,v1,t1);
    return Create(b1,v1,t1);
  end Sample;

  function Sample ( file : file_type; output : boolean;
                    b,v,t : Vector ) return QuadDobl_Sample_List is

    b1 : constant Vector(b'range) := Random_Vector(b'first,b'last);
    v1 : constant Vector(v'range) := Random_Vector(v'first,v'last);
    t1 : Vector(t'range) := t;

  begin
    Reporting_Hypersurface_Sampler(file,ep.all,b,v,b1,v1,output,t1);
    return Create(b1,v1,t1);
  end Sample;

  function Parallel_Sample ( b,v,t : Vector ) return QuadDobl_Sample_List is

    b1 : constant Vector(b'range) := Random_Vector(b'first,b'last);
    t1 : Vector(t'range) := t;

  begin
    Silent_Hypersurface_Sampler(ep.all,b,v,b1,v,t1);
    return Create(b1,v,t1);
  end Parallel_Sample;

  function Parallel_Sample ( file : file_type; output : boolean;
                             b,v,t : Vector ) return QuadDobl_Sample_List is

    b1 : constant Vector(b'range) := Random_Vector(b'first,b'last);
    t1 : Vector(t'range) := t;

  begin
    Reporting_Hypersurface_Sampler(file,ep.all,b,v,b1,v,output,t1);
    return Create(b1,v,t1);
  end Parallel_Sample;

  function Parallel_Sample1 ( b,v,t : Vector ) return QuadDobl_Sample_List is

    sli : VecVec(b'first..b'last-1) := Equations1(b,v); 
    b1 : Vector(b'range);
    t1 : Vector(t'range) := t;

  begin
    sli(b'first)(0) := Random1;
    b1 := OffSet(sli,b,v);
    Clear(sli);
    Silent_Hypersurface_Sampler(ep.all,b,v,b1,v,t1);
    return Create(b1,v,t1);
  end Parallel_Sample1;

  function Parallel_Sample1 ( file : file_type; output : boolean;
                              b,v,t : Vector ) return QuadDobl_Sample_List is

    sli : VecVec(b'first..b'last-1) := Equations1(b,v); 
    b1 : Vector(b'range);
    t1 : Vector(t'range) := t;

  begin
    sli(b'first)(0) := Random1;
    b1 := OffSet(sli,b,v);
    Clear(sli);
    Reporting_Hypersurface_Sampler(file,ep.all,b,v,b1,v,output,t1);
    return Create(b1,v,t1);
  end Parallel_Sample1;

  function Sample ( b,v,t : Vector; k : integer32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..k);

  begin
    res(0) := Create(b,v,t);
    for i in 1..k loop
      res(i) := Sample(b,v,t);
    end loop;
    return res;
  end Sample;

  function Sample ( file : file_type; output : boolean;
                    b,v,t : Vector; k : integer32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..k);

  begin
    res(0) := Create(b,v,t);
    for i in 1..k loop
      res(i) := Sample(file,output,b,v,t);
    end loop;
    return res;
  end Sample;

  function Parallel_Sample ( b,v,t : Vector; k : integer32 )
                           return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..k);

  begin
    res(0) := Create(b,v,t);
    for i in 1..k loop
      res(i) := Parallel_Sample(b,v,t);
    end loop;
    return res;
  end Parallel_Sample;

  function Parallel_Sample ( file : file_type; output : boolean;
                             b,v,t : Vector; k : integer32 )
                           return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..k);

  begin
    res(0) := Create(b,v,t);
    for i in 1..k loop
      res(i) := Parallel_Sample(file,output,b,v,t);
    end loop;
    return res;
  end Parallel_Sample;

  function Parallel_Sample1 ( b,v,t : Vector; k : integer32 )
                            return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..k);

  begin
    res(0) := Create(b,v,t);
    for i in 1..k loop
      res(i) := Parallel_Sample1(b,v,t);
    end loop;
    return res;
  end Parallel_Sample1;

  function Parallel_Sample1 ( file : file_type; output : boolean;
                              b,v,t : Vector; k : integer32 )
                            return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..k);

  begin
    res(0) := Create(b,v,t);
    for i in 1..k loop
      res(i) := Parallel_Sample1(file,output,b,v,t);
    end loop;
    return res;
  end Parallel_Sample1;

  function Sample ( k,n,d : integer32; b,v,t : Vector;
                    sli : VecVec; rancff : Matrix ) 
                  return Stacked_Sample_Grid is

  -- NOTE :
  --   Recursive version of the creator; k is the index to the varying slice.
  --   See Standard_Stacked_Sample_Grids.Create.

    res : Stacked_Sample_Grid(k,d);
    len : integer32;

  begin
    res.hyp := sli(1..k);
    res.n := natural32(n);
    if k = 1 then
      len := t'length;
      res.g := new Array_of_QuadDobl_Sample_Lists(0..len);
      res.g(0) := Create(b,v,t);
      for i in 1..len loop
        declare
          newsli : VecVec(sli'range) := Deep_Copy(sli);
          b1 : Vector(b'range);
          t1 : Vector(t'range) := t;
        begin
          newsli(k)(0) := rancff(k,i);
          b1 := Offset(newsli,b,v);
          Clear(newsli);
          Silent_Hypersurface_Sampler(ep.all,b,v,b1,v,t1);
          res.g(i) := Create(b1,v,t1);
        end;
      end loop;
    else
      for i in reverse 1..d loop
        declare
          newsli : VecVec(sli'range);
          gridai : Stacked_Sample_Grid(k-1,d);
        begin
          if i = d then
            gridai := Sample(k-1,n,d,b,v,t,sli,rancff);
          else
            declare
              b1 : Vector(b'range);
              t1 : Vector(t'range) := t;
            begin
              newsli := Deep_Copy(sli);
              newsli(k)(0) := rancff(k,d-i);
              b1 := OffSet(newsli,b,v);
              Silent_Hypersurface_Sampler(ep.all,b,v,b1,v,t1);
              gridai := Sample(k-1,n,d,b1,v,t1,newsli,rancff);
            end;
          end if;
          res.a(i) := new Stacked_Sample_Grid'(gridai);
        end;
      end loop;
      res.pts(0) := -res.hyp(k)(0);
      for i in 1..d loop
        res.pts(i) := -rancff(k,i);
      end loop;
      declare
        newsli : constant VecVec(sli'range) := Deep_Copy(sli);
        b1 : Vector(b'range);
        t1 : Vector(t'first..t'first) := t(t'first..t'first);
      begin
        newsli(k)(0) := rancff(k,d);
        b1 := Offset(newsli,b,v);
        Silent_Hypersurface_Sampler(ep.all,b,v,b1,v,t1);
        res.spt := Create(b1,v,t1(t1'first));
      end;
    end if;
    return res;
  end Sample;

  function Full_Sample ( k,n,d : integer32; b,v,t : Vector;
                         sli : VecVec; rancff : Matrix )
                       return Stacked_Sample_Grid is

  -- NOTE :       
  --   Recursive version of the sampler for the full grid;
  --   k is the index to the varying slice.
  --   See Standard_Stacked_Sample_Grids.Create_Full.

    res : Stacked_Sample_Grid(sli'last-k+1,d);
    len : integer32;

  begin
    for i in k..sli'last loop
      res.hyp(i-k+1) := sli(i);
    end loop;
    res.n := natural32(n);
    if k = sli'last then
      len := t'length;
      res.g := new Array_of_QuadDobl_Sample_Lists(0..len);
      res.g(0) := Create(b,v,t);
      for i in 1..len loop
        declare
          newsli : VecVec(sli'range) := Deep_Copy(sli);
          b1 : Vector(b'range);
          t1 : Vector(t'range) := t;
        begin
          newsli(k)(0) := rancff(k,i);
          b1 := Offset(newsli,b,v);
          Clear(newsli);
          Silent_Hypersurface_Sampler(ep.all,b,v,b1,v,t1);
          res.g(i) := Create(b1,v,t1);
        end;
      end loop;
    else
      for i in 0..d loop
        declare
          newsli : VecVec(sli'range);
          gridai : Stacked_Sample_Grid(sli'last-k,d);
        begin
          if i = 0 then
            gridai := Full_Sample(k+1,n,d,b,v,t,sli,rancff);
          else
            declare
              b1 : Vector(b'range);
              t1 : Vector(t'range) := t;
            begin
              newsli := Deep_Copy(sli);
              newsli(k)(0) := rancff(k,i);
              b1 := OffSet(newsli,b,v);
              Silent_Hypersurface_Sampler(ep.all,b,v,b1,v,t1);
              gridai := Full_Sample(k+1,n,d,b1,v,t1,newsli,rancff);
            end;
          end if;
          res.a(i) := new Stacked_Sample_Grid'(gridai);
        end;
      end loop;
    end if;
    res.pts(0) := -sli(k)(0);
    for i in 1..d loop
      res.pts(i) := -rancff(k,i);
    end loop;
    return res;
  end Full_Sample;

  function Sample ( b,v,t : Vector ) return Stacked_Sample_Grid is

    n : constant integer32 := b'length;
    k : constant integer32 := n-1;
    d : constant integer32 := t'length;
    rancff : constant Matrix(1..k,1..d)
           := Random_Matrix(natural32(k),natural32(d));
    sli : constant VecVec := Equations1(b,v); 

  begin
    return Sample(k,n,d,b,v,t,sli,rancff);
  end Sample;

  function Sample ( file : file_type;
                    b,v,t : Vector ) return Stacked_Sample_Grid is

    n : constant integer32 := b'length;
    k : constant integer32 := n-1;
    d : constant integer32 := t'length;
    rancff : constant Matrix(1..k,1..d)
           := Random_Matrix(natural32(k),natural32(d));
    sli : constant VecVec := Equations1(b,v);

  begin
    put_line(file,"The matrix of random constant coefficients : ");
    put(file,rancff,3);
    return Sample(k,n,d,b,v,t,sli,rancff);
  end Sample;

  function Full_Sample ( b,v,t : Vector ) return Stacked_Sample_Grid is

    n : constant integer32 := b'length;
    k : constant integer32 := n-1;
    d : constant integer32 := t'length;
    rancff : constant Matrix(1..k,1..d)
           := Random_Matrix(natural32(k),natural32(d));
    sli : constant VecVec := Equations1(b,v);

  begin
    return Full_Sample(1,n,d,b,v,t,sli,rancff);
  end Full_Sample;

  function Full_Sample ( file : file_type;
                         b,v,t : Vector ) return Stacked_Sample_Grid is

    n : constant integer32 := b'length;
    k : constant integer32 := n-1;
    d : constant integer32 := t'length;
    rancff : constant Matrix(1..k,1..d)
           := Random_Matrix(natural32(k),natural32(d));
    sli : constant VecVec := Equations1(b,v);

  begin
    put_line(file,"The matrix of random constant coefficients : ");
    put(file,rancff,3);
    return Full_Sample(1,n,d,b,v,t,sli,rancff);
  end Full_Sample;

  procedure Clear is
  begin
    Clear(ep);
  end Clear;

end QuadDobl_Gridded_Hypersurfaces;
