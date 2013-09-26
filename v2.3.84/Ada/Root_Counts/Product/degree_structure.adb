with unchecked_deallocation;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
--with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
--with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
--with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
--with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Degrees_in_Sets_of_Unknowns;        use Degrees_in_Sets_of_Unknowns;
with Generate_Unions;

package body Degree_Structure is

-- DECLARATIONS :

  type zd ( m : natural32; im : integer32 ) is record
    z : Partition(1..m); -- requirement: im = m
    d : Standard_Natural_Vectors.Vector(1..im);
  end record;

  type Link_To_zd is access zd;

  procedure free is new unchecked_deallocation(zd,Link_To_zd);

  type dgst is array ( natural32 range <> ) of Link_To_zd;

  type Link_To_dgst is access dgst;

  procedure free is new unchecked_deallocation(dgst,Link_To_dgst);

-- INTERNAL DATA :

  ds : Link_To_dgst;

-- CREATORS :

  procedure Find_Partition
              ( p : in Poly; 
                z : in out Partition; m : in out natural32; 
                dg : in out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   This routine finds an good partition for the polynomial p 

  -- ON ENTRY :
  --   p        a polynomial.

  -- ON RETURN :
  --   z        a partition of the set of unknowns of p;
  --   m        the number of sets in the partition z;
  --   dg       the degrees of the polynomial p in the sets of z,
  --            dg(i) = Degree(p,z(i)).

    n : constant natural32 := Number_Of_Unknowns(p);
    di : integer32;
    added : boolean;

  begin
    for i in 1..n loop
      di := Degree(p,integer32(i));
      if di > 0 then
        added := false;
        for j in 1..m loop
          if di = dg(integer32(j)) then
            Add(z(j),i);
            if Degree(p,z(j)) = dg(integer32(j))
             then added := true;
             else Remove(z(j),i);
            end if;
          end if;
          exit when added;
        end loop;
        if not added then
          m := m + 1;
          Add(z(m),i);
          dg(integer32(m)) := di;
        end if;
      end if;
    end loop;
  end Find_Partition;

  procedure Create ( p : in Poly_Sys ) is

    n : constant natural32 := natural32(p'length);
    z : Partition(1..n);
    d : Standard_Integer_Vectors.Vector(1..integer32(n))
      := (1..integer32(n) => 0);
    m : natural32;

  begin
    if ds /= null
     then Clear;
    end if;
    ds := new dgst(1..n);
    for i in 1..integer32(n) loop
      m := 0;
      Create(z,n);
      Find_Partition(p(i),z,m,d);
      ds(natural32(i)) := new zd(m,integer32(m));
      for j in 1..m loop
        ds(natural32(i)).z(j) := Create(z(j));
        ds(natural32(i)).d(integer32(j)) := natural32(d(integer32(j)));
      end loop;
      Clear(z);
    end loop;
  end Create;

  procedure Put ( p : in Poly_Sys;
                  i,m : in natural32; z : in Partition ) is

    n : constant natural32 := natural32(p'length);

  begin
    if ds = null
     then ds := new dgst(1..n);
    end if;
    ds(i) := new zd(m,integer32(m));
    for j in 1..m loop
      ds(i).z(j) := Create(z(j));
      ds(i).d(integer32(j)) := natural32(Degree(p(integer32(i)),z(j)));
    end loop;
  end Put;

-- SELECTORS :

  function Empty return boolean is
  begin
    return (ds = null);
  end Empty;

  function Get ( i : natural32 ) return natural32 is
  begin
    return ds(i).m;
  end Get;

  procedure Get ( i : in natural32; z : in out Partition;
                  d : out Standard_Natural_Vectors.Vector ) is
  begin
    for j in 1..ds(i).m loop
      z(j) := Create(ds(i).z(j));
      d(integer32(j)) := natural32(ds(i).d(integer32(j)));
    end loop;
  end Get;

-- COMPUTING THE GENERALIZED BEZOUT NUMBER :

 -- function Matrix_Criterion ( z : Partition ) return boolean is

  -- DESCRIPTION : 
  --   This is the Matrix criterion for testing if
  --   a product is admissible or not.

  --  n : constant integer32 := integer32(z'last - z'first + 1);
  --  mat : Matrix(1..n,1..n);
  --  ipvt : Standard_Integer_Vectors.Vector(1..n);
  --  eps : constant double_float := 10.0**(-10);
  --  r : Complex_Number;
  --  rcond : double_float;

  --begin
  --  for i in 1..n loop
  --    r := Create(double_float(i+1));
  --    for j in 1..n loop
  --      if Is_In(z(natural32(i)),natural32(j))
  --       then mat(i,j) := r; r := r*r;
  --       else mat(i,j) := Create(0.0);
  --      end if;
  --    end loop;
  --  end loop;
  --  lufco(mat,n,ipvt,rcond);
  --  return (abs(rcond) > eps);
  --exception
  --  when others => return false;
  --end Matrix_Criterion;
  
  function Admissible ( z : Partition; n : natural32 ) return boolean is

    temp : Partition(1..n);
    admis : boolean := true;

  begin
    temp(1) := Create(z(1));
    for i in 2..(n-1) loop
      temp(i) := Create(z(i));
      admis := Admissible(temp,i,z(i+1));
      exit when not admis;
    end loop;
    Clear(temp);
    return admis;
  end Admissible;

  function Admissible ( z : Partition; n : natural32; s : Set )
                      return boolean is
  begin
    for k in 1..n loop
      if not Admissible(z,k,n,s)
       then return false;
      end if;
    end loop;
    return true;
  end Admissible;

  function Admissible ( z : Partition; k,n : natural32; s : Set )
                      return boolean is

    type arr is array ( natural32 range <> ) of boolean;
    admis : boolean := true;

    procedure check ( a : in arr; continue : out boolean ) is

      u : Set := Create(s);

    begin
      for i in a'range loop
        if a(i)
         then Union(u,z(i));
        end if;
      end loop;
      admis := ( Extent_Of(u) >= k+1 );
      continue := admis;
      Clear(u);
    end check;

    procedure gen is new Generate_Unions(arr,check);

  begin
    gen(k,1,n);
    return admis;
  end Admissible;

  procedure Compute ( i,n,sum : in natural32; res : in out natural32;
                      z : in out Partition ) is
  begin
    if i > n then
      res := res + sum;
    else -- Pick out a set and check if it is allowed :
      for j in 1..ds(i).m loop
        if ds(i).d(integer32(j)) /= 0
            and then Admissible(z,i-1,ds(i).z(j)) then
          z(i) := Create(ds(i).z(j));
          Compute(i+1,n,sum*ds(i).d(integer32(j)),res,z);
          Clear(z(i));
        end if;
      end loop;
    end if;
  end Compute;

  function Generalized_Bezout_Number return natural32 is

    res : natural32 := 0;
    n : constant natural32 := ds'length;
    z : Partition(1..n);

  begin
    Compute(1,n,1,res,z);
    return res;
  end Generalized_Bezout_Number;

  function Generalized_Bezout_Number ( p : in Poly_Sys ) return natural32 is
  begin
    Create(p);
    return Generalized_Bezout_Number;
  end Generalized_Bezout_Number;

-- DESTRUCTOR :

  procedure Clear is
  begin
    if ds /= null
     then for i in ds'range loop
            free(ds(i));
          end loop;
          free(ds);
    end if;
  end Clear;

end Degree_Structure;
