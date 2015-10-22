with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials; 
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;
with Evaluated_Minors;                   use Evaluated_Minors;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Numeric_Minor_Equations;            use Numeric_Minor_Equations;
with Plane_Representations;              use Plane_Representations;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;

package body Determinantal_Systems is

-- AUXILIARY TO EVALUATION OF MAXIMAL MINORS :

  function Number_of_Maximal_Minors
             ( nrows,ncols : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of maximal minors of a matrix with nrows and ncols.

  -- REQUIRED : nrows > ncols.

  begin
    if ncols = 1
     then return nrows;
     else return Number_of_Maximal_Minors(nrows-1,ncols-1)*nrows/ncols;
    end if;
  end Number_of_Maximal_Minors;

-- LOCALIZATION MAPS :

  function Localize ( locmap : Standard_Natural_Matrices.Matrix;
                      t : Term ) return Term is

  -- DESCRIPTION :
  --   Applies the localization map to the term, eliminating those xij's
  --   xij for which the corresponding entry in locmap is either 0 or 1.

  -- NOTE : 
  --   This localization assumes that t.dg(k) = 0 with k for which the
  --   corresponding (i,j) with locmap(i,j) = 0.

    res : Term;
    ndg : Standard_Natural_Vectors.Vector(t.dg'range);
    cnt : integer32 := t.dg'first-1;
    ind : integer32 := cnt;

  begin
    for i in locmap'range(1) loop       -- lexicographic order of variables
      for j in locmap'range(2) loop
        ind := ind+1;
        if locmap(i,j) = 2 then
          cnt := cnt + 1;
          ndg(cnt) := t.dg(ind);
        end if;
      end loop;
    end loop;
    for i in ind+1..t.dg'last loop      -- leave the lifting !
      cnt := cnt+1;
      ndg(cnt) := t.dg(i);
    end loop;
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector'(ndg(1..cnt));
    return res;
  end Localize;

  function Localize ( locmap : Standard_Natural_Matrices.Matrix;
                      p : Poly ) return Poly is

  -- DESCRIPTION :
  --   Applies the localization map to the polynomial, eliminating
  --   those xij's for which locmap(i,j) is either 0 or 1.

    res : Poly := Null_Poly;

    procedure Localize_Term ( t : in Term; continue : out boolean ) is

      lt : Term := Localize(locmap,t);

    begin
      Add(res,lt);
      Clear(lt.dg);
      continue := true;
    end Localize_Term;
    procedure Localize_Terms is new Visiting_Iterator(Localize_Term);

  begin
    Localize_Terms(p);
    return res;
  end Localize;

-- TARGET ROUTINES :

  function Standard_Coordinate_Frame
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               plane : Standard_Complex_Matrices.Matrix )
             return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(x'range(1),x'range(2));
    tol : constant double_float := 10.0**(-10);
    first : boolean;

  begin
    for j in res'range(2) loop
      first := true;
      for i in res'range(1) loop
        if x(i,j) = Null_Poly then
          res(i,j) := 0;
        elsif (first and (AbsVal(plane(i,j)) > tol)) then
          res(i,j) := 1; first := false;
        else
          res(i,j) := 2;
        end if;
      end loop;
    end loop;
    return res;
  end Standard_Coordinate_Frame;

  function Maximal_Coordinate_Frame
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               plane : Standard_Complex_Matrices.Matrix )
             return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(x'range(1),x'range(2));
    max,tmp : double_float;
    ind : integer32;

  begin
    for j in res'range(2) loop
      for i in res'range(1) loop
        if x(i,j) = Null_Poly
         then res(i,j) := 0;
         else res(i,j) := 2;
        end if;
      end loop;
      max := 0.0; ind := 0;
      for i in res'range(1) loop
        tmp := AbsVal(plane(i,j));
        if tmp > max
         then max := tmp; ind := i;
        end if;
      end loop;
      res(ind,j) := 1;
    end loop;
    return res;
  end Maximal_Coordinate_Frame;

  function Localize ( locmap : Standard_Natural_Matrices.Matrix;
                      p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Localize(locmap,p(i));
    end loop;
    return res;
  end Localize;

-- CONSTRUCT THE POLYNOMIAL EQUATIONS :

  function Polynomial_Equations
              ( L : Standard_Complex_Matrices.Matrix;
                x : Standard_Complex_Poly_Matrices.Matrix ) return Poly_Sys is

    n : constant natural32 := natural32(x'length(1));
    p : constant natural32 := natural32(x'length(2));
    kd : constant natural32 := natural32(p + L'length(2));
    bm : Bracket_Monomial := Maximal_Minors(n,kd);
    bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(kd,kd-p,bm);
    res : constant Poly_Sys(bs'first+1..bs'last) := Expanded_Minors(L,x,bs);

  begin
    Clear(bm); Clear(bs);
    return res;
  end Polynomial_Equations;

  procedure Concat ( L : in out Link_to_Poly_Sys; p : Poly_Sys ) is
  begin
    if L = null then
      declare
        newsys : Poly_Sys(p'range);
        cnt : integer32 := p'first-1;
      begin
        for i in p'range loop
          if p(i) /= Null_Poly then
            cnt := cnt+1;
            newsys(cnt) := p(i);
          end if;
        end loop;
        L := new Poly_Sys'(newsys(p'first..cnt));
      end;
    else
      declare
        newsys : Poly_Sys(L'first..L'last+p'length);
        ind : integer32 := L'last;
      begin
        newsys(l'range) := L(L'range);
        for i in p'range loop
          if p(i) /= Null_Poly then
            ind := ind+1;
            newsys(ind) := p(i);
          end if;
        end loop;
        L := new Poly_Sys'(newsys(L'first..ind));
      end;
    end if;
  end Concat;

  function Polynomial_Equations
              ( L : Standard_Complex_VecMats.VecMat;
                x : Standard_Complex_Poly_Matrices.Matrix ) return Poly_Sys is

    n : constant natural32 := natural32(x'length(1));
    p : constant natural32 := natural32(x'length(2));
    kd : constant natural32 := p + natural32(L(L'first)'length(2));
    bm : Bracket_Monomial := Maximal_Minors(n,kd);
    bs : Bracket_System(0..integer32(Number_of_Brackets(bm)))
       := Minor_Equations(kd,kd-p,bm);
    nb : constant integer32 := L'length;
    res : Poly_Sys(bs'first+1..nb*bs'last);
    cnt : integer32 := bs'first;

  begin
    for i in 1..nb loop
      declare 
        sys : constant Poly_Sys := Expanded_Minors(L(i).all,x,bs);
      begin
        for j in sys'range loop
          cnt := cnt + 1;
          res(cnt) := sys(j);
        end loop;
      end;
    end loop;
    Clear(bm); Clear(bs);
    return res;
  end Polynomial_Equations;

-- EVALUATORS AND DIFFERENTIATORS :

  function Eval ( L,x : Standard_Complex_Matrices.Matrix )
                return Complex_Number is

    wrk : Standard_Complex_Matrices.Matrix(x'range(1),x'range(1));

  begin
    for j in L'range(2) loop
      for i in L'range(1) loop
        wrk(i,j) := L(i,j);
      end loop;
    end loop;
    for j in x'range(2) loop
      for i in x'range(1) loop
        wrk(i,L'last(2)+j) := x(i,j);
      end loop;
    end loop;
    return Determinant(wrk);
  end Eval;

  function Eval ( L,x : Standard_Complex_Matrices.Matrix ) return Vector is

    n : constant integer32 := L'length(2) + x'length(2);
    dim : constant integer32 := Number_of_Maximal_Minors(L'length(1),n);
    res : Standard_Complex_Vectors.Vector(1..dim);
      -- := (1..dim => Create(0.0));
    cnt : integer32 := 0;
    rws : Standard_Integer_Vectors.Vector(1..n);
    wrk : Standard_Complex_Matrices.Matrix(1..n,1..n);

    procedure Choose_next_Row ( k,start : in integer32 ) is

    -- DESCRIPTION :
    --   Chooses next k-th row within the range start..l'last(1), if k <= n.
    --   If k > n, then the minor defined by the current row selection
   --    is computed and added to the result.

    begin
      if k > n then
        for j in L'range(2) loop              -- select rows
          for i in 1..n loop
            wrk(i,j) := L(rws(i),j);
          end loop;
        end loop;
        for j in x'range(2) loop
          for i in 1..n loop
            wrk(i,L'last(2)+j) := x(rws(i),j);
          end loop;
        end loop;
        cnt := cnt + 1;
        res(cnt) := Determinant(wrk);
      else
        for i in start..L'last(1) loop
          rws(k) := i;
          Choose_next_Row(k+1,i+1);
        end loop;
      end if;
    end Choose_next_Row;

  begin
    Choose_next_Row(1,1);
    return res;
  end Eval;

  function Minor ( L,x : Standard_Complex_Matrices.Matrix;
                  row,col : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the minor [L|x] with row and column removed and with
  --   the appropriate sign.

    wrk : Standard_Complex_Matrices.Matrix
            (x'first(1)..x'last(1)-1,x'first(1)..x'last(1)-1);
    ind : integer32;

  begin
    for j in L'range(2) loop
      for i in L'range(1) loop
        if i < row then
          wrk(i,j) := L(i,j);
        elsif i > row then
          wrk(i-1,j) := L(i,j);
        end if;
      end loop;
    end loop;
    for j in x'range(2) loop
      if j < col then
        ind := l'last(2) + j;
      elsif j > col then
        ind := l'last(2) + j - 1;
      end if;
      if j /= col then
        for i in x'range(1) loop
          if i < row then
            wrk(i,ind) := x(i,j);
          elsif i > row then
            wrk(i-1,ind) := x(i,j);
          end if;
        end loop;
      end if;
    end loop;
    if (row + L'last(2) + col) mod 2 = 0
	 then return Determinant(wrk);
	 else return -Determinant(wrk);
    end if;
  end Minor;

  function Diff ( L,x : Standard_Complex_Matrices.Matrix;
                  i : integer32 ) return Complex_Number is

    p : constant integer32 := x'length(2);
    row,col : integer32;

  begin
    row := i/p + 1;
    col := i mod p;
    if col = 0 then
      col := p;
      row := row - 1;
    end if;
    return Minor(L,x,row,col);
  end Diff;

  function Diff ( L,x : Standard_Complex_Matrices.Matrix;
                  locmap : Standard_Natural_Matrices.Matrix;
                  i : integer32 ) return Complex_Number is

    row,col : integer32;
    cnt : integer32 := 0;

  begin
    for k1 in locmap'range(1) loop
      for k2 in locmap'range(2) loop
        if locmap(k1,k2) = 2 then
          cnt := cnt+1;
          if cnt = i then
            row := k1;
            col := k2;
          end if;
        end if;
        exit when (cnt = i);
      end loop;
      exit when (cnt = i);
    end loop;
    return Minor(l,x,row,col);
  end Diff;

  function Eval ( L : Standard_Complex_VecMats.VecMat;
                  x : Standard_Complex_Matrices.Matrix )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(l'range);

  begin
    for i in L'range loop
      res(i) := Eval(L(i).all,x);
    end loop;
    return res;
  end Eval;

  function Diff ( L : Standard_Complex_VecMats.VecMat;
                  x : Standard_Complex_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(l'range,1..x'last(1)*x'last(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Diff(l(i).all,x,j);
      end loop;
    end loop;
    return res;
  end Diff;

--  function Old_Diff ( l : Standard_Complex_VecMats.VecMat;
--                      x : Standard_Complex_Matrices.Matrix; nvars : natural;
--                      locmap : Standard_Natural_Matrices.Matrix )
--                    return Standard_Complex_Matrices.Matrix is
--
--    res : Standard_Complex_Matrices.Matrix(l'range,1..nvars);
--
--  begin
--    for i in res'range(1) loop
--      for j in res'range(2) loop
--        res(i,j) := Diff(l(i).all,x,locmap,j);
--      end loop;
--    end loop;
--    return res;
--  end Old_Diff;

  function Diff ( L : Standard_Complex_VecMats.VecMat;
                  x : Standard_Complex_Matrices.Matrix; nvars : integer32;
                  locmap : Standard_Natural_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix is

  -- NOTE :
  --   This Diff is organized according to the localization map,
  --   to avoid multiple searches.

    res : Standard_Complex_Matrices.Matrix(L'range,1..nvars);
    ind : integer32 := 0;

  begin
    for i in locmap'range(1) loop
      for j in locmap'range(2) loop
        if locmap(i,j) = 2 then
          ind := ind+1;
          for k in L'range loop
            res(k,ind) := Minor(l(k).all,x,i,j);
          end loop;
        end if;
      end loop;
    end loop;
    return res;
  end Diff;

-- SOLUTIONS AND SYSTEMS FOR QUANTUM PIERI :

  function Solution_Plane
              ( top,bottom : Bracket;
                locmap : Standard_Natural_Matrices.Matrix;
                mat : Standard_Complex_Matrices.Matrix )
              return Solution is

    solloc : constant Standard_Complex_Vectors.Vector
           := Column_Vector_Rep(locmap,Localize(locmap,mat));
    sol : Solution(solloc'length);

  begin
    sol.m := 1;
    sol.t := Create(0.0);
    sol.res := 0.0;
    sol.err := 0.0;
    sol.rco := 0.0;
    sol.v := solloc;
    return sol;
  end Solution_Plane;

  function Solution_Planes 
             ( top,bottom : Bracket;
               locmap : Standard_Natural_Matrices.Matrix;
               vm : Standard_Complex_VecMats.VecMat ) return Solution_List is

    res,res_last : Solution_List;

  begin
    for i in vm'range loop
      Append(res,res_last,Solution_Plane(top,bottom,locmap,vm(i).all));
    end loop;
    return res;
  end Solution_Planes;

  function Create_Polynomial_System
             ( top,bottom : Bracket;
               locmap : Standard_Natural_Matrices.Matrix;
               xpm : Standard_Complex_Poly_Matrices.Matrix;
               svals : Standard_Complex_Vectors.Vector;
               planes : Standard_Complex_VecMats.VecMat ) return Poly_Sys is

    res,wrksys : Poly_Sys(svals'range);

  begin
    for i in svals'range loop
      declare
        eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
            := Elim(xpm,svals(i),Create(1.0));
        wrk : constant Poly_Sys := Polynomial_Equations(planes(i).all,eva);
      begin
        wrksys(i) := wrk(wrk'first);
        Standard_Complex_Poly_Matrices.Clear(eva);
      end;
    end loop;
    res := Column_Localize(top,bottom,locmap,wrksys);
    Clear(wrksys);
    return res;
  end Create_Polynomial_System;

end Determinantal_Systems;
