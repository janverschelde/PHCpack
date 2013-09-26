with unchecked_deallocation;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with QuadDobl_Complex_VecMats;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Row_Reduction;     use QuadDobl_Complex_Row_Reduction;
with Generic_Lists;
with Lexicographic_Root_Enumeration;

package body QuadDobl_Linear_Product_System is

-- DATA STRUCTURES :

  package List_of_Vectors is
    new Generic_Lists(QuadDobl_Complex_Vectors.Link_to_Vector);
  type Equation_List is new List_of_Vectors.List;
  
  type Equation is record
    first,last : Equation_List;
  end record;

  type Equations is array ( natural32 range <> ) of Equation;
  type Link_To_Equations is access Equations;

-- INTERNAL DATA :
 
  rps : Link_To_Equations;

  Bad_Condition : constant double_float := 10.0**(-12);

  getpiv : Standard_Integer_VecVecs.Link_to_VecVec;
  getsys : QuadDobl_Complex_VecMats.Link_to_VecMat;
  getpos : Standard_Natural_Vectors.Link_to_Vector;
  getdeg : Standard_Natural_Vectors.Link_to_Vector;

-- CONSTRUCTORS :

  procedure Init ( n : in natural32 ) is
  begin
    rps := new Equations(1..n);
  end Init;

  procedure Add_Hyperplane 
              ( i : in natural32; h : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Vectors;

    eqi : Equation renames rps(i);
    lh : constant Link_To_Vector := new Vector'(h);

  begin
    if Is_Null(eqi.first) then
      Construct(lh,eqi.first);
      eqi.last := eqi.first;
    else
      declare 
        temp : Equation_List;
      begin
        Construct(lh,temp);
        Swap_Tail(eqi.last,temp);
        eqi.last := Tail_Of(eqi.last);
      end;
    end if;
  end Add_Hyperplane;

  function Dimension return natural32 is
  begin
    if rps = null
     then return 0;
     else return rps'last;
    end if;
  end Dimension;

  function Number_of_Hyperplanes ( i : natural32 ) return natural32 is
  begin
    if rps = null
     then return 0;
     else return Length_Of(rps(i).first);
    end if;
  end Number_Of_Hyperplanes;
 
  function Get_Hyperplane
             ( i,j : in natural32 )
             return QuadDobl_Complex_Vectors.Link_to_Vector is

    use QuadDobl_Complex_Vectors;

    nulvec : constant Link_to_Vector := null;

  begin
    if rps = null then
      return nulvec;
    else
      declare
        eqi : Equation_List := rps(i).first;
        count : natural32 := 1;
      begin
        while not Is_Null(eqi) loop
          if count = j
           then return Head_Of(eqi);
           else count := count + 1; eqi := Tail_Of(eqi);
          end if;
        end loop;
      end;
      return nulvec;
    end if;
  end Get_Hyperplane;

  function Get_Hyperplane
             ( i,j : in natural32 ) return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Vectors;

    lres : constant Link_to_Vector := Get_Hyperplane(i,j);
    nulvec : constant Vector(0..0) := (0..0 => Create(integer(0)));
 
  begin
    if lres = null
     then return nulvec;
     else return lres.all;
    end if;
  end Get_Hyperplane;

  procedure Enumerate_Hyperplanes ( i : in natural32 ) is
  begin
    if rps /= null then
      declare
        eqi : Equation_List := rps(i).first;
        cont : boolean;
      begin
        while not Is_Null(eqi) loop
          process(Head_Of(eqi).all,cont);
          exit when not cont;
          eqi := Tail_Of(eqi);
        end loop;
      end;
    end if;
  end Enumerate_Hyperplanes;

-- SOLVERS :

  procedure Linear_System 
              ( s : in Standard_Natural_Vectors.Vector;
                fail : out boolean;
                A : out QuadDobl_Complex_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector ) is
  begin
    if rps = null then
      fail := true;
    else
      if s'first /= integer32(rps'first) 
       or s'last /= integer32(rps'last) then
        fail := true;
      else 
        fail := false;
        declare
          use QuadDobl_Complex_Vectors;
          lv : Link_to_Vector;
        begin
          for i in s'range loop
            lv := Get_Hyperplane(natural32(i),s(i));
            if lv = null
             then fail := true; exit;
            end if;
            b(i) := -lv(0);
            for j in 1..lv'last loop
              A(i,j) := lv(j);
            end loop;
          end loop;
        end;
      end if;
    end if;
  end Linear_System;

  procedure Solve ( s : in Standard_Natural_Vectors.Vector;
                    tol : in double_float; rcond : out quad_double;
                    fail : out boolean;
                    v : out QuadDobl_Complex_Vectors.Vector ) is
  begin
    if rps = null then
      fail := true;
    else
      declare
        A : QuadDobl_Complex_Matrices.Matrix(s'range,s'range);
        ipvt : Standard_Integer_Vectors.Vector(s'range);
      begin
        Linear_System(s,fail,A,v);
        if not fail then
          lufco(A,A'last(1),ipvt,rcond);
          if rcond < tol
           then fail := true;
           else fail := false; lusolve(A,A'last(1),ipvt,v);
          end if;
        end if;
      end;
    end if;
  end Solve;

  procedure Solve ( i,n : in integer32;
                    sols,sols_last : in out Solution_List;
                    a : in out QuadDobl_Complex_Matrices.Matrix;
                    b : in out QuadDobl_Complex_Vectors.Vector; 
                    nb : in out natural32 ) is
  begin
    if i > n then
      declare
        aa : QuadDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));
        bb : QuadDobl_Complex_Vectors.Vector(b'range);
        rcond : quad_double;
        ipvt : Standard_Integer_Vectors.Vector(b'range);
      begin
        for k in a'range(1) loop
          for l in a'range(2) loop
            aa(k,l) := a(k,l);
          end loop;
          bb(k) := b(k);
        end loop;
        lufco(aa,n,ipvt,rcond);
        nb := nb + 1;
        if abs(rcond) > Bad_Condition then
          lusolve(aa,n,ipvt,bb);
          declare
            s : Solution(n);
          begin
            s.t := Create(integer(0)); s.m := 1; s.v := bb;
            s.err := create(integer(0));
            s.rco := rcond;
            s.res := create(integer(0));
            Append(sols,sols_last,s);
          end;
        end if;
      exception
        when constraint_error => return;
      end;
    else
      declare
        eqi : Equation_List := rps(natural32(i)).first;
        h : QuadDobl_Complex_Vectors.Vector(0..n);
        count : natural32 := 0;
      begin
        while not Is_Null(eqi) loop
          count := count + 1;
          h := Head_Of(eqi).all;
          b(i) := -h(0);
          for j in 1..n loop
            a(i,j) := h(j);
          end loop;
          Solve(i+1,n,sols,sols_last,a,b,nb);
          eqi := Tail_Of(eqi);
        end loop;
      end;
    end if;
  end Solve;

  procedure Solve ( sols : in out Solution_List; nl : out natural32 ) is

    n : constant integer32 := integer32(rps'last);
    m : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : QuadDobl_Complex_Vectors.Vector(1..n);
    num : natural32 := 0;
    last : Solution_List;

  begin
    for i in 1..n loop
      for j in 1..n loop
        m(i,j) := Create(integer(0));
      end loop;
      v(i) := Create(integer(0));
    end loop;
    Solve(1,n,sols,last,m,v,num);
    nl := num;
  end Solve;

  procedure Solve ( sols : in out Solution_List; nl : out natural32;
                    l : in List ) is

    n : constant integer32 := integer32(rps'last);
    m : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : QuadDobl_Complex_Vectors.Vector(1..n);
    num : natural32 := 0;
    temp : List := l;
    pos : Standard_Integer_Vectors.Vector(1..n);
    stop : boolean := false;
    last : Solution_List;

    procedure PSolve ( i,n : in integer32;
                       sols,sols_last : in out Solution_List;
                       a : in out QuadDobl_Complex_Matrices.Matrix;
                       b : in out QuadDobl_Complex_Vectors.Vector;
		       nb : in out natural32 ) is
    begin
      if i > n then
        declare
          aa : QuadDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));
          bb : QuadDobl_Complex_Vectors.Vector(b'range);
          rcond : quad_double;
          ipvt : Standard_Integer_Vectors.Vector(b'range);
        begin
          for k in a'range(1) loop
            for l in a'range(2) loop
              aa(k,l) := a(k,l);
            end loop;
            bb(k) := b(k);
          end loop;
          lufco(aa,n,ipvt,rcond);
          nb := nb + 1;
          if abs(rcond) > Bad_Condition then
            lusolve(aa,n,ipvt,bb);
            declare
              s : Solution(n);
            begin
              s.t := Create(integer(0)); s.m := 1; s.v := bb;
              s.err := create(0.0);
              s.rco := rcond;
              s.res := create(0.0);
              Append(sols,sols_last,s);
            end;
          end if;
          if Is_Null(temp)
           then stop := true;
           else pos := Head_Of(temp).all;
                temp := Tail_Of(temp);
          end if;
        exception
          when constraint_error => return;
        end;
      else
        declare
          eqi : Equation_List := rps(natural32(i)).first;
          h : QuadDobl_Complex_Vectors.Vector(0..n);
          count : natural32 := 0;
        begin
          while not Is_Null(eqi) loop
            count := count + 1;
            if integer32(count) = pos(i) then
              h := Head_Of(eqi).all;
              b(i) := -h(0);
              for j in 1..n loop
                a(i,j) := h(j);
              end loop;
             --put("eq"); put(i,1); put(count,1); put(" ");
              PSolve(i+1,n,sols,sols_last,a,b,nb);
            end if;
            exit when stop;
            eqi := Tail_Of(eqi);
          end loop;
        end;
      end if;
    end PSolve;

  begin
    if not Is_Null(temp) then
      pos := Head_Of(temp).all;
      temp := Tail_Of(temp);
      for i in 1..n loop
        for j in 1..n loop
          m(i,j) := Create(integer(0));
        end loop;
        v(i) := Create(integer(0));
      end loop;
      PSolve(1,n,sols,last,m,v,num);
      nl := num;
    end if;
  end Solve;

-- ENUMERATORS :

  procedure Enumerate_Solutions
              ( tol : in double_float; rc : out natural32 ) is

    n : constant integer32
      := integer32(QuadDobl_Linear_Product_System.Dimension);
    pos : Standard_Natural_Vectors.Vector(1..n);
    hyp : QuadDobl_Complex_Vectors.Vector(0..n);
    singular : boolean;
    cont : boolean := true;

    procedure Count ( row : in integer32;
                      piv : in Standard_Integer_Vectors.Vector;
                      sys : in QuadDobl_Complex_Matrices.Matrix ) is

      newpiv : Standard_Integer_Vectors.Vector(piv'range);
      newsys : QuadDobl_Complex_Matrices.Matrix(sys'range(1),sys'range(2));
   
    begin
      if row > n then
        rc := rc + 1;
        Process(pos,cont);
      else
        for k in 1..Number_of_Hyperplanes(natural32(row)) loop
          pos(row) := k;
          hyp := QuadDobl_Linear_Product_System.Get_Hyperplane
                   (natural32(row),k);
          for j in sys'range(2) loop
            for i in sys'first(1)..(row-1) loop
              newsys(i,j) := sys(i,j);
            end loop;
            newsys(row,j) := hyp(piv(j));
          end loop;
          newpiv := piv;
          Reduce_Row(newsys,row,newpiv,tol,singular);
          if not singular
           then Count(row+1,newpiv,newsys);
          end if;
          exit when not cont;
        end loop;
      end if;
    end Count;

  begin
    rc := 0;
    declare
      piv : constant Standard_Integer_Vectors.Vector(1..n) := Start_Pivots(n);
      sys : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    begin
      sys(1,1) := Create(integer(0)); -- to avoid a compiler warning ...
      Count(1,piv,sys);
    end;
  end Enumerate_Solutions;

  function Count_All_Solutions ( tol : double_float ) return natural32 is

   -- res : natural := 0;
    rc : natural32 := 0;

    procedure Just_Continue ( s : in Standard_Natural_Vectors.Vector;
                              c : out boolean ) is
    begin
     -- res := res + 1;
     -- put(s); put(" : + 1 = "); put(res,1); new_line;
      c := true;
    end Just_Continue;
    procedure Count is new Enumerate_Solutions(Just_Continue);

  begin
    Count(tol,rc);
    return rc;
  end Count_All_Solutions;

  procedure Write_Hyp ( file : in file_type;
                        h : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the patterns of the coefficients in h
  --   to check whether the pivoting information is correct.

  begin
    for i in 1..h'last loop
      if AbsVal(h(i)) < 1.0E-8
       then put(" 0");
       else put(" *");
      end if;
    end loop;
  end Write_Hyp;

  procedure Get_First
               ( tol : in double_float;
                 s : out Standard_Natural_Vectors.Vector;
                 fail : out boolean ) is

    n : constant integer32
      := integer32(QuadDobl_Linear_Product_System.Dimension);
    pos : Standard_Natural_Vectors.Vector(1..n);
    hyp : QuadDobl_Complex_Vectors.Vector(0..n);
    singular : boolean;
    stop : boolean := false;

    procedure Count ( row : in integer32 ) is

      newpiv : Standard_Integer_Vectors.Vector(1..n);
      newsys : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
   
    begin
      if row > n then
        s := pos;
        stop := true; fail := false;
      else
        for k in 1..Number_of_Hyperplanes(natural32(row)) loop
          pos(row) := k;
          hyp := QuadDobl_Linear_Product_System.Get_Hyperplane
                   (natural32(row),k);
          for j in 1..n loop
            for i in 1..(row-1) loop
              newsys(i,j) := getsys(row-1)(i,j);
            end loop;
            newsys(row,j) := hyp(getpiv(row)(j));
          end loop;
          newpiv := getpiv(row).all;
          Reduce_Row(newsys,row,newpiv,tol,singular);
          if not singular then
            if row < n then
              for i in 1..n loop
                getpiv(row+1)(i) := newpiv(i);
                for j in 1..n loop
                  getsys(row)(i,j) := newsys(i,j);
                end loop;
              end loop;
            end if;
            Count(row+1);
          end if;
          exit when stop;
        end loop;
      end if;
    end Count;

  begin
    fail := true;
    getpiv := new Standard_Integer_VecVecs.VecVec(1..n);
    getpiv(1) := new Standard_Integer_Vectors.Vector'(Start_Pivots(n));
    for i in 2..n loop
      getpiv(i) := new Standard_Integer_Vectors.Vector(1..n);
    end loop;
    getsys := new QuadDobl_Complex_VecMats.VecMat(1..n);
    for i in 1..n loop
      declare
        sys : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
      begin
        sys(1,1) := Create(integer(0));  -- to avoid a compiler warning ...
        getsys(i) := new QuadDobl_Complex_Matrices.Matrix'(sys);
      end;
    end loop;
    Count(1);
    getpos := new Standard_Natural_Vectors.Vector'(pos);
    getdeg := new Standard_Natural_Vectors.Vector(1..n);
    for i in 1..n loop
      getdeg(i) := Number_of_Hyperplanes(natural32(i));
    end loop;
  end Get_First;

  procedure Get_First
              ( file : in file_type; tol : in double_float;
                s : out Standard_Natural_Vectors.Vector;
                fail : out boolean ) is

    n : constant integer32
      := integer32(QuadDobl_Linear_Product_System.Dimension);
    pos : Standard_Natural_Vectors.Vector(1..n);
    hyp : QuadDobl_Complex_Vectors.Vector(0..n);
    singular : boolean;
    stop : boolean := false;

    procedure Count ( row : in integer32 ) is

      newpiv : Standard_Integer_Vectors.Vector(1..n);
      newsys : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
   
    begin
      put(file,"in Count with row = ");
      put(file,row,1); put_line(file," :");
      if row > n then
        s := pos;
        stop := true; fail := false;
      else
        for k in 1..Number_of_Hyperplanes(natural32(row)) loop
          put(file,"in loop with k = "); put(file,k,1);
          put(file," and row = "); put(file,row,1); new_line(file);
          pos(row) := k;
          hyp := QuadDobl_Linear_Product_System.Get_Hyperplane
                   (natural32(row),k);
          put(file,"Hyperplane has pattern ");
          Write_Hyp(file,hyp); new_line(file);
          for j in 1..n loop
            for i in 1..(row-1) loop
              newsys(i,j) := getsys(row-1)(i,j);
            end loop;
            newsys(row,j) := hyp(getpiv(row)(j));
          end loop;
          newpiv := getpiv(row).all;
          put(file,"Pivots from row "); put(file,row,1);
          put(file," : "); put(file,newpiv);
          Reduce_Row(newsys,row,newpiv,tol,singular);
          put(file," after reduce pivots : ");
          put(file,newpiv); new_line(file);
          if not singular then
            if row < n then
              for i in 1..n loop
                getpiv(row+1)(i) := newpiv(i);
                for j in 1..n loop
                  getsys(row)(i,j) := newsys(i,j);
                end loop;
              end loop;
            end if;
            put(file," not singular at row = ");
            put(file,row,1); new_line(file);
            Count(row+1);
          else
            put(file,"singular after reduced move to next k...");
          end if;
          exit when stop;
        end loop;
      end if;
    end Count;

  begin
    fail := true;
    getpiv := new Standard_Integer_VecVecs.VecVec(1..n);
    getpiv(1) := new Standard_Integer_Vectors.Vector'(Start_Pivots(n));
    for i in 2..n loop
      getpiv(i) := new Standard_Integer_Vectors.Vector(1..n);
    end loop;
    getsys := new QuadDobl_Complex_VecMats.VecMat(1..n);
    for i in 1..n loop
      declare
        sys : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
      begin
        sys(1,1) := Create(integer(0));  -- to avoid a compiler warning ...
        getsys(i) := new QuadDobl_Complex_Matrices.Matrix'(sys);
      end;
    end loop;
    Count(1);
    put_line(file,"Contents of getpiv after Get_First : ");
    for i in getpiv'range loop
      put(file,i,1); put(file," : ");
      put(file,getpiv(i).all); new_line(file);
    end loop;
  end Get_First;

  procedure Get_Next
              ( tol : in double_float;
                d : in Standard_Natural_Vectors.Vector;
                s : in out Standard_Natural_Vectors.Vector;
                fail : out boolean ) is

    row : integer32;
    n : constant integer32 := s'last;
    hyp : QuadDobl_Complex_Vectors.Vector(0..n);
    piv : Standard_Integer_Vectors.Vector(1..n);
    sys : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    singular : boolean;

  begin
    Lexicographic_Root_Enumeration.Add_One(d,s,row,fail);
    while not fail loop
      hyp := QuadDobl_Linear_Product_System.Get_Hyperplane
               (natural32(row),s(row));
      for i in 1..row-1 loop
        for j in 1..n loop
          sys(i,j) := getsys(row-1)(i,j);
        end loop;
      end loop;
      piv := getpiv(row).all;
      for j in 1..n loop
        sys(row,j) := hyp(piv(j));
      end loop;
      Reduce_Row(sys,row,piv,tol,singular);
      exit when (row = n) and not singular;
      if singular then
        Lexicographic_Root_Enumeration.Add_One(d(1..row),s(1..row),row,fail);
        for i in row+1..s'last loop
          s(i) := 1;
        end loop;
      else
        getsys(row).all := sys;
        getpiv(row+1).all := piv;
       -- put("s ="); put(s); put(" + 1 =");
       -- Lexicographic_Root_Enumeration.Add_One(d,s,row,fail);
        row := row + 1;
      end if;
    end loop;
  end Get_Next;

  procedure Get_Next
              ( file : in file_type; tol : in double_float;
                d : in Standard_Natural_Vectors.Vector;
                s : in out Standard_Natural_Vectors.Vector;
                fail : out boolean ) is

    row : integer32;
    n : constant integer32 := s'last;
    hyp : QuadDobl_Complex_Vectors.Vector(0..n);
    piv : Standard_Integer_Vectors.Vector(1..n);
    sys : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    singular : boolean;

  begin
    put(file,"Calling Get_Next with s ="); put(file,s); put(file," + 1 =");
    Lexicographic_Root_Enumeration.Add_One(d,s,row,fail);
    put(file,s); put(file," row = "); put(file,row,1); new_line(file);
    while not fail loop
      hyp := QuadDobl_Linear_Product_System.Get_Hyperplane
               (natural32(row),s(row));
      put(file,"Retrieving coefficients of system of row = ");
      put(file,row-1,1); new_line(file);
      for i in 1..row-1 loop
        for j in 1..n loop
          sys(i,j) := getsys(row-1)(i,j);
        end loop;
      end loop;
      piv := getpiv(row).all;
      put(file,"Pivots used from row = "); put(file,row,1);
      put(file," : "); put(file,piv); new_line(file);
      for j in 1..n loop
        sys(row,j) := hyp(piv(j));
      end loop;
      Reduce_Row(sys,row,piv,tol,singular);
      exit when (row = n) and not singular;
      if singular then
        put(file,"singular after reduce s =");
        put(file,s); put(file," + 1 =");
        Lexicographic_Root_Enumeration.Add_One(d(1..row),s(1..row),row,fail);
        for i in row+1..s'last loop
          s(i) := 1;
        end loop;
        put(file,s); put(file," row = "); put(file,row,1); new_line(file);
      else
        put(file,"not singular but row = "); put(file,row,1);
        put(file," < "); put(file,n,1); put_line(file," storing data");
        getsys(row).all := sys;
        getpiv(row+1).all := piv;
       -- put("s ="); put(s); put(" + 1 =");
       -- Lexicographic_Root_Enumeration.Add_One(d,s,row,fail);
        row := row + 1;
        put(file,s); put(file," row = "); put(file,row,1); new_line(file);
      end if;
    end loop;
  end Get_Next;

  procedure Get_Next 
              ( tol : in double_float;
                s : out Standard_Natural_Vectors.Vector;
                fail : out boolean ) is
  begin
    Get_Next(tol,getdeg.all,getpos.all,fail);
    s := getpos.all;
  end Get_Next;

  procedure Get_Clear is
  begin
    Standard_Integer_VecVecs.Deep_Clear(getpiv);
    QuadDobl_Complex_VecMats.Deep_Clear(getsys);
    Standard_Natural_Vectors.Clear(getpos);
    Standard_Natural_Vectors.Clear(getdeg);
  end Get_Clear;

-- EXPAND LINEAR-PRODUCT SYSTEM INTO POLYNOMIAL SYSTEM :

  function Polynomial ( h : QuadDobl_Complex_Vectors.Vector ) return Poly is
 
    res : Poly;
    t : Term;
    n : constant integer32 := h'last;

  begin
    for j in 0..n loop
      if h(j) /= Create(integer(0)) then
        t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
        t.cf := h(j);
        if j /= 0
         then t.dg(j) := 1;
        end if;
        Add(res,t);
        Clear(t.dg);
      end if;
    end loop;
    return res;
  end Polynomial;

  function Create ( i : in natural32 ) return Poly is

    eql : Equation_List := rps(i).first;
    hyp,res : Poly := Null_Poly;

  begin
    while not Is_Null(eql) loop
      hyp := Polynomial(Head_Of(eql).all);
      if res = Null_Poly
       then Copy(hyp,res);
       else Mul(res,hyp);
      end if;
      Clear(hyp);
      eql := Tail_Of(eql);
    end loop;
    return res;
  end Create;

  function Polynomial_System return Poly_Sys is

    res : Poly_Sys(integer32(rps'first)..integer32(rps'last));

  begin
    for i in rps'range loop
      res(integer32(i)) := Create(i);
    end loop;
    return res;
  end Polynomial_System;

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation(Equations,Link_To_Equations);

  procedure Clear ( eql : in out Equation_List ) is

    tmp : Equation_List := eql;
    lv : QuadDobl_Complex_Vectors.Link_To_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      QuadDobl_Complex_Vectors.Clear(lv);
      tmp := Tail_of(tmp);
    end loop;
    List_Of_Vectors.Clear(List_Of_Vectors.List(eql));
  end Clear;

  procedure Clear ( eq : in out Equation ) is
  begin
    Clear(eq.first);
    -- eq.last is just a pointer to the last element of eq.first;
    -- if eq.first disappears, then also eq.last does
  end Clear;

  procedure Clear ( eqs : in out Equations ) is
  begin
    for i in eqs'range loop
      Clear(eqs(i));
    end loop;
  end Clear;

  procedure Clear is
  begin
    if rps /= null then
      for i in rps'range loop
        Clear(rps(i));
      end loop;
      free(rps);
    end if;
  end Clear;

end QuadDobl_Linear_Product_System;
