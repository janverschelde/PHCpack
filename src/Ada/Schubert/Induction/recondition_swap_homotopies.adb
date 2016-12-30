with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Double_Double_Numbers;               use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors;
with Checker_Localization_Patterns;

package body Recondition_Swap_Homotopies is

  procedure Insert_One_Variable
              ( t : in out Standard_Complex_Polynomials.Term;
                k : in integer32 ) is
  
    newdeg : Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);

  begin
    for i in t.dg'first..(k-1) loop
      newdeg(i) := t.dg(i);
    end loop;
    newdeg(k) := 0;
    for i in k..t.dg'last loop
      newdeg(i+1) := t.dg(i);
    end loop;
    Standard_Complex_Polynomials.Clear(t.dg);
    t.dg := new Standard_Natural_Vectors.Vector'(newdeg);
  end Insert_One_Variable;

  procedure Insert_One_Variable
              ( t : in out DoblDobl_Complex_Polynomials.Term;
                k : in integer32 ) is
  
    newdeg : Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);

  begin
    for i in t.dg'first..(k-1) loop
      newdeg(i) := t.dg(i);
    end loop;
    newdeg(k) := 0;
    for i in k..t.dg'last loop
      newdeg(i+1) := t.dg(i);
    end loop;
    DoblDobl_Complex_Polynomials.Clear(t.dg);
    t.dg := new Standard_Natural_Vectors.Vector'(newdeg);
  end Insert_One_Variable;

  procedure Insert_One_Variable
              ( t : in out QuadDobl_Complex_Polynomials.Term;
                k : in integer32 ) is
  
    newdeg : Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);

  begin
    for i in t.dg'first..(k-1) loop
      newdeg(i) := t.dg(i);
    end loop;
    newdeg(k) := 0;
    for i in k..t.dg'last loop
      newdeg(i+1) := t.dg(i);
    end loop;
    QuadDobl_Complex_Polynomials.Clear(t.dg);
    t.dg := new Standard_Natural_Vectors.Vector'(newdeg);
  end Insert_One_Variable;

  procedure Remove_One_Variable
              ( t : in out Standard_Complex_Polynomials.Term;
                k : in integer32 ) is

    newdeg : Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);

  begin
    for i in t.dg'first..(k-1) loop
      newdeg(i) := t.dg(i);
    end loop;
    for i in (k+1)..t.dg'last loop
      newdeg(i-1) := t.dg(i);
    end loop;
    Standard_Complex_Polynomials.Clear(t.dg);
    t.dg := new Standard_Natural_Vectors.Vector'(newdeg);
  end Remove_One_Variable;

  procedure Remove_One_Variable
              ( t : in out DoblDobl_Complex_Polynomials.Term;
                k : in integer32 ) is

    newdeg : Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);

  begin
    for i in t.dg'first..(k-1) loop
      newdeg(i) := t.dg(i);
    end loop;
    for i in (k+1)..t.dg'last loop
      newdeg(i-1) := t.dg(i);
    end loop;
    DoblDobl_Complex_Polynomials.Clear(t.dg);
    t.dg := new Standard_Natural_Vectors.Vector'(newdeg);
  end Remove_One_Variable;

  procedure Remove_One_Variable
              ( t : in out QuadDobl_Complex_Polynomials.Term;
                k : in integer32 ) is

    newdeg : Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);

  begin
    for i in t.dg'first..(k-1) loop
      newdeg(i) := t.dg(i);
    end loop;
    for i in (k+1)..t.dg'last loop
      newdeg(i-1) := t.dg(i);
    end loop;
    QuadDobl_Complex_Polynomials.Clear(t.dg);
    t.dg := new Standard_Natural_Vectors.Vector'(newdeg);
  end Remove_One_Variable;

  procedure Insert_One_Variable
              ( p : in out Standard_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use Standard_Complex_Polynomials;

    procedure Insert_to_Term ( t : in out Term; c : out boolean ) is
    begin
      Insert_One_Variable(t,k);
      c := true;
    end Insert_to_Term;
    procedure Insert_to_Terms is new Changing_Iterator(Insert_to_Term);

  begin
    if p /= Null_Poly
     then Insert_to_Terms(p);
    end if;
  end Insert_One_Variable;

  procedure Insert_One_Variable
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use DoblDobl_Complex_Polynomials;

    procedure Insert_to_Term ( t : in out Term; c : out boolean ) is
    begin
      Insert_One_Variable(t,k);
      c := true;
    end Insert_to_Term;
    procedure Insert_to_Terms is new Changing_Iterator(Insert_to_Term);

  begin
    if p /= Null_Poly
     then Insert_to_Terms(p);
    end if;
  end Insert_One_Variable;

  procedure Insert_One_Variable
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use QuadDobl_Complex_Polynomials;

    procedure Insert_to_Term ( t : in out Term; c : out boolean ) is
    begin
      Insert_One_Variable(t,k);
      c := true;
    end Insert_to_Term;
    procedure Insert_to_Terms is new Changing_Iterator(Insert_to_Term);

  begin
    if p /= Null_Poly
     then Insert_to_Terms(p);
    end if;
  end Insert_One_Variable;

  procedure Remove_One_Variable
              ( p : in out Standard_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use Standard_Complex_Polynomials;

    procedure Remove_from_Term ( t : in out Term; c : out boolean ) is
    begin
      Remove_One_Variable(t,k);
      c := true;
    end Remove_from_Term;
    procedure Remove_from_Terms is new Changing_Iterator(Remove_from_Term);

  begin
    if p /= Null_Poly
     then Remove_from_Terms(p);
    end if;
  end Remove_One_Variable;

  procedure Remove_One_Variable
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use DoblDobl_Complex_Polynomials;

    procedure Remove_from_Term ( t : in out Term; c : out boolean ) is
    begin
      Remove_One_Variable(t,k);
      c := true;
    end Remove_from_Term;
    procedure Remove_from_Terms is new Changing_Iterator(Remove_from_Term);

  begin
    if p /= Null_Poly
     then Remove_from_Terms(p);
    end if;
  end Remove_One_Variable;

  procedure Remove_One_Variable
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use QuadDobl_Complex_Polynomials;

    procedure Remove_from_Term ( t : in out Term; c : out boolean ) is
    begin
      Remove_One_Variable(t,k);
      c := true;
    end Remove_from_Term;
    procedure Remove_from_Terms is new Changing_Iterator(Remove_from_Term);

  begin
    if p /= Null_Poly
     then Remove_from_Terms(p);
    end if;
  end Remove_One_Variable;

  procedure Insert_One_Variable
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                k : in integer32 ) is
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        Insert_One_Variable(x(i,j),k);
      end loop;
    end loop;
  end Insert_One_Variable;

  procedure Insert_One_Variable
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 ) is
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        Insert_One_Variable(x(i,j),k);
      end loop;
    end loop;
  end Insert_One_Variable;

  procedure Insert_One_Variable
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 ) is
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        Insert_One_Variable(x(i,j),k);
      end loop;
    end loop;
  end Insert_One_Variable;

  procedure Remove_One_Variable
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                k : in integer32 ) is
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        Remove_One_Variable(x(i,j),k);
      end loop;
    end loop;
  end Remove_One_Variable;

  procedure Remove_One_Variable
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 ) is
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        Remove_One_Variable(x(i,j),k);
      end loop;
    end loop;
  end Remove_One_Variable;

  procedure Remove_One_Variable
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                k : in integer32 ) is
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        Remove_One_Variable(x(i,j),k);
      end loop;
    end loop;
  end Remove_One_Variable;

  procedure Insert_Variable_Pivot
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                i,j,k : in integer32 ) is

    use Standard_Complex_Polynomials;

    procedure Insert_to_Term ( t : in out Term; c : out boolean ) is
    begin
      t.dg(k) := 1;
      c := true;
    end Insert_to_Term;
    procedure Insert_to_Terms is new Changing_Iterator(Insert_to_Term);

  begin
    if x(i,j) /= Null_Poly
     then Insert_to_Terms(x(i,j));
    end if;
  end Insert_Variable_Pivot;

  procedure Insert_Variable_Pivot
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                i,j,k : in integer32 ) is

    use DoblDobl_Complex_Polynomials;

    procedure Insert_to_Term ( t : in out Term; c : out boolean ) is
    begin
      t.dg(k) := 1;
      c := true;
    end Insert_to_Term;
    procedure Insert_to_Terms is new Changing_Iterator(Insert_to_Term);

  begin
    if x(i,j) /= Null_Poly
     then Insert_to_Terms(x(i,j));
    end if;
  end Insert_Variable_Pivot;

  procedure Insert_Variable_Pivot
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                i,j,k : in integer32 ) is

    use QuadDobl_Complex_Polynomials;

    procedure Insert_to_Term ( t : in out Term; c : out boolean ) is
    begin
      t.dg(k) := 1;
      c := true;
    end Insert_to_Term;
    procedure Insert_to_Terms is new Changing_Iterator(Insert_to_Term);

  begin
    if x(i,j) /= Null_Poly
     then Insert_to_Terms(x(i,j));
    end if;
  end Insert_Variable_Pivot;

  procedure Recondition
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                dim,s : in integer32 ) is

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);

  begin
    Insert_One_Variable(x,dim+1);
    Insert_Variable_Pivot(x,rowpiv,s+1,dim+1);
  end Recondition;

  procedure Recondition
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                dim,s : in integer32 ) is

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);

  begin
    Insert_One_Variable(x,dim+1);
    Insert_Variable_Pivot(x,rowpiv,s+1,dim+1);
  end Recondition;

  procedure Recondition
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                dim,s : in integer32 ) is

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);

  begin
    Insert_One_Variable(x,dim+1);
    Insert_Variable_Pivot(x,rowpiv,s+1,dim+1);
  end Recondition;

  function Random_Linear_Equation
              ( x : Standard_Complex_Poly_Matrices.Matrix;
                s : integer32 )
              return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    rnd : Complex_Number;
    acc : Poly;

  begin
    for i in x'range(1) loop
      -- avoid nonlinear terms!
      --if x(i,s) /= Null_Poly then
      --  rnd := Standard_Random_Numbers.Random1;     
      --  acc := rnd*x(i,s);
      --  Add(res,acc); Clear(acc);
      --end if;
      if x(i,s+1) /= Null_Poly then
        rnd := Standard_Random_Numbers.Random1;     
        acc := rnd*x(i,s+1);
        Add(res,acc); Clear(acc);
      end if;
    end loop;
    return res;
  end Random_Linear_Equation;

  function Random_Linear_Equation
              ( x : DoblDobl_Complex_Poly_Matrices.Matrix;
                s : integer32 )
              return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    rnd : Complex_Number;
    acc : Poly;

  begin
    for i in x'range(1) loop
      -- avoid nonlinear terms!
      --if x(i,s) /= Null_Poly then
      --  rnd := Standard_Random_Numbers.Random1;     
      --  acc := rnd*x(i,s);
      --  Add(res,acc); Clear(acc);
      --end if;
      if x(i,s+1) /= Null_Poly then
        rnd := DoblDobl_Random_Numbers.Random1;     
        acc := rnd*x(i,s+1);
        Add(res,acc); Clear(acc);
      end if;
    end loop;
    return res;
  end Random_Linear_Equation;

  function Random_Linear_Equation
              ( x : QuadDobl_Complex_Poly_Matrices.Matrix;
                s : integer32 )
              return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    rnd : Complex_Number;
    acc : Poly;

  begin
    for i in x'range(1) loop
      -- avoid nonlinear terms!
      --if x(i,s) /= Null_Poly then
      --  rnd := Standard_Random_Numbers.Random1;     
      --  acc := rnd*x(i,s);
      --  Add(res,acc); Clear(acc);
      --end if;
      if x(i,s+1) /= Null_Poly then
        rnd := QuadDobl_Random_Numbers.Random1;     
        acc := rnd*x(i,s+1);
        Add(res,acc); Clear(acc);
      end if;
    end loop;
    return res;
  end Random_Linear_Equation;

  procedure Set_Exponent_to_Zero
              ( p : in out Standard_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use Standard_Complex_Polynomials;

    procedure Set_in_Term ( t : in out Term; c : out boolean ) is
    begin
      t.dg(k) := 0;
      c := true;
    end Set_in_Term;

    procedure Set_in_Terms is new Changing_Iterator(Set_in_Term);

  begin
    Set_in_Terms(p);
  end Set_Exponent_to_Zero;

  procedure Set_Exponent_to_Zero
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use DoblDobl_Complex_Polynomials;

    procedure Set_in_Term ( t : in out Term; c : out boolean ) is
    begin
      t.dg(k) := 0;
      c := true;
    end Set_in_Term;

    procedure Set_in_Terms is new Changing_Iterator(Set_in_Term);

  begin
    Set_in_Terms(p);
  end Set_Exponent_to_Zero;

  procedure Set_Exponent_to_Zero
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                k : in integer32 ) is

    use QuadDobl_Complex_Polynomials;

    procedure Set_in_Term ( t : in out Term; c : out boolean ) is
    begin
      t.dg(k) := 0;
      c := true;
    end Set_in_Term;

    procedure Set_in_Terms is new Changing_Iterator(Set_in_Term);

  begin
    Set_in_Terms(p);
  end Set_Exponent_to_Zero;

  procedure Add_Random_Constant
              ( p : in out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.cf := Standard_Random_Numbers.Random1;
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    Add(p,t);
    Clear(t);
  end Add_Random_Constant;

  procedure Add_Random_Constant
              ( p : in out DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.cf := DoblDobl_Random_Numbers.Random1;
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    Add(p,t);
    Clear(t);
  end Add_Random_Constant;

  procedure Add_Random_Constant
              ( p : in out QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    t : Term;

  begin
    t.cf := QuadDobl_Random_Numbers.Random1;
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    Add(p,t);
    Clear(t);
  end Add_Random_Constant;

  function Recondition_Target_Equation
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               s,t : integer32 )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;

  begin
    res := Random_Linear_Equation(x,s);
    Set_Exponent_to_Zero(res,t);
    Add_Random_Constant(res);
    return res;
  end Recondition_Target_Equation;

  function Recondition_Target_Equation
             ( x : DoblDobl_Complex_Poly_Matrices.Matrix;
               s,t : integer32 )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly;

  begin
    res := Random_Linear_Equation(x,s);
    Set_Exponent_to_Zero(res,t);
    Add_Random_Constant(res);
    return res;
  end Recondition_Target_Equation;

  function Recondition_Target_Equation
             ( x : QuadDobl_Complex_Poly_Matrices.Matrix;
               s,t : integer32 )
             return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly;

  begin
    res := Random_Linear_Equation(x,s);
    Set_Exponent_to_Zero(res,t);
    Add_Random_Constant(res);
    return res;
  end Recondition_Target_Equation;

  function Recondition_Start_Equation
             ( n,k : integer32 )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;
    trm : Standard_Complex_Polynomials.term;

  begin
    trm.cf := Standard_Complex_Numbers.Create(1.0);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := Standard_Complex_Polynomials.Create(trm);
    trm.dg(k) := 0;
    Standard_Complex_Polynomials.Sub(res,trm);
    Standard_Complex_Polynomials.Clear(trm);
    return res;
  end Recondition_Start_Equation;

  function Recondition_Start_Equation
             ( n,k : integer32 )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly;
    trm : DoblDobl_Complex_Polynomials.term;

  begin
    trm.cf := DoblDobl_Complex_Numbers.Create(integer32(1));
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := DoblDobl_Complex_Polynomials.Create(trm);
    trm.dg(k) := 0;
    DoblDobl_Complex_Polynomials.Sub(res,trm);
    DoblDobl_Complex_Polynomials.Clear(trm);
    return res;
  end Recondition_Start_Equation;

  function Recondition_Start_Equation
             ( n,k : integer32 )
             return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly;
    trm : QuadDobl_Complex_Polynomials.term;

  begin
    trm.cf := QuadDobl_Complex_Numbers.Create(integer32(1));
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := QuadDobl_Complex_Polynomials.Create(trm);
    trm.dg(k) := 0;
    QuadDobl_Complex_Polynomials.Sub(res,trm);
    QuadDobl_Complex_Polynomials.Clear(trm);
    return res;
  end Recondition_Start_Equation;

  function Recondition_Equation
             ( x : Standard_Complex_Poly_Matrices.Matrix;
               s,t,k : integer32 )
             return Standard_Complex_Polynomials.Poly is

    target : Standard_Complex_Polynomials.Poly
           := Recondition_Target_Equation(x,s,t);
    start : Standard_Complex_Polynomials.Poly
          := Recondition_Start_Equation(t,k);
    trm : Standard_Complex_Polynomials.Term;

  begin
    trm.cf := Standard_Complex_Numbers.Create(1.0);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..t => 0);
    trm.dg(t) := 1; -- trm is term with the continuation parameter t
    Standard_Complex_Polynomials.Mul(target,trm);   -- t*target
    Standard_Complex_Polynomials.Add(target,start); -- t*target + start
    Standard_Complex_Polynomials.Mul(start,trm);    
    Standard_Complex_Polynomials.Sub(target,start); -- t*target + (1-t)*start
    Standard_Complex_Polynomials.Clear(start);
    Standard_Complex_Polynomials.Clear(trm);
    return target;
  end Recondition_Equation;

  function Recondition_Equation
             ( x : DoblDobl_Complex_Poly_Matrices.Matrix;
               s,t,k : integer32 )
             return DoblDobl_Complex_Polynomials.Poly is

    target : DoblDobl_Complex_Polynomials.Poly
           := Recondition_Target_Equation(x,s,t);
    start : DoblDobl_Complex_Polynomials.Poly
          := Recondition_Start_Equation(t,k);
    trm : DoblDobl_Complex_Polynomials.Term;

  begin
    trm.cf := DoblDobl_Complex_Numbers.Create(integer32(1));
    trm.dg := new Standard_Natural_Vectors.Vector'(1..t => 0);
    trm.dg(t) := 1; -- trm is term with the continuation parameter t
    DoblDobl_Complex_Polynomials.Mul(target,trm);   -- t*target
    DoblDobl_Complex_Polynomials.Add(target,start); -- t*target + start
    DoblDobl_Complex_Polynomials.Mul(start,trm);    
    DoblDobl_Complex_Polynomials.Sub(target,start); -- t*target + (1-t)*start
    DoblDobl_Complex_Polynomials.Clear(start);
    DoblDobl_Complex_Polynomials.Clear(trm);
    return target;
  end Recondition_Equation;

  function Recondition_Equation
             ( x : QuadDobl_Complex_Poly_Matrices.Matrix;
               s,t,k : integer32 )
             return QuadDobl_Complex_Polynomials.Poly is

    target : QuadDobl_Complex_Polynomials.Poly
           := Recondition_Target_Equation(x,s,t);
    start : QuadDobl_Complex_Polynomials.Poly
          := Recondition_Start_Equation(t,k);
    trm : QuadDobl_Complex_Polynomials.Term;

  begin
    trm.cf := QuadDobl_Complex_Numbers.Create(integer32(1));
    trm.dg := new Standard_Natural_Vectors.Vector'(1..t => 0);
    trm.dg(t) := 1; -- trm is term with the continuation parameter t
    QuadDobl_Complex_Polynomials.Mul(target,trm);   -- t*target
    QuadDobl_Complex_Polynomials.Add(target,start); -- t*target + start
    QuadDobl_Complex_Polynomials.Mul(start,trm);    
    QuadDobl_Complex_Polynomials.Sub(target,start); -- t*target + (1-t)*start
    QuadDobl_Complex_Polynomials.Clear(start);
    QuadDobl_Complex_Polynomials.Clear(trm);
    return target;
  end Recondition_Equation;

  function Recondition_Solution_Vector
             ( x : Standard_Complex_Vectors.Vector; k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Vectors.Vector(x'first..x'last+1);
    fac : constant Complex_Number := 1.0/x(k);
    r,c : integer32;

  begin
    for i in x'range loop
      Checker_Localization_Patterns.Position(locmap,i,r,c);
      if c = s+1 then
        res(i) := fac*x(i);
      elsif c = s then
        if Degree(xp(r,c),k) = 0   -- x(ctr+1,s+1) does not occur
         then res(i) := fac*x(i);
         else res(i) := x(i);
        end if;
      else
        res(i) := x(i);
      end if;
    end loop;
    res(res'last) := fac;
    return res;
  end Recondition_Solution_Vector;

  function Recondition_Solution_Vector
             ( x : DoblDobl_Complex_Vectors.Vector; k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);
    one : constant double_double := create(1.0);
    fac : constant Complex_Number := one/x(k);
    r,c : integer32;

  begin
    for i in x'range loop
      Checker_Localization_Patterns.Position(locmap,i,r,c);
      if c = s+1 then
        res(i) := fac*x(i);
      elsif c = s then
        if Degree(xp(r,c),k) = 0   -- x(ctr+1,s+1) does not occur
         then res(i) := fac*x(i);
         else res(i) := x(i);
        end if;
      else
        res(i) := x(i);
      end if;
    end loop;
    res(res'last) := fac;
    return res;
  end Recondition_Solution_Vector;

  function Recondition_Solution_Vector
             ( x : QuadDobl_Complex_Vectors.Vector; k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);
    one : constant quad_double := create(1.0);
    fac : constant Complex_Number := one/x(k);
    r,c : integer32;

  begin
    for i in x'range loop
      Checker_Localization_Patterns.Position(locmap,i,r,c);
      if c = s+1 then
        res(i) := fac*x(i);
      elsif c = s then
        if Degree(xp(r,c),k) = 0   -- x(ctr+1,s+1) does not occur
         then res(i) := fac*x(i);
         else res(i) := x(i);
        end if;
      else
        res(i) := x(i);
      end if;
    end loop;
    res(res'last) := fac;
    return res;
  end Recondition_Solution_Vector;

  function Recondition_Solution
             ( sol : Standard_Complex_Solutions.Solution;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(sol.n+1);

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.v := Recondition_Solution_Vector(sol.v,k,s,locmap,xp);
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    return res;
  end Recondition_Solution;

  function Recondition_Solution
             ( sol : DoblDobl_Complex_Solutions.Solution;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(sol.n+1);

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.v := Recondition_Solution_Vector(sol.v,k,s,locmap,xp);
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    return res;
  end Recondition_Solution;

  function Recondition_Solution
             ( sol : QuadDobl_Complex_Solutions.Solution;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(sol.n+1);

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.v := Recondition_Solution_Vector(sol.v,k,s,locmap,xp);
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    return res;
  end Recondition_Solution;

  function Recondition_Solutions
             ( sols : Standard_Complex_Solutions.Solution_List;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Recondition_Solution(ls.all,k,s,locmap,xp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Recondition_Solutions;

  function Recondition_Solutions
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Recondition_Solution(ls.all,k,s,locmap,xp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Recondition_Solutions;

  function Recondition_Solutions
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               k,s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Recondition_Solution(ls.all,k,s,locmap,xp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Recondition_Solutions;

  function Rescale_Solution_Vector
             ( x : Standard_Complex_Vectors.Vector; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'first..x'last-1);
    r,c : integer32;

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

  begin
    for k in res'range loop
      Checker_Localization_Patterns.Position(locmap,k,r,c);
      if c = s+1 then
        res(k) := x(k)/x(x'last);
      elsif c = s then
        if Degree(xp(r,c),k) = 0       -- x(ctr+1,s+1) does not occur
         then res(k) := x(k)/x(x'last);
         else res(k) := x(k);
        end if;
      else
        res(k) := x(k);
      end if;
    end loop;
    return res;
  end Rescale_Solution_Vector;

  function Rescale_Solution_Vector
             ( x : DoblDobl_Complex_Vectors.Vector; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(x'first..x'last-1);
    r,c : integer32;

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

  begin
    for k in res'range loop
      Checker_Localization_Patterns.Position(locmap,k,r,c);
      if c = s+1 then
        res(k) := x(k)/x(x'last);
      elsif c = s then
        if Degree(xp(r,c),k) = 0       -- x(ctr+1,s+1) does not occur
         then res(k) := x(k)/x(x'last);
         else res(k) := x(k);
        end if;
      else
        res(k) := x(k);
      end if;
    end loop;
    return res;
  end Rescale_Solution_Vector;

  function Rescale_Solution_Vector
             ( x : QuadDobl_Complex_Vectors.Vector; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(x'first..x'last-1);
    r,c : integer32;

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

  begin
    for k in res'range loop
      Checker_Localization_Patterns.Position(locmap,k,r,c);
      if c = s+1 then
        res(k) := x(k)/x(x'last);
      elsif c = s then
        if Degree(xp(r,c),k) = 0       -- x(ctr+1,s+1) does not occur
         then res(k) := x(k)/x(x'last);
         else res(k) := x(k);
        end if;
      else
        res(k) := x(k);
      end if;
    end loop;
    return res;
  end Rescale_Solution_Vector;

  function Rescale_Solution
             ( sol : Standard_Complex_Solutions.Solution; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(sol.n-1);

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.v := Rescale_Solution_Vector(sol.v,s,locmap,xp,pivot);
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    return res;
  end Rescale_Solution;

  function Rescale_Solution
             ( sol : DoblDobl_Complex_Solutions.Solution; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(sol.n-1);

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.v := Rescale_Solution_Vector(sol.v,s,locmap,xp,pivot);
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    return res;
  end Rescale_Solution;

  function Rescale_Solution
             ( sol : QuadDobl_Complex_Solutions.Solution; s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(sol.n-1);

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.v := Rescale_Solution_Vector(sol.v,s,locmap,xp,pivot);
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    return res;
  end Rescale_Solution;

  function Rescale_Solutions
             ( sols : Standard_Complex_Solutions.Solution_List;
               s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : Standard_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Rescale_Solution(ls.all,s,locmap,xp,pivot));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Rescale_Solutions;

  function Rescale_Solutions
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : DoblDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Rescale_Solution(ls.all,s,locmap,xp,pivot));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Rescale_Solutions;

  function Rescale_Solutions
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               s : integer32;
               locmap : Standard_Natural_Matrices.Matrix;
               xp : QuadDobl_Complex_Poly_Matrices.Matrix; pivot : integer32 )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Rescale_Solution(ls.all,s,locmap,xp,pivot));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Rescale_Solutions;

end Recondition_Swap_Homotopies;
