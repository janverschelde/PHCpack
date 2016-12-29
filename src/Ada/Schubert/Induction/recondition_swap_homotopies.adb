with Standard_Natural_Vectors;
with Checker_Localization_Patterns;

package body Recondition_Swap_Homotopies is

  procedure Insert_One_Variable
              ( k : in integer32;
                t : in out Standard_Complex_Polynomials.Term ) is
  
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
              ( k : in integer32;
                p : in out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;

    procedure Insert_to_Term ( t : in out Term; c : out boolean ) is
    begin
      Insert_One_Variable(k,t);
      c := true;
    end Insert_to_Term;
    procedure Insert_to_Terms is new Changing_Iterator(Insert_to_Term);

  begin
    if p /= Null_Poly
     then Insert_to_Terms(p);
    end if;
  end Insert_One_Variable;

  procedure Insert_One_Variable
              ( k : in integer32;
                x : in out Standard_Complex_Poly_Matrices.Matrix ) is
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        Insert_One_Variable(k,x(i,j));
      end loop;
    end loop;
  end Insert_One_Variable;

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

  procedure Recondition
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                dim,s : in integer32 ) is

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);

  begin
    Insert_One_Variable(dim+1,x);
    Insert_Variable_Pivot(x,rowpiv,s+1,dim+1);
  end Recondition;

end Recondition_Swap_Homotopies;
