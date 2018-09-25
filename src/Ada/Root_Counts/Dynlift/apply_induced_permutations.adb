with Supports_of_Polynomial_Systems;
with Induced_Permutations;
with Black_Mixed_Volume_Computations;

package body Apply_Induced_Permutations is

  procedure Compute_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector ) is

    mix : Standard_Integer_Vectors.Vector(1..r);

    use Black_Mixed_Volume_Computations;

  begin
    for i in 1..r loop
      mix(i) := mtype(i-1);
    end loop;
    Make_Induced_Permutation(sup,mix,mcc,iprm);
  end Compute_Permutation;

  procedure Compute_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector ) is

    mix : Standard_Integer_Vectors.Vector(1..r);

    use Black_Mixed_Volume_Computations;

  begin
    for i in 1..r loop
      mix(i) := mtype(i-1);
    end loop;
    Make_Induced_Permutation(sup,stlb,mix,mcc,iprm);
  end Compute_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

-- WITH STABLE LIFTING BOUND :

  procedure Apply_Induced_Permutation
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,stlb,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,stlb,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

  procedure Apply_Induced_Permutation
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                stlb : in double_float; r : in integer32;
                perm,mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

    iprm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

    use Standard_Integer_Vectors;

  begin
    if perm /= null then
      sup := Supports_of_Polynomial_Systems.Create(p);
      Compute_Permutation(sup,stlb,r,perm,mtype,mcc,iprm);
      Induced_Permutations.Permute(iprm.all,p);
      Standard_Integer_Vectors.Clear(iprm);
      Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    end if;
  end Apply_Induced_Permutation;

end Apply_Induced_Permutations;
