with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Degrees_in_Sets_of_Unknowns;
with m_Homogeneous_Bezout_Numbers;
with Standard_Linear_Product_System;

package body m_homogeneous_permanent_factors is

  procedure Permanent
              ( row : in integer32;
                deg : in Standard_Integer_Matrices.Matrix;
                cols,crd : in out Standard_Integer_Vectors.Vector;
                per : in out integer32;
                ind,ind_last : in out Lists_of_integer_Vectors.List ) is

    use Lists_of_Integer_Vectors;
    acc : integer32;

  begin
    if row > deg'last(1) then
      acc := 1;
      for i in cols'range loop
        acc := acc*deg(i,cols(i));
      end loop;
      per := per + acc;
      Lists_of_Integer_Vectors.Append(ind,ind_last,cols);
    else
      for col in crd'range loop
        if crd(col) > 0 and deg(row,col) /= 0 then
          crd(col) := crd(col) - 1;
          cols(row) := col;
          Permanent(row+1,deg,cols,crd,per,ind,ind_last);
          crd(col) := crd(col) + 1;
        end if;
      end loop;
    end if;
  end Permanent;

  procedure Split_Indices
              ( ind : integer32;
                accu : in out Standard_Integer_Vectors.Vector;
                deg : in Standard_Integer_Matrices.Matrix;
                cols,base : in Standard_Integer_Vectors.Vector;
                res,res_last : in out Lists_of_Integer_Vectors.List ) is
  begin
    if ind > accu'last then
      Lists_of_Integer_Vectors.Append(res,res_last,accu);
    else
      for k in 1..deg(ind,cols(ind)) loop
        accu(ind) := base(ind) + k;
        Split_Indices(ind+1,accu,deg,cols,base,res,res_last);
      end loop;
    end if;
  end Split_Indices;

  function Split_Column_Indices
             ( deg : Standard_Integer_Matrices.Matrix;
               ind : Lists_of_Integer_Vectors.List )
             return Lists_of_Integer_Vectors.List is

    use Lists_of_Integer_Vectors;
    tmp : List := ind;
    lvc : Standard_Integer_Vectors.Link_to_Vector;
    res,res_last : List;
    accu,base : Standard_Integer_Vectors.Vector(deg'range(1));

  begin
    while not Is_Null(tmp) loop
      lvc := Head_Of(tmp);
      for i in lvc'range loop
        base(i) := 0;
        for k in deg'first(2)..lvc(i)-1 loop
          base(i) := base(i) + deg(i,k);
        end loop;
      end loop;
      Split_Indices(lvc'first,accu,deg,lvc.all,base,res,res_last);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Split_Column_Indices;

  procedure Permanent_Factors
              ( p : in Poly_Sys; z : in Partition;
                sols : out Lists_of_Integer_Vectors.List ) is

    k : Standard_Integer_Vectors.Vector(integer32(z'first)..integer32(z'last))
      := m_Homogeneous_Bezout_Numbers.Cardinalities(z);
    d : constant Standard_Integer_Matrices.Matrix
      := Degrees_in_Sets_of_Unknowns.Degree_Table(p,z);
    c : Standard_Integer_Vectors.Vector(d'range(1));
    prm : integer32 := 0;
    cols,cols_last : Lists_of_Integer_Vectors.List;

  begin
   -- put("The cardinalities : "); put(k); new_line;
   -- put_line("The degree table : "); put(d);
    Permanent(d'first(1),d,c,k,prm,cols,cols_last);
   -- put("The permanent of the degree table : "); put(prm,1); put_line(".");
   -- put_line("The list of column indices : "); put(cols);
   -- put("The length of the list of indices : ");
   -- put(Lists_of_Integer_Vectors.Length_Of(cols),1); put_line(".");
    sols := Split_Column_Indices(d,cols);
   -- put_line("The splitted list of indices : "); put(sols);
   -- put("The length of the splitted list of indices : ");
   -- put(Lists_of_Integer_Vectors.Length_Of(sols),1); put_line(".");
  end Permanent_Factors;

  function Solve_Linear_System
             ( ind : Standard_Integer_Vectors.Vector ) return Solution is

    res : Solution(ind'last);
    mat : Standard_Complex_Matrices.Matrix(ind'range,ind'range);
    rhs : Standard_Complex_Vectors.Vector(ind'range);
    cff : Standard_Complex_Vectors.Link_to_Vector;
    ipvt : Standard_Integer_Vectors.Vector(ind'range);
    info : integer32;

  begin
    for row in ind'range loop
      cff := Standard_Linear_Product_System.Get_Hyperplane
               (natural32(row),natural32(ind(row)));
      rhs(row) := -cff(0);
      for col in 1..cff'last loop
        mat(row,col) := cff(col);
      end loop;
    end loop; 
    lufac(mat,mat'last(1),ipvt,info);
    lusolve(mat,mat'last(1),ipvt,rhs);
    res.v := rhs;
    res.t := Create(0.0);
    res.m := 1;
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Solve_Linear_System; 

  procedure Solve_m_Homogeneous_Start_System
              ( ind_sols : in Lists_of_Integer_Vectors.List;
                q : out Poly_Sys; qsols : out Solution_List ) is

    use Lists_of_Integer_Vectors;
    tmp : List := ind_sols;
    lvc : Standard_Integer_Vectors.Link_to_Vector;
    qsols_last : Solution_List;

  begin
    q := Standard_Linear_Product_System.Polynomial_System;
    while not Is_Null(tmp) loop
      lvc := Head_Of(tmp);
      Append(qsols,qsols_last,Solve_Linear_System(lvc.all));
      tmp := Tail_Of(tmp);
    end loop;
  end Solve_m_Homogeneous_Start_System;

end m_homogeneous_permanent_factors;
