with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Evaluated_Minors;                   use Evaluated_Minors;

procedure Verify_Solution_Maps
            ( file : in file_type; pts : in Vector; planes : in VecMat;
              solmaps : in Array_of_Polynomial_Matrices;
              tol : in double_float; fail : out boolean ) is

  procedure Verify ( k : in integer32; pm : in Polynomial_Matrix ) is

    eva : Matrix(pm'range(1),pm'range(2));
    mat : Matrix(pm'range(1),pm'range(1));
    det : Complex_Number;

  begin
    for i in pts'range loop
      eva := Eval(pm,pts(i));
      put_line(file,"The evaluated matrix : "); put(file,eva);
      for j1 in eva'range(1) loop
        for j2 in eva'range(2) loop
          mat(j1,j2) := eva(j1,j2);
        end loop;
        for j2 in planes(i)'range(2) loop
          mat(j1,eva'last(2)+j2) := planes(i)(j1,j2);
        end loop;
      end loop;
      det := Determinant(mat);
      put(file," curve "); put(file,k,1);
      put(file," condition "); put(file,i,1); 
      put(file," determinant : ");
      put(file,det,3); new_line(file);
      if AbsVal(det) > tol
       then fail := true;
      end if;
    end loop;
  end Verify;

begin
  fail := false;
  for i in solmaps'range loop
    Verify(i,solmaps(i).all);
  end loop;
end Verify_Solution_Maps;
