with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Poly_Matrices_io; use Standard_Complex_Poly_Matrices_io;
with Curves_into_Grassmannian;          use Curves_into_Grassmannian;
with Evaluated_Minors;                  use Evaluated_Minors;

package body Verification_with_Determinants is

  function Determinant
              ( m1,m2 : Standard_Complex_Matrices.Matrix )
              return Complex_Number is

    res : Complex_Number;
    ind : integer32;
    n1 : constant integer32 := m1'length(2);
    n2 : constant integer32 := m2'length(2);
    mat : Standard_Complex_Matrices.Matrix(m1'range(1),1..n1+n2);

  begin
    for i in mat'range(1) loop
      ind := mat'first(2);
      for j in m1'range(2) loop
        mat(i,ind) := m1(i,j);
        ind := ind + 1;
      end loop;
      for j in m2'range(2) loop
        mat(i,ind) := m2(i,j);
        ind := ind + 1;
      end loop;
    end loop;
    res := Determinant(mat);
    return res;
  end Determinant;

  procedure Determinant
              ( file : in file_type;
                m1,m2 : in Standard_Complex_Matrices.Matrix ) is

    d : constant Complex_Number := Determinant(m1,m2);

  begin
    put(file,"The determinant : "); put(file,d); new_line(file);
  end Determinant;

  procedure Verify_Determinants
              ( file : in file_type; dim : in integer32; nd : in Node;
                xpm : in Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                locsols : in Standard_Complex_VecVecs.VecVec;
                s : in Standard_Complex_Vectors.Vector; ip : in VecMat ) is

   -- addmix : constant natural := Mixed_Type(nd);
    i : constant integer32 := locsols'last;
   -- dim : natural := locsols(i)'last-2-addmix;
   -- dim : constant natural := locsols'last;
    lxpm : constant Standard_Complex_Poly_Matrices.Matrix
                       (xpm'range(1),xpm'range(2))
         := Column_Localize(nd.top,nd.bottom,locmap,xpm);
    xs : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2))
       := Substitute(lxpm,locsols(i)(1..dim));
    eva : Standard_Complex_Poly_Matrices.Matrix(xpm'range(1),xpm'range(2));
    mat : Standard_Complex_Matrices.Matrix(eva'range(1),eva'range(2));

  begin
   -- if dim >= locsols'last
   --  then dim := locsols'last;
   -- end if;
    put(file,"The dimension is "); put(file,dim,1); new_line(file);
    put(file,"The number of solutions : ");
    put(file,locsols'last,1); new_line(file);
    put_line(file,"The map on entry : "); put(file,xpm);
    put_line(file,"The map before substitution : "); put(file,lxpm);
    put_line(file,"The localized solutions : ");
    for j in locsols'range loop
      put(file,"solution "); put(file,j,1); put_line(file," :");
      put_line(file,locsols(j).all);
      xs := Substitute(lxpm,locsols(j)(1..dim));
      put_line("The map after substitution : "); put(file,xs);
      for k in 1..dim loop
        eva := Eval(xs,s(k),Create(1.0));
        put_line(file,"The evaluated plane : "); put(file,eva);
        mat := Convert(eva);
        put_line(file,"The plane as a matrix : "); put(file,mat);
        put(file,j,1); put(file," "); put(file,k,1); put(file," : ");
        Determinant(file,mat,ip(k).all);
      end loop;
    end loop;
  end Verify_Determinants;

  procedure Verify_Determinants
              ( file : in file_type; nd : in Node;
                xpm : in Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                locsols : in Solution_List;
                s : in Standard_Complex_Vectors.Vector; ip : in VecMat ) is

    len : constant integer32 := integer32(Length_Of(locsols));
    sols : Standard_Complex_VecVecs.VecVec(1..len);
    tmp : Solution_List := locsols;
    ls : Link_to_Solution;

  begin
    put(file,"The number of solution curves : ");
    put(file,len,1); new_line(file);
    for i in sols'range loop
      ls := Head_Of(tmp);
      sols(i) := new Standard_Complex_Vectors.Vector'(ls.v);
      tmp := Tail_Of(tmp);
    end loop;
    Verify_Determinants(file,Head_Of(locsols).n,nd,xpm,locmap,sols,s,ip);
    Standard_Complex_VecVecs.Clear(sols);
  end Verify_Determinants;

end Verification_with_Determinants;
