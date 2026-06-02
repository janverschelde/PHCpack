with Test_Real_Powered_Series;

package body Double_Real_Power_Series_IO is

  procedure Write ( x : in Double_rpSeries_Vectors.Vector ) is
  begin
    for i in x'range loop
      Test_Real_Powered_Series.Write(x(i).cff,x(i).pwt);
    end loop;
  end Write;

  procedure Write ( A : in Double_rpSeries_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        Test_Real_Powered_Series.Write(A(i,j).cff,A(i,j).pwt);
      end loop;
    end loop;
  end Write;

end Double_Real_Power_Series_IO;
