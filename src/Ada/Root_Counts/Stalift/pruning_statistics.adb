with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

procedure Pruning_Statistics
              ( file : in file_type;
                nbsucc,nbfail : in Standard_Floating_Vectors.Vector ) is

  totsucc,totfail : double_float := 0.0;

begin
  new_line(file);
  put_line(file,"STATISTICS OF #FACE-FACE COMBINATIONS :");
  new_line(file);
  put_line(file,"           #Success      #Fail       Ratio ");
  new_line(file);
  for i in nbsucc'range loop
    put(file,"          ");
    put(file,nbsucc(i),2,3,3); put(file,"   ");
    put(file,nbfail(i),2,3,3); put(file,"   ");
    if nbsucc(i) + 1.0 /= 1.0
     then put(file,nbfail(i)/nbsucc(i),2,3,3);
     else put(file,"1/0");
    end if;
    new_line(file);
    totsucc := totsucc + nbsucc(i);
    totfail := totfail + nbfail(i);
  end loop;
  put_line(file," ----------------------------------------------");
  put(file,"  Total : "); 
    put(file,totsucc,2,3,3); put(file," + ");
    put(file,totfail,2,3,3); put(file," = ");
    put(file,totsucc+totfail,2,3,3); new_line(file);
  new_line(file);
  put_line(file,"  Success = successful face-face combinations");
  put_line(file,"  Fail    = unsuccessful face-face combinations");
  put_line(file,"  Ratio   = #Fail / #Success");
  put_line(file,"  Total   = total number of Linear-Programming problems");
end Pruning_Statistics;
