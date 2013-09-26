with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Minkowski_Polynomials;              use Minkowski_Polynomials;

procedure Driver_for_Minkowski_Polynomials
                ( file : in file_type;
                  n : in integer32; mix : in Vector; t : in Triangulation;
                  alltri : in boolean; mixsub : out Mixed_Subdivision ) is

  mp : Poly;
  timer : timing_widget;

  procedure Coefficient_Subdivision
                 ( submix : in Vector; sub : in Mixed_Subdivision;
                   vol : out natural32 ) is

    len : constant integer32 := integer32(Length_Of(sub));
    wrksub : Mixed_Subdivision := sub;
    volu : natural32;

  begin
    put(file,"#Cells of type "); put(file,submix); put(file," is " );
    put(file,len,1); new_line(file);
    if len = 0
     then volu := 0;
     else put(file,natural32(n),submix,wrksub,volu);
    end if;
    put(file,"The volume equals : "); put(file,volu,1); new_line(file);
    new_line(file);
    vol := volu;
  end Coefficient_Subdivision;
  procedure Coefficient_Subdivisions is 
    new Minkowski_Polynomial_Subdivisions(Coefficient_Subdivision);

  procedure Write_Minkowski_Polynomial ( file : in file_type; p : in Poly ) is

    first : boolean;
    cnt : natural32 := 0;

    procedure Minkowski_Term ( t : in Term; cont : out boolean ) is
    begin
      if first
       then first := false;
       else put(file," + "); cnt := cnt + 3;
      end if;
      put(file,integer32(REAL_PART(t.cf)),1); cnt := cnt + 2;
      for i in t.dg'range loop
        if t.dg(i) /= 0 then
          put(file,"*l"); put(file,i,1); cnt := cnt + 3;
          if t.dg(i) /= 1 then
            put(file,"^"); put(file,t.dg(i),1);
            cnt := cnt + 2;
          end if;
        end if;
      end loop;
      if cnt > 60
       then new_line(file); cnt := 0;
      end if;
      cont := true;
    end Minkowski_Term;
    procedure Minkowski_Terms is new Visiting_Iterator(Minkowski_Term);

  begin
    Minkowski_Terms(p);
  end Write_Minkowski_Polynomial;

begin
  tstart(timer);
  mp := Minkowski_Polynomial(natural32(n),natural32(mix'last));
  if not alltri
   then Minkowski_Polynomial(mp,t,natural32(n),mix,mixsub);
   else Coefficient_Subdivisions(mp,t,natural32(n),mix,mixsub);
  end if;
  tstop(timer);
  put_line(file,"the Minkowski-polynomial : ");
  Write_Minkowski_Polynomial(file,mp);
  new_line(file);
  new_line(file);
  print_times(file,timer,"computing the Minkowski-polynomial");
  Clear(mp);
end Driver_for_Minkowski_Polynomials;
