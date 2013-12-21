with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Multprec_Complex_VecVecs;
with Multprec_Complex_VecVecs_io;        use Multprec_Complex_VecVecs_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;

package body Sample_Points_io is

  procedure get ( n,k : in integer32; sample : out Standard_Sample ) is
  begin
    get(Standard_Input,n,k,sample);
  end get;

  procedure get ( file : in file_type;
                  n,k : in integer32; sample : out Standard_Sample ) is

    c : character;
    sol : Standard_Complex_Solutions.Solution(n);
    hyp : Standard_Complex_VecVecs.VecVec(1..k);

  begin
    get(file,sol);
    get(file,c); skip_line(file);   -- skip status bar
    get(file,natural32(n+1),hyp);
    sample := Create(sol,hyp);
  end get;

  procedure get ( n,k : in integer32; sample : out Multprec_Sample ) is
  begin
    get(Standard_Input,n,k,sample);
  end get;

  procedure get ( file : in file_type;
                  n,k : in integer32; sample : out Multprec_Sample ) is

    c : character;
    sol : Multprec_Complex_Solutions.Solution(n);
    hyp : Multprec_Complex_VecVecs.VecVec(1..k);

  begin
    get(file,sol);
    get(file,c); skip_line(file);   -- skip status bar
    get(file,natural32(n+1),hyp);
    sample := Create(sol,hyp);
  end get;

  procedure put ( sample : in Standard_Sample ) is
  begin
    put(Standard_Output,sample);
  end put;

  procedure put ( file : in file_type; sample : in Standard_Sample ) is
  begin
    put(file,Sample_Point(sample));
    new_line(file);
    put_line(file,Hyperplane_Sections(sample));
  end put;

  procedure put ( sample : in Multprec_Sample ) is
  begin
    put(Standard_Output,sample);
  end put;

  procedure put ( file : in file_type; sample : in Multprec_Sample ) is
  begin
    put(file,Sample_Point(sample));
    new_line(file);
    put_line(file,Hyperplane_Sections(sample));
  end put;

end Sample_Points_io;
