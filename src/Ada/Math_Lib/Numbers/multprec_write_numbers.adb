with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;

package body Multprec_Write_Numbers is

  function Is_Imag ( c : Complex_Number ) return boolean is

    rlp : Floating_Number := REAL_PART(c);
    res : constant boolean := Equal(rlp,0.0);
 
  begin
    Clear(rlp);
    return res;
  end Is_Imag;

  function Is_Real ( c : Complex_Number ) return boolean is

    imp : Floating_Number := IMAG_PART(c);
    res : constant boolean := Equal(imp,0.0);   

  begin
    Clear(imp);
    return res;
  end Is_Real;

  function Is_Positive_Real ( c : Complex_Number ) return boolean is

    imp : Floating_Number := IMAG_PART(c);
    res : boolean := Equal(imp,0.0);
    rlp : Floating_Number;

  begin
    if res then
      rlp := REAL_PART(c);
      res := (rlp > 0.0);
      Clear(rlp);
    end if;
    Clear(imp);
    return res;
  end Is_Positive_Real;

  function Is_Positive_Imag ( c : Complex_Number ) return boolean is

    rlp : Floating_Number := REAL_PART(c);
    res : boolean := Equal(rlp,0.0);
    imp : Floating_Number;

  begin
    if res then
      imp := IMAG_PART(c);
      res := (imp > 0.0);
      Clear(imp);
    end if;
    Clear(rlp);
    return res;
  end Is_Positive_Imag;

  procedure Write_Number ( file : in file_type; c : in Complex_Number ) is

    re,im : Floating_Number;

  begin
    if Is_Real(c) then
      re := REAL_PART(c);
      put(file,re); 
      Clear(re);
    elsif Is_Imag(c) then
      im := IMAG_PART(c);
      put(file,im);
      Clear(im);
      put(file,"*i");
    else
      re := REAL_PART(c);
      im := IMAG_PART(c);
      put(file,"(");
      put(file,re);
      if im > 0.0
       then put(file,"+");
      end if;
      put(file,im);
      put(file,"*i)");
      Clear(re); Clear(im);
    end if;
  end Write_Number;

  procedure Write_Coefficient ( file : in file_type; c : in Complex_Number ) is

    realzero : constant Floating_Number := Create(integer(0));
    realone : constant Floating_Number := Create(integer(1));
    realmin1 : constant Floating_Number := Create(integer(-1));
    compone : constant Complex_Number := Create(realone);
    compmin1 : constant Complex_Number := Create(realmin1);
    imagunit : constant Complex_Number := Create(realzero,realone);
    minimagu : constant Complex_Number := Create(realzero,realmin1);

  begin
    if Equal(c,compmin1) then
      put(file,'-');
    elsif Equal(c,imagunit) then
      put(file,"i*");
    elsif Equal(c,minimagu) then
      put(file,"-i*");
    elsif not Equal(c,compone) then
      Write_Number(file,c); put(file,'*');
    else -- the case where c = +1, nothing is written to file
      null;
    end if;
  end Write_Coefficient;

  procedure Write_Plus ( file : in file_type; c : in Complex_Number ) is
  begin
    if (Is_Positive_Real(c) or else Is_Positive_Imag(c))
      or else (not Is_Real(c) and not Is_Imag(c))
     then put(file,'+');
    end if;
  end Write_Plus;

end Multprec_Write_Numbers;
