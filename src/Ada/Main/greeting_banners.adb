package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.09 released 2016-01-02";

  begin
    return res;
  end Version;

end Greeting_Banners;
