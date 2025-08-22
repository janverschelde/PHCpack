package demics_global_constants is

-- DESCRIPTION :
--   Defines the global constants as in global.h,
--   translation initiated via g++ -c -fdump-ada-spec global.h.

  ITER_BLAND : constant := 25;  --  global.h:41
  ITER : constant := 1000;  --  global.h:42
  REINV : constant := 0;  --  global.h:43

  EXCESS : constant := 5;  --  global.h:45
  STRLENGTH : constant := 128;  --  global.h:46

  FALSE : constant := 0;  --  global.h:48
  TRUE : constant := 1;  --  global.h:49

  LENGTHOFTOKEN : constant := 128;  --  global.h:51
  PLUSZERO : constant := 1.0E-8;  --  global.h:52
  MINUSZERO : constant := -1.0E-8;  --  global.h:53
  BIGDOUBLE : constant := 1.0E+16;  --  global.h:54
  SMALLDOUBLE : constant := -1.0E+16;  --  global.h:55
  BIGINT : constant := 1000000000;  --  global.h:56

  FREE : constant := 2;  --  global.h:58
  NONNEGATIVE : constant := 3;  --  global.h:59
  POSITIVE : constant := 25;  --  global.h:60
  NEGATIVE : constant := 26;  --  global.h:61

  POSTHETA : constant := 6;  --  global.h:63
  NEGTHETA : constant := 7;  --  global.h:64

  OPT : constant := 4;  --  global.h:66
  UNBOUNDED : constant := 8;  --  global.h:67
  FEASIBLE : constant := 10;  --  global.h:68
  INFEASIBLE : constant := 11;  --  global.h:69

  CONTINUE : constant := 9;  --  global.h:71
  STOP : constant := 14;  --  global.h:72
  PIVOT_IN : constant := 28;  --  global.h:73

  ICHECK : constant := 20;  --  global.h:79
  MCHECK : constant := 21;  --  global.h:80

  TRIANGLE : constant := 28;  --  global.h:82
  SQUARE : constant := 29;  --  global.h:83

  SLIDE : constant := 16;  --  global.h:85
  DOWN : constant := 17;  --  global.h:86

  NODE : constant := 22;  --  global.h:88
  FN : constant := 23;  --  global.h:89
  FNN : constant := 24;  --  global.h:90

  ON : constant := 30;  --  global.h:92
  OFF : constant := 31;  --  global.h:93

  UNB_TAR : constant := 32;  --  global.h:95
  UNB_COR : constant := 33;  --  global.h:96

  ERROR_ITER : constant := 27;  --  global.h:98

-- for debugging ...

  DBG_INFO : constant := 0;  --  global.h:101

  DBG_NODE : constant := 0;  --  global.h:103
  DBG_FEA : constant := 0;  --  global.h:104
  DBG_SUC : constant := 0;  --  global.h:105
  DBG_TMP : constant := 0;  --  global.h:106
  DEBUG : constant := 0;  --  global.h:107

-- for relation table
  DBG_REL_TABLE : constant := 0;  --  global.h:110

-- for initial one-point test
  DBG_INI_CUR_INFO : constant := 0;  --  global.h:113

-- for one-point test
  DBG_PRE_INFO : constant := 0;  --  global.h:116
  DBG_CUR_INFO : constant := 0;  --  global.h:117

-- for suc-one-point test
  DBG_S_PRE_INFO : constant := 0;  --  global.h:120
  DBG_S_CUR_INFO : constant := 0;  --  global.h:121

-- for solLP
  DBG_TSOLLP : constant := 0;  --  global.h:124
  DBG_SOLLP : constant := 0;  --  global.h:125
  DBG_ISOLLP : constant := 0;  --  global.h:126
  DBG_MSOLLP : constant := 0;  --  global.h:127

-- for chooseSup
  DBG_CHOOSESUP : constant := 0;  --  global.h:130
  DBG_FINDUNB : constant := 0;  --  global.h:131

-- for modify p_sol
  DBG_MODIFY : constant := 0;  --  global.h:134

end demics_global_constants;
