# This makefile is organized along the library structure of PHCpack.
# For every library, there are a number of executables (to test the
# modules in the library) which can be made.

# swap the # symbol in the following lines to use makefile_windows
# or to use makefile_mac on a Mac OS x computer

# MAKEFILE = makefile_unix
# MAKEFILE = makefile_mac
MAKEFILE = makefile_windows

# the main executable is phc

phc:
	make -f $(MAKEFILE) phc

phc_noqd:
	make -f $(MAKEFILE) phc_noqd

# for phcpy for Python 2.6 and 2.7, make phcpy2c2.so

phcpy2c2.so:
	make -f $(MAKEFILE) phcpy2c2.so

# for use with Sage :

sage_phcpy2c2.so:
	make -f $(MAKEFILE) sage_phcpy2c2.so

# for phcpy on Python 3.5, make phcpy2c3.so
# of for 3.6 and 3.7, do make phcpy2c36.so and make phcpy2c37.so

phcpy2c3.so:
	make -f $(MAKEFILE) phcpy2c3.so

phcpy2c36.so:
	make -f $(MAKEFILE) phcpy2c36.so

phcpy2c37.so:
	make -f $(MAKEFILE) phcpy2c37.so

# without the dependencies on GNAT, on Linux:

phcpy2c2.a:
	make -f $(MAKEFILE) phcpy2c2.a

# with the double version of the Path library :

phcpy2cpath_d.so:
	make -f $(MAKEFILE) phcpy2cpath_d.so

# all interactive test programs, lib programs and MPI executables :

testall: test_system test_math_lib test_deformations \
         test_root_counts test_schubert test_components \
         test_main test_ctophc test_lib test_tasking test_mpi

# Test program to see if compiler works, "hello world" :

hello:
	make -f $(MAKEFILE) hello

ada_hello:
	make -f $(MAKEFILE) ada_hello

# Test programs for the System library :

test_system: ts_time ts_clock ts_stamps ts_timer ts_ctimer ts_syscall \
             ts_mach ts_cmdline ts_byemess ts_strsplit ts_commuser \
             ts_filescan ts_striscan ts_threads

ts_time:
	make -f $(MAKEFILE) ts_time

ts_clock:
	make -f $(MAKEFILE) ts_clock

ts_stamps:
	make -f $(MAKEFILE) ts_stamps

ts_timer:
	make -f $(MAKEFILE) ts_timer

ts_ctimer:
	make -f $(MAKEFILE) ts_ctimer

ts_syscall:
	make -f $(MAKEFILE) ts_syscall

ts_mach:
	make -f $(MAKEFILE) ts_mach

ts_cmdline:
	make -f $(MAKEFILE) ts_cmdline

ts_byemess:
	make -f $(MAKEFILE) ts_byemess

ts_strsplit:
	make -f $(MAKEFILE) ts_strsplit

ts_commuser:
	make -f $(MAKEFILE) ts_commuser

ts_filescan:
	make -f $(MAKEFILE) ts_filescan

ts_striscan:
	make -f $(MAKEFILE) ts_striscan

ts_threads:
	make -f $(MAKEFILE) ts_threads

# Test programs for the library Math_Lib :

test_math_lib: test_numbers test_qd test_vectors test_matrices test_divisors \
               test_reduction test_polynomials test_functions test_supports \
               test_circuits test_series test_laurent test_ad

# Test programs for the library Math_Lib/Numbers :

test_numbers: ts_strnum ts_printf ts_natnum ts_intnum ts_fltnum ts_cmpnum \
              ts_natdiv ts_longflt ts_plrexp ts_random ts_extran ts_matfun \
              ts_mpcst ts_numbio ts_intval

ts_strnum:
	make -f $(MAKEFILE) ts_strnum

ts_printf:
	make -f $(MAKEFILE) ts_printf

ts_intnum:
	make -f $(MAKEFILE) ts_intnum

ts_natnum:
	make -f $(MAKEFILE) ts_natnum

ts_natdiv:
	make -f $(MAKEFILE) ts_natdiv

ts_longflt:
	make -f $(MAKEFILE) ts_longflt

ts_fltnum:
	make -f $(MAKEFILE) ts_fltnum

ts_cmpnum:
	make -f $(MAKEFILE) ts_cmpnum

ts_plrexp:
	make -f $(MAKEFILE) ts_plrexp

ts_random:
	make -f $(MAKEFILE) ts_random

ts_extran:
	make -f $(MAKEFILE) ts_extran

ts_matfun:
	make -f $(MAKEFILE) ts_matfun

ts_intval:
	make -f $(MAKEFILE) ts_intval

ts_mpcst:
	make -f $(MAKEFILE) ts_mpcst

ts_numbio:
	make -f $(MAKEFILE) ts_numbio

# Test programs for Math_Lib/QD :

test_qd: ts_dbldbl ts_qddbl ts_ddcmplx ts_qdcmplx ts_mpddcv ts_mpqdcv \
         ts_ddfun ts_qdfun ts_cmpexp ts_tridbl ts_octdbl ts_tdcmplx \
         ts_odcmplx ts_pentadbl ts_decadbl ts_pdcmplx ts_dacmplx ts_tdfun \
	 ts_pdfun ts_odfun ts_dafun ts_mpcasts

ts_dbldbl:
	make -f $(MAKEFILE) ts_dbldbl

ts_qddbl:
	make -f $(MAKEFILE) ts_qddbl

ts_ddcmplx:
	make -f $(MAKEFILE) ts_ddcmplx

ts_qdcmplx:
	make -f $(MAKEFILE) ts_qdcmplx

ts_mpddcv:
	make -f $(MAKEFILE) ts_mpddcv

ts_mpqdcv:
	make -f $(MAKEFILE) ts_mpqdcv

ts_ddfun:
	make -f $(MAKEFILE) ts_ddfun

ts_qdfun:
	make -f $(MAKEFILE) ts_qdfun

ts_tdfun:
	make -f $(MAKEFILE) ts_tdfun

ts_pdfun:
	make -f $(MAKEFILE) ts_pdfun

ts_odfun:
	make -f $(MAKEFILE) ts_odfun

ts_dafun:
	make -f $(MAKEFILE) ts_dafun

ts_cmpexp:
	make -f $(MAKEFILE) ts_cmpexp

ts_isnan:
	make -f $(MAKEFILE) ts_isnan

ts_tridbl:
	make -f $(MAKEFILE) ts_tridbl

ts_tdcmplx:
	make -f $(MAKEFILE) ts_tdcmplx

ts_octdbl:
	make -f $(MAKEFILE) ts_octdbl

ts_odcmplx:
	make -f $(MAKEFILE) ts_odcmplx

ts_pentadbl:
	make -f $(MAKEFILE) ts_pentadbl

ts_decadbl:
	make -f $(MAKEFILE) ts_decadbl

ts_hexdbl:
	make -f $(MAKEFILE) ts_hexdbl

ts_hdfun:
	make -f $(MAKEFILE) ts_hdfun

ts_pdcmplx:
	make -f $(MAKEFILE) ts_pdcmplx

ts_dacmplx:
	make -f $(MAKEFILE) ts_dacmplx

ts_hdcmplx:
	make -f $(MAKEFILE) ts_hdcmplx

ts_mpcasts:
	make -f $(MAKEFILE) ts_mpcasts

# Test programs for Math_Lib/Vectors :

test_vectors: ts_natvec ts_intvec ts_fltvec ts_cmpvec ts_ddvec ts_qdvec \
              ts_nesvec ts_strvec ts_perfddvc ts_perfqdvc ts_tdvec ts_pdvec \
	      ts_odvec ts_davec ts_hdvec

ts_natvec:
	make -f $(MAKEFILE) ts_natvec

ts_intvec:
	make -f $(MAKEFILE) ts_intvec

ts_fltvec:
	make -f $(MAKEFILE) ts_fltvec

ts_cmpvec:
	make -f $(MAKEFILE) ts_cmpvec

ts_nesvec:
	make -f $(MAKEFILE) ts_nesvec

ts_ddvec:
	make -f $(MAKEFILE) ts_ddvec

ts_qdvec:
	make -f $(MAKEFILE) ts_qdvec

ts_strvec:
	make -f $(MAKEFILE) ts_strvec

ts_perfddvc:
	make -f $(MAKEFILE) ts_perfddvc

ts_perfqdvc:
	make -f $(MAKEFILE) ts_perfqdvc

ts_tdvec:
	make -f $(MAKEFILE) ts_tdvec

ts_pdvec:
	make -f $(MAKEFILE) ts_pdvec

ts_odvec:
	make -f $(MAKEFILE) ts_odvec

ts_davec:
	make -f $(MAKEFILE) ts_davec

ts_hdvec:
	make -f $(MAKEFILE) ts_hdvec

# Test programs for Math_Lib/Matrices :

test_matrices: ts_natmat ts_intmat ts_fltmat ts_fltdls ts_cmpmat \
               ts_ddmat ts_qdmat ts_tdmat ts_pdmat ts_odmat ts_damat ts_hdmat

ts_natmat:
	make -f $(MAKEFILE) ts_natmat

ts_intmat:
	make -f $(MAKEFILE) ts_intmat

ts_fltmat:
	make -f $(MAKEFILE) ts_fltmat

ts_fltdls:
	make -f $(MAKEFILE) ts_fltdls

ts_cmpmat:
	make -f $(MAKEFILE) ts_cmpmat

ts_ddmat:
	make -f $(MAKEFILE) ts_ddmat

ts_qdmat:
	make -f $(MAKEFILE) ts_qdmat

ts_tdmat:
	make -f $(MAKEFILE) ts_tdmat

ts_pdmat:
	make -f $(MAKEFILE) ts_pdmat

ts_odmat:
	make -f $(MAKEFILE) ts_odmat

ts_damat:
	make -f $(MAKEFILE) ts_damat

ts_hdmat:
	make -f $(MAKEFILE) ts_hdmat

# Test programs for Math_Lib/Divisors :

test_divisors: ts_gcd ts_intlina ts_smith ts_intinv ts_circuit

ts_gcd:
	make -f $(MAKEFILE) ts_gcd

ts_intlina:
	make -f $(MAKEFILE) ts_intlina

ts_smith:
	make -f $(MAKEFILE) ts_smith

ts_intinv:
	make -f $(MAKEFILE) ts_intinv

ts_circuit:
	make -f $(MAKEFILE) ts_circuit

# Test programs for Math_Lib/Reduction :

test_reduction: ts_fltlu ts_cmplu ts_ddlu ts_qdlu ts_qrd ts_svd ts_eigval \
                ts_mgs ts_mgsrnk ts_backsubs ts_perflu ts_vvlu ts_vmplu \
                ts_perfdlu ts_perfdsvd ts_tdlu ts_pdlu ts_odlu ts_dalu \
                ts_ddnorm ts_tdnorm ts_qdnorm ts_pdnorm ts_odnorm ts_danorm \
                ts_dqrd ts_ddqrd ts_tdqrd ts_qdqrd ts_pdqrd ts_odqrd ts_daqrd \
                ts_mpqrd ts_dsvd ts_ddsvd ts_tdsvd ts_qdsvd ts_pdsvd ts_odsvd \
                ts_dasvd ts_mpsvd ts_errfree ts_hdlu ts_hdqrd ts_hdsvd

ts_fltlu:
	make -f $(MAKEFILE) ts_fltlu

ts_cmplu:
	make -f $(MAKEFILE) ts_cmplu

ts_ddlu:
	make -f $(MAKEFILE) ts_ddlu

ts_qdlu:
	make -f $(MAKEFILE) ts_qdlu

ts_qrd:
	make -f $(MAKEFILE) ts_qrd

ts_dqrd:
	make -f $(MAKEFILE) ts_dqrd

ts_ddqrd:
	make -f $(MAKEFILE) ts_ddqrd

ts_tdqrd:
	make -f $(MAKEFILE) ts_tdqrd

ts_qdqrd:
	make -f $(MAKEFILE) ts_qdqrd

ts_pdqrd:
	make -f $(MAKEFILE) ts_pdqrd

ts_odqrd:
	make -f $(MAKEFILE) ts_odqrd

ts_daqrd:
	make -f $(MAKEFILE) ts_daqrd

ts_hdqrd:
	make -f $(MAKEFILE) ts_hdqrd

ts_mpqrd:
	make -f $(MAKEFILE) ts_mpqrd

ts_svd:
	make -f $(MAKEFILE) ts_svd

ts_dsvd:
	make -f $(MAKEFILE) ts_dsvd

ts_ddsvd:
	make -f $(MAKEFILE) ts_ddsvd

ts_tdsvd:
	make -f $(MAKEFILE) ts_tdsvd

ts_qdsvd:
	make -f $(MAKEFILE) ts_qdsvd

ts_pdsvd:
	make -f $(MAKEFILE) ts_pdsvd

ts_odsvd:
	make -f $(MAKEFILE) ts_odsvd

ts_dasvd:
	make -f $(MAKEFILE) ts_dasvd

ts_hdsvd:
	make -f $(MAKEFILE) ts_hdsvd

ts_mpsvd:
	make -f $(MAKEFILE) ts_mpsvd

ts_rowred:
	make -f $(MAKEFILE) ts_rowred

ts_eigval:
	make -f $(MAKEFILE) ts_eigval

ts_mgs:
	make -f $(MAKEFILE) ts_mgs

ts_mgsrnk:
	make -f $(MAKEFILE) ts_mgsrnk

ts_backsubs:
	make -f $(MAKEFILE) ts_backsubs

ts_perflu:
	make -f $(MAKEFILE) ts_perflu

ts_perfdlu:
	make -f $(MAKEFILE) ts_perfdlu

ts_perfdsvd:
	make -f $(MAKEFILE) ts_perfdsvd

ts_vvlu:
	make -f $(MAKEFILE) ts_vvlu

ts_vmplu:
	make -f $(MAKEFILE) ts_vmplu

ts_tdlu:
	make -f $(MAKEFILE) ts_tdlu

ts_pdlu:
	make -f $(MAKEFILE) ts_pdlu

ts_odlu:
	make -f $(MAKEFILE) ts_odlu

ts_dalu:
	make -f $(MAKEFILE) ts_dalu

ts_hdlu:
	make -f $(MAKEFILE) ts_hdlu

ts_ddnorm:
	make -f $(MAKEFILE) ts_ddnorm

ts_tdnorm:
	make -f $(MAKEFILE) ts_tdnorm

ts_qdnorm:
	make -f $(MAKEFILE) ts_qdnorm

ts_pdnorm:
	make -f $(MAKEFILE) ts_pdnorm

ts_odnorm:
	make -f $(MAKEFILE) ts_odnorm

ts_danorm:
	make -f $(MAKEFILE) ts_danorm

ts_errfree:
	make -f $(MAKEFILE) ts_errfree

# Test programs for Math_Lib/Polynomials :

test_polynomials: ts_poly ts_laurpoly ts_laurlist ts_rpol ts_ddpol ts_qdpol \
                  ts_polmat ts_jaco ts_parse ts_subs ts_tabform ts_randpoly \
                  ts_str2pol ts_termlist ts_rwspol ts_nbrvar ts_hessian \
		  ts_tdpol ts_pdpol ts_odpol ts_dapol

ts_poly:
	make -f $(MAKEFILE) ts_poly

ts_laurpoly:
	make -f $(MAKEFILE) ts_laurpoly

ts_laurlist:
	make -f $(MAKEFILE) ts_laurlist

ts_rpol:
	make -f $(MAKEFILE) ts_rpol

ts_ddpol:
	make -f $(MAKEFILE) ts_ddpol

ts_tdpol:
	make -f $(MAKEFILE) ts_tdpol

ts_qdpol:
	make -f $(MAKEFILE) ts_qdpol

ts_pdpol:
	make -f $(MAKEFILE) ts_pdpol

ts_odpol:
	make -f $(MAKEFILE) ts_odpol

ts_dapol:
	make -f $(MAKEFILE) ts_dapol

ts_polmat:
	make -f $(MAKEFILE) ts_polmat

ts_jaco:
	make -f $(MAKEFILE) ts_jaco

ts_parse:
	make -f $(MAKEFILE) ts_parse

ts_subs:
	make -f $(MAKEFILE) ts_subs

ts_tabform:
	make -f $(MAKEFILE) ts_tabform

ts_randpoly:
	make -f $(MAKEFILE) ts_randpoly

ts_str2pol:
	make -f $(MAKEFILE) ts_str2pol

ts_termlist:
	make -f $(MAKEFILE) ts_termlist

ts_rwspol:
	make -f $(MAKEFILE) ts_rwspol

ts_nbrvar:
	make -f $(MAKEFILE) ts_nbrvar

ts_hessian:
	make -f $(MAKEFILE) ts_hessian

# Test programs for Math_Lib/Functions :

test_functions: ts_evaline ts_evalpoly ts_evalform ts_lineva \
                ts_expvec ts_polyflat ts_perfeval ts_evddpol

ts_evaline:
	make -f $(MAKEFILE) ts_evaline

ts_evalpoly:
	make -f $(MAKEFILE) ts_evalpoly

ts_evalform:
	make -f $(MAKEFILE) ts_evalform

ts_lineva:
	make -f $(MAKEFILE) ts_lineva

ts_polyflat:
	make -f $(MAKEFILE) ts_polyflat

ts_perfeval:
	make -f $(MAKEFILE) ts_perfeval

ts_expvec:
	make -f $(MAKEFILE) ts_expvec

ts_evddpol:
	make -f $(MAKEFILE) ts_evddpol

# Test programs for Math_Lib/Circuits :

test_circuits: ts_cffsup ts_speel ts_cycfun ts_speelsys ts_gradcirc \
               ts_vmpeval ts_vmpdiff ts_speelser ts_speelcnv ts_fabry \
               ts_evalcnv ts_hesscnv ts_jacocnv ts_shiftcnv ts_perfconv \
               ts_perfcirc ts_perfade

ts_cffsup:
	make -f $(MAKEFILE) ts_cffsup

ts_gradcirc:
	make -f $(MAKEFILE) ts_gradcirc

ts_speel:
	make -f $(MAKEFILE) ts_speel

ts_speelsys:
	make -f $(MAKEFILE) ts_speelsys

ts_cycfun:
	make -f $(MAKEFILE) ts_cycfun

ts_vmpeval:
	make -f $(MAKEFILE) ts_vmpeval

ts_vmpdiff:
	make -f $(MAKEFILE) ts_vmpdiff

ts_speelser:
	make -f $(MAKEFILE) ts_speelser

ts_speelcnv:
	make -f $(MAKEFILE) ts_speelcnv

ts_fabry:
	make -f $(MAKEFILE) ts_fabry

ts_evalcnv:
	make -f $(MAKEFILE) ts_evalcnv

ts_hesscnv:
	make -f $(MAKEFILE) ts_hesscnv

ts_jacocnv:
	make -f $(MAKEFILE) ts_jacocnv

ts_shiftcnv:
	make -f $(MAKEFILE) ts_shiftcnv

ts_perfconv:
	make -f $(MAKEFILE) ts_perfconv

ts_perfcirc:
	make -f $(MAKEFILE) ts_perfcirc

ts_perfade:
	make -f $(MAKEFILE) ts_perfade

ts_perfhess:
	make -f $(MAKEFILE) ts_perfhess

# Test programs for Math_Lib/Supports :

test_supports: ts_basex ts_feasi ts_enurs ts_diclp ts_givrot ts_lstvec \
               ts_fvector ts_giftwrap

ts_basex:
	make -f $(MAKEFILE) ts_basex

ts_feasi:
	make -f $(MAKEFILE) ts_feasi

ts_enurs:
	make -f $(MAKEFILE) ts_enurs

ts_diclp:
	make -f $(MAKEFILE) ts_diclp

ts_givrot:
	make -f $(MAKEFILE) ts_givrot

ts_lstvec:
	make -f $(MAKEFILE) ts_lstvec

ts_fvector:
	make -f $(MAKEFILE) ts_fvector

ts_giftwrap:
	make -f $(MAKEFILE) ts_giftwrap

# Test programs for Math_Lib/Series :

test_series: ts_pstrunc ts_bincff ts_series ts_series2 ts_cseries \
             ts_servec ts_servec2 ts_cservec ts_sermat ts_sermat2 ts_csermat \
             ts_seritp ts_seritp2 ts_cseritp ts_sersin ts_sersin2 ts_csersin \
             ts_serlin ts_serlin2 ts_cserlin ts_serinv ts_cserinv \
             ts_serpol ts_cserpol ts_sersys ts_csersys ts_sercffpol \
             ts_perfserlin

ts_pstrunc:
	make -f $(MAKEFILE) ts_pstrunc

ts_bincff:
	make -f $(MAKEFILE) ts_bincff

ts_series:
	make -f $(MAKEFILE) ts_series

ts_series2:
	make -f $(MAKEFILE) ts_series2

ts_cseries:
	make -f $(MAKEFILE) ts_cseries

ts_servec:
	make -f $(MAKEFILE) ts_servec

ts_servec2:
	make -f $(MAKEFILE) ts_servec2

ts_cservec:
	make -f $(MAKEFILE) ts_cservec

ts_sermat:
	make -f $(MAKEFILE) ts_sermat

ts_sermat2:
	make -f $(MAKEFILE) ts_sermat2

ts_csermat:
	make -f $(MAKEFILE) ts_csermat

ts_seritp:
	make -f $(MAKEFILE) ts_seritp

ts_seritp2:
	make -f $(MAKEFILE) ts_seritp2

ts_cseritp:
	make -f $(MAKEFILE) ts_cseritp

ts_sersin:
	make -f $(MAKEFILE) ts_sersin

ts_sersin2:
	make -f $(MAKEFILE) ts_sersin2

ts_csersin:
	make -f $(MAKEFILE) ts_csersin

ts_serlin:
	make -f $(MAKEFILE) ts_serlin

ts_serlin2:
	make -f $(MAKEFILE) ts_serlin2

ts_cserlin:
	make -f $(MAKEFILE) ts_cserlin

ts_perfserlin:
	make -f $(MAKEFILE) ts_perfserlin

ts_serinv:
	make -f $(MAKEFILE) ts_serinv

ts_cserinv:
	make -f $(MAKEFILE) ts_cserinv

ts_serpol:
	make -f $(MAKEFILE) ts_serpol

ts_cserpol:
	make -f $(MAKEFILE) ts_cserpol

ts_sersys:
	make -f $(MAKEFILE) ts_sersys

ts_csersys:
	make -f $(MAKEFILE) ts_csersys

ts_sercffpol:
	make -f $(MAKEFILE) ts_sercffpol

# Test programs for the Math_Lib/Laurent library :

test_laurent: ts_pslaur ts_laurmat ts_lserpol ts_lsernew

ts_pslaur:
	make -f $(MAKEFILE) ts_pslaur

ts_laurmat:
	make -f $(MAKEFILE) ts_laurmat

ts_lserpol:
	make -f $(MAKEFILE) ts_lserpol

ts_lsernew:
	make -f $(MAKEFILE) ts_lsernew

# Test programs for the Math_Lib/AD library :

test_ad: ts_monom ts_monvec ts_polvec ts_sys2vec

ts_monom:
	make -f $(MAKEFILE) ts_monom

ts_monvec:
	make -f $(MAKEFILE) ts_monvec

ts_polvec:
	make -f $(MAKEFILE) ts_polvec

ts_sys2vec:
	make -f $(MAKEFILE) ts_sys2vec

# Test programs for the Deformations library :

test_deformations: test_solutions test_homotopy test_newton test_curves \
                   test_end_games test_trackers test_continuation

# Test programs for the Deformations/Solutions library :

test_solutions: ts_sols ts_ddsols ts_qdsols ts_solfilt ts_solar ts_solstr \
                ts_mapsols_io ts_dicsols_io ts_grepsols ts_selsol \
                ts_ptlists ts_quadpart ts_contab ts_quad ts_soldrop \
                ts_trasols_io ts_mixres ts_getstart

ts_sols:
	make -f $(MAKEFILE) ts_sols

ts_ddsols:
	make -f $(MAKEFILE) ts_ddsols

ts_qdsols:
	make -f $(MAKEFILE) ts_qdsols

ts_solfilt:
	make -f $(MAKEFILE) ts_solfilt

ts_solar:
	make -f $(MAKEFILE) ts_solar

ts_mapsols_io:
	make -f $(MAKEFILE) ts_mapsols_io

ts_dicsols_io:
	make -f $(MAKEFILE) ts_dicsols_io

ts_grepsols:
	make -f $(MAKEFILE) ts_grepsols

ts_selsol:
	make -f $(MAKEFILE) ts_selsol

ts_solstr:
	make -f $(MAKEFILE) ts_solstr

ts_contab:
	make -f $(MAKEFILE) ts_contab

ts_quadpart:
	make -f $(MAKEFILE) ts_quadpart

ts_quad:
	make -f $(MAKEFILE) ts_quad

ts_soldrop:
	make -f $(MAKEFILE) ts_soldrop

ts_ptlists:
	make -f $(MAKEFILE) ts_ptlists

ts_trasols_io:
	make -f $(MAKEFILE) ts_trasols_io

ts_mixres:
	make -f $(MAKEFILE) ts_mixres

ts_getstart:
	make -f $(MAKEFILE) ts_getstart

# Test programs for the Deformations/Homotopy library :

test_homotopy: ts_quad ts_mainscal ts_scal ts_reduce ts_homotopy \
               ts_homline ts_linsol ts_qdhom ts_poldrop ts_evalhomt ts_vmphom \
               ts_scalplane ts_multproj ts_homcnv

ts_mainscal:
	make -f $(MAKEFILE) ts_mainscal

ts_scal:
	make -f $(MAKEFILE) ts_scal

ts_reduce:
	make -f $(MAKEFILE) ts_reduce

ts_homotopy:
	make -f $(MAKEFILE) ts_homotopy

ts_homline:
	make -f $(MAKEFILE) ts_homline

ts_linsol:
	make -f $(MAKEFILE) ts_linsol

ts_qdhom:
	make -f $(MAKEFILE) ts_qdhom

ts_poldrop:
	make -f $(MAKEFILE) ts_poldrop

ts_evalhomt:
	make -f $(MAKEFILE) ts_evalhomt

ts_vmphom:
	make -f $(MAKEFILE) ts_vmphom

ts_scalplane:
	make -f $(MAKEFILE) ts_scalplane

ts_multproj:
	make -f $(MAKEFILE) ts_multproj

ts_homcnv:
	make -f $(MAKEFILE) ts_homcnv

# Test programs for the Deformations/Newton library :

test_newton: ts_mreseva ts_rootrefi ts_newton ts_realnewt ts_numdif \
             ts_deftrees ts_jactrees ts_defmat ts_deflate ts_multip ts_stair \
             ts_deflop ts_numdif ts_straight ts_ddnewt ts_qdnewt ts_vmpnewt \
	     ts_newcirc

ts_mreseva:
	make -f $(MAKEFILE) ts_mreseva

ts_rootrefi:
	make -f $(MAKEFILE) ts_rootrefi

ts_newton:
	make -f $(MAKEFILE) ts_newton

ts_realnewt:
	make -f $(MAKEFILE) ts_realnewt

ts_numdif:
	make -f $(MAKEFILE) ts_numdif

ts_jactrees:
	make -f $(MAKEFILE) ts_jactrees

ts_deftrees:
	make -f $(MAKEFILE) ts_deftrees

ts_deflate:
	make -f $(MAKEFILE) ts_deflate

ts_straight:
	make -f $(MAKEFILE) ts_straight

ts_stair:
	make -f $(MAKEFILE) ts_stair

ts_defmat:
	make -f $(MAKEFILE) ts_defmat

ts_multip:
	make -f $(MAKEFILE) ts_multip

ts_muldef:
	make -f $(MAKEFILE) ts_muldef

ts_deflop:
	make -f $(MAKEFILE) ts_deflop

ts_ddnewt:
	make -f $(MAKEFILE) ts_ddnewt

ts_qdnewt:
	make -f $(MAKEFILE) ts_qdnewt

ts_vmpnewt:
	make -f $(MAKEFILE) ts_vmpnewt

ts_newcirc:
	make -f $(MAKEFILE) ts_newcirc

# Test programs for the Deformations/Curves library :

test_curves: ts_durker ts_extrapol ts_extrapts \
             ts_sernew ts_sersol ts_serhom ts_serpred ts_serpath ts_serpade \
             ts_padepred ts_jacrabin ts_nxtpadsol ts_hesscrit ts_hesspred \
             ts_predstep ts_sernewcnv ts_ratapp ts_padepcnv ts_hesspcnv \
             ts_corcnv ts_pcscnv ts_scalecnv ts_wrappade ts_fabryhom

ts_durker:
	make -f $(MAKEFILE) ts_durker

ts_extrapol:
	make -f $(MAKEFILE) ts_extrapol

ts_extrapts:
	make -f $(MAKEFILE) ts_extrapts

ts_sernew:
	make -f $(MAKEFILE) ts_sernew

ts_sersol:
	make -f $(MAKEFILE) ts_sersol

ts_serhom:
	make -f $(MAKEFILE) ts_serhom

ts_serpred:
	make -f $(MAKEFILE) ts_serpred

ts_serpath:
	make -f $(MAKEFILE) ts_serpath

ts_serpade:
	make -f $(MAKEFILE) ts_serpade

ts_padepred:
	make -f $(MAKEFILE) ts_padepred

ts_padepcnv:
	make -f $(MAKEFILE) ts_padepcnv

ts_jacrabin:
	make -f $(MAKEFILE) ts_jacrabin

ts_nxtpadsol:
	make -f $(MAKEFILE) ts_nxtpadsol

ts_hesscrit:
	make -f $(MAKEFILE) ts_hesscrit

ts_hesspred:
	make -f $(MAKEFILE) ts_hesspred

ts_hesspcnv:
	make -f $(MAKEFILE) ts_hesspcnv

ts_predstep:
	make -f $(MAKEFILE) ts_predstep

ts_sernewcnv:
	make -f $(MAKEFILE) ts_sernewcnv

ts_ratapp:
	make -f $(MAKEFILE) ts_ratapp

ts_corcnv:
	make -f $(MAKEFILE) ts_corcnv

ts_pcscnv:
	make -f $(MAKEFILE) ts_pcscnv

ts_scalecnv:
	make -f $(MAKEFILE) ts_scalecnv

ts_wrappade:
	make -f $(MAKEFILE) ts_wrappade

ts_fabryhom:
	make -f $(MAKEFILE) ts_fabryhom

# Test programs for the Deformations/End_Games library :

test_end_games: ts_vlprs postpoco valipoco

ts_vlprs:
	make -f $(MAKEFILE) ts_vlprs

postpoco:
	make -f $(MAKEFILE) postpoco

valipoco:
	make -f $(MAKEFILE) valipoco

# Test programs for the Deformations/Sweep library :

test_sweep: ts_parcon ts_realcont

ts_parcon:
	make -f $(MAKEFILE) ts_parcon

ts_realcont:
	make -f $(MAKEFILE) ts_realcont

# Test programs for the Deformations/Trackers library :

test_trackers: ts_preco ts_correct ts_ortocor ts_track ts_path ts_nxtsol

ts_preco:
	make -f $(MAKEFILE) ts_preco

ts_correct:
	make -f $(MAKEFILE) ts_correct

ts_ortocor:
	make -f $(MAKEFILE) ts_ortocor

ts_path:
	make -f $(MAKEFILE) ts_path

ts_track:
	make -f $(MAKEFILE) ts_track

ts_nxtsol:
	make -f $(MAKEFILE) ts_nxtsol

# Test programs for the Deformations/Continuation library :

test_continuation: ts_pactun ts_poco ts_bbpoco ts_exp ts_autotune

ts_pactun:
	make -f $(MAKEFILE) ts_pactun

ts_poco:
	make -f $(MAKEFILE) ts_poco

ts_bbpoco:
	make -f $(MAKEFILE) ts_bbpoco

ts_exp:
	make -f $(MAKEFILE) ts_exp

ts_autotune:
	make -f $(MAKEFILE) ts_autotune

mypoco:
	make -f $(MAKEFILE) mypoco

# Test programs for the Root_Counts library :

test_root_counts: test_product test_binomials test_implift test_mixedvol \
                  test_static test_dynlift test_symmetry # test_puiseux

# Test programs for the Root_Counts/Product homotopies library :

test_product: ts_startsys ts_permanent ts_prodsys ts_mainmhom ts_mainmuho \
              ts_mainsets ts_strset ts_strpart ts_permstar ts_pivsel

ts_startsys:
	make -f $(MAKEFILE) ts_startsys

ts_permanent:
	make -f $(MAKEFILE) ts_permanent

ts_prodsys:
	make -f $(MAKEFILE) ts_prodsys

ts_strpart:
	make -f $(MAKEFILE) ts_strpart

ts_mainmhom:
	make -f $(MAKEFILE) ts_mainmhom

ts_mainmuho:
	make -f $(MAKEFILE) ts_mainmuho

ts_mainsets:
	make -f $(MAKEFILE) ts_mainsets

ts_strset:
	make -f $(MAKEFILE) ts_strset

ts_permstar:
	make -f $(MAKEFILE) ts_permstar

ts_pivsel:
	make -f $(MAKEFILE) ts_pivsel

# Test programs in the Root_Counts/Binomials library :

test_binomials: ts_binsys ts_kernel ts_binset ts_incols ts_monmap ts_binsol \
                ts_affpol ts_simsys ts_simposol ts_powtrans ts_binpser

ts_binsys:
	make -f $(MAKEFILE) ts_binsys

ts_kernel:
	make -f $(MAKEFILE) ts_kernel

ts_binset:
	make -f $(MAKEFILE) ts_binset

ts_incols:
	make -f $(MAKEFILE) ts_incols

ts_monmap:
	make -f $(MAKEFILE) ts_monmap

ts_binsol:
	make -f $(MAKEFILE) ts_binsol

ts_affpol:
	make -f $(MAKEFILE) ts_affpol

ts_simsys:
	make -f $(MAKEFILE) ts_simsys

ts_simposol:
	make -f $(MAKEFILE) ts_simposol

ts_powtrans:
	make -f $(MAKEFILE) ts_powtrans

ts_binpser:
	make -f $(MAKEFILE) ts_binpser

ts_cspsol:
	make -f $(MAKEFILE) ts_cspsol

ts_unimod:
	make -f $(MAKEFILE) ts_unimod

# Test programs for the Root_Counts/Implift library :

test_implift: ts_supports ts_supspan ts_transfo ts_mainvpts ts_impvol \
              ts_drivimpl

ts_supports:
	make -f $(MAKEFILE) ts_supports

ts_supspan:
	make -f $(MAKEFILE) ts_supspan

ts_transfo:
	make -f $(MAKEFILE) ts_transfo

ts_mainvpts:
	make -f $(MAKEFILE) ts_mainvpts

ts_impvol:
	make -f $(MAKEFILE) ts_impvol

ts_drivimpl:
	make -f $(MAKEFILE) ts_drivimpl

# Test programs for the Root_Counts/Stalift library :

test_static: ts_rndcff ts_supsub ts_conint ts_mixture ts_mixvol ts_mixlab \
             ts_jpmix ts_drivcrit ts_drivstal ts_mvasta ts_exptrack \
             ts_stable ts_balance ts_userlift ts_demics

ts_rndcff:
	make -f $(MAKEFILE) ts_rndcff

ts_supsub:
	make -f $(MAKEFILE) ts_supsub

ts_conint:
	make -f $(MAKEFILE) ts_conint

ts_mixture:
	make -f $(MAKEFILE) ts_mixture

ts_mixvol:
	make -f $(MAKEFILE) ts_mixvol

ts_mixlab:
	make -f $(MAKEFILE) ts_mixlab

ts_jpmix:
	make -f $(MAKEFILE) ts_jpmix

ts_drivcrit:
	make -f $(MAKEFILE) ts_drivcrit

ts_drivstal:
	make -f $(MAKEFILE) ts_drivstal

ts_mvasta:
	make -f $(MAKEFILE) ts_mvasta

ts_exptrack:
	make -f $(MAKEFILE) ts_exptrack

ts_stable:
	make -f $(MAKEFILE) ts_stable

ts_balance:
	make -f $(MAKEFILE) ts_balance

ts_userlift:
	make -f $(MAKEFILE) ts_userlift

ts_demics:
	make -f $(MAKEFILE) ts_demics

# Test programs for the Root_Counts/Dynlift library :

test_dynlift: ts_dyntri ts_drivmink ts_drivdynl

ts_dyntri:
	make -f $(MAKEFILE) ts_dyntri

ts_drivmink:
	make -f $(MAKEFILE) ts_drivmink

ts_drivdynl:
	make -f $(MAKEFILE) ts_drivdynl

# Test programs for the Root_Counts/Symmetry library :

test_symmetry: ts_equpol ts_drivsss

ts_group:
	make -f $(MAKEFILE) ts_group

ts_equpol:
	make -f $(MAKEFILE) ts_equpol

ts_drivsss:
	make -f $(MAKEFILE) ts_drivsss

# Test programs for the Root_Counts/MixedVol library :

test_mixedvol: ts_convsub ts_mv2c ts_mv ts_mva

ts_convsub:
	make -f $(MAKEFILE) ts_convsub

ts_mv2c:
	make -f $(MAKEFILE) ts_mv2c

ts_mv:
	make -f $(MAKEFILE) ts_mv

ts_mva:
	make -f $(MAKEFILE) ts_mva

# Test programs for the Root_Counts/DEMiCs library :

test_demics: ts_calldemics demics ts_inputData ts_demicsrun ts_outputData

ts_calldemics:
	make -f $(MAKEFILE) ts_calldemics

demics:
	make -f $(MAKEFILE) demics

ts_inputData:
	make -f $(MAKEFILE) ts_inputData

ts_demicsrun:
	make -f $(MAKEFILE) ts_demicsrun

ts_outdata:
	make -f $(MAKEFILE) ts_outdata

ts_outputData:
	make -f $(MAKEFILE) ts_outputData

# Test programs for the Puiseux library :

test_puiseux: ts_degcurve ts_puiseux

ts_degcurve:
	make -f $(MAKEFILE) ts_degcurve

ts_puiseux:
	make -f $(MAKEFILE) ts_puiseux

# Test programs for the Schubert library :

test_schubert: test_sagbi test_pieri test_induction

# Test programs for SAGBI library in numerical Schubert calculus :

test_sagbi: ts_subsets ts_brackets ts_brackmons ts_brackpols ts_straighten \
            ts_expand ts_local ts_mathom ts_cheby ts_shapiro ts_eremenko \
            ts_detrock ts_sagbi

ts_subsets:
	make -f $(MAKEFILE) ts_subsets

ts_brackets:
	make -f $(MAKEFILE) ts_brackets

ts_brackmons:
	make -f $(MAKEFILE) ts_brackmons

ts_brackpols:
	make -f $(MAKEFILE) ts_brackpols

ts_straighten:
	make -f $(MAKEFILE) ts_straighten

ts_expand:
	make -f $(MAKEFILE) ts_expand

ts_local:
	make -f $(MAKEFILE) ts_local

ts_mathom:
	make -f $(MAKEFILE) ts_mathom

ts_cheby:
	make -f $(MAKEFILE) ts_cheby

ts_shapiro:
	make -f $(MAKEFILE) ts_shapiro

ts_eremenko:
	make -f $(MAKEFILE) ts_eremenko

ts_detrock:
	make -f $(MAKEFILE) ts_detrock

ts_sagbi:
	make -f $(MAKEFILE) ts_sagbi

ts_topos:
	make -f $(MAKEFILE) ts_topos

# Test programs for Pieri library in numerical Schubert calculus :

test_pieri: ts_canocurv ts_matinv ts_pieri ts_posets \
            ts_detsys ts_defpos ts_piroco # ts_org_pieri

ts_org_pieri:
	make -f $(MAKEFILE) ts_org_pieri

ts_canocurv:
	make -f $(MAKEFILE) ts_canocurv

ts_matinv:
	make -f $(MAKEFILE) ts_matinv

ts_pieri:
	make -f $(MAKEFILE) ts_pieri

ts_posets:
	make -f $(MAKEFILE) ts_posets

ts_detsys:
	make -f $(MAKEFILE) ts_detsys

ts_defpos:
	make -f $(MAKEFILE) ts_defpos

ts_piroco:
	make -f $(MAKEFILE) ts_piroco

# Test programs for Induction library in numerical Schubert calculus :

test_induction: ts_checkers ts_induce ts_flagcond ts_lrhom

ts_checkers:
	make -f $(MAKEFILE) ts_checkers

ts_induce:
	make -f $(MAKEFILE) ts_induce

ts_flagcond:
	make -f $(MAKEFILE) ts_flagcond

ts_lrhom:
	make -f $(MAKEFILE) ts_lrhom

# Test programs for the evaluators :

ts_cyclic7:
	make -f $(MAKEFILE) ts_cyclic7

solve_cyclic7:
	make -f $(MAKEFILE) solve_cyclic7

# Test programs for the Components library :

test_components: test_samplers test_interpolators test_factorization \
                 test_decomposition test_solver test_tropical

# Test programs for the Components/Samplers library :

test_samplers: ts_roots ts_hypsam ts_hyprts ts_hyppts ts_unisam \
               ts_planes ts_flag ts_coords \
               ts_itrack ts_ispan ts_inewt ts_iddnewt ts_iwset ts_ilocal \
               ts_sample ts_endgm ts_locdim ts_ddqdsam

ts_roots:
	make -f $(MAKEFILE) ts_roots

ts_hypsam:
	make -f $(MAKEFILE) ts_hypsam

ts_hyprts:
	make -f $(MAKEFILE) ts_hyprts

ts_hyppts:
	make -f $(MAKEFILE) ts_hyppts

ts_unisam:
	make -f $(MAKEFILE) ts_unisam

ts_planes:
	make -f $(MAKEFILE) ts_planes

ts_flag:
	make -f $(MAKEFILE) ts_flag

ts_coords:
	make -f $(MAKEFILE) ts_coords

ts_itrack:
	make -f $(MAKEFILE) ts_itrack

ts_ispan:
	make -f $(MAKEFILE) ts_ispan

ts_inewt:
	make -f $(MAKEFILE) ts_inewt

ts_iddnewt:
	make -f $(MAKEFILE) ts_iddnewt

ts_iwset:
	make -f $(MAKEFILE) ts_iwset

ts_ilocal:
	make -f $(MAKEFILE) ts_ilocal

ts_sampar:
	make -f $(MAKEFILE) ts_sampar

ts_sample:
	make -f $(MAKEFILE) ts_sample

ts_mulsam:
	make -f $(MAKEFILE) ts_mulsam

ts_endgm:
	make -f $(MAKEFILE) ts_endgm

ts_locdim:
	make -f $(MAKEFILE) ts_locdim

ts_ddqdsam:
	make -f $(MAKEFILE) ts_ddqdsam

# Test programs for the Components/Interpolators libary :

test_interpolators: ts_stpolin ts_mupolin ts_adpolin ts_project ts_span \
                    ts_filter ts_nbterms ts_divdif ts_newint \
                    ts_powtrc ts_traces ts_trapol ts_rewrite ts_elim \
                    ts_embed ts_square ts_squemb

ts_stpolin:
	make -f $(MAKEFILE) ts_stpolin

ts_mupolin:
	make -f $(MAKEFILE) ts_mupolin

ts_adpolin:
	make -f $(MAKEFILE) ts_adpolin

ts_project:
	make -f $(MAKEFILE) ts_project

ts_span:
	make -f $(MAKEFILE) ts_span

ts_filter:
	make -f $(MAKEFILE) ts_filter

ts_nbterms:
	make -f $(MAKEFILE) ts_nbterms

ts_divdif:
	make -f $(MAKEFILE) ts_divdif

ts_dvdexp:
	make -f $(MAKEFILE) ts_dvdexp

ts_newint:
	make -f $(MAKEFILE) ts_newint

ts_powtrc:
	make -f $(MAKEFILE) ts_powtrc

ts_traces:
	make -f $(MAKEFILE) ts_traces

ts_trapol:
	make -f $(MAKEFILE) ts_trapol

ts_rewrite:
	make -f $(MAKEFILE) ts_rewrite

ts_elim:
	make -f $(MAKEFILE) ts_elim

ts_embed:
	make -f $(MAKEFILE) ts_embed

ts_square:
	make -f $(MAKEFILE) ts_square

ts_squemb:
	make -f $(MAKEFILE) ts_squemb

# Test programs in the Components/Factorization library :

test_factorization: ts_irco ts_combfac ts_newfac ts_factor \
                    ts_hypfac ts_monfac ts_serwit # ts_rwpoly ts_ifactor

ts_rwpoly:
	make -f $(MAKEFILE) ts_rwpoly

ts_irco:
	make -f $(MAKEFILE) ts_irco

ts_combfac:
	make -f $(MAKEFILE) ts_combfac

ts_newfac:
	make -f $(MAKEFILE) ts_newfac

ts_factor:
	make -f $(MAKEFILE) ts_factor

ts_ifactor:
	make -f $(MAKEFILE) ts_ifactor

ts_hypfac:
	make -f $(MAKEFILE) ts_hypfac

ts_monfac:
	make -f $(MAKEFILE) ts_monfac

ts_serwit:
	make -f $(MAKEFILE) ts_serwit

# Test programs in the Components/Decomposition library :

test_decomposition: ts_breakup ts_cascade ts_filfac ts_irdeco ts_mbthom \
                    ts_jmpdia ts_intdia ts_witsols \
                    # ts_intcas ts_diahom ts_reorder ts_squeeze ts_dimred

ts_irdeco:
	make -f $(MAKEFILE) ts_irdeco

ts_breakup:
	make -f $(MAKEFILE) ts_breakup

ts_cascade:
	make -f $(MAKEFILE) ts_cascade

ts_intcas:
	make -f $(MAKEFILE) ts_intcas

ts_intdia:
	make -f $(MAKEFILE) ts_intdia

ts_reorder:
	make -f $(MAKEFILE) ts_reorder

ts_filfac:
	make -f $(MAKEFILE) ts_filfac

ts_squeeze:
	make -f $(MAKEFILE) ts_squeeze

ts_dimred:
	make -f $(MAKEFILE) ts_dimred

ts_mbthom:
	make -f $(MAKEFILE) ts_mbthom

ts_diahom:
	make -f $(MAKEFILE) ts_diahom

ts_fillit:
	make -f $(MAKEFILE) ts_fillit

ts_jmpdia:
	make -f $(MAKEFILE) ts_jmpdia

ts_witsols:
	make -f $(MAKEFILE) ts_witsols

# Test programs in the Components/Solver library :

test_solver: ts_hypwit ts_solver ts_rdisol

ts_hypwit:
	make -f $(MAKEFILE) ts_hypwit

ts_solver:
	make -f $(MAKEFILE) ts_solver

ts_rdisol:
	make -f $(MAKEFILE) ts_rdisol

testnine:
	make -f $(MAKEFILE) testnine

# Test programs in the Components/Tropical library

test_tropical: ts_inform ts_comfac ts_tropisms ts_pretrop \
               ts_intorto ts_tropawit ts_convhull ts_ptlstr

ts_tropawit:
	make -f $(MAKEFILE) ts_tropawit

ts_tropisms:
	make -f $(MAKEFILE) ts_tropisms

ts_inform:
	make -f $(MAKEFILE) ts_inform

ts_comfac:
	make -f $(MAKEFILE) ts_comfac

ts_intorto:
	make -f $(MAKEFILE) ts_intorto

ts_convhull:
	make -f $(MAKEFILE) ts_convhull

ts_ptlstr:
	make -f $(MAKEFILE) ts_ptlstr

ts_pretrop:
	make -f $(MAKEFILE) ts_pretrop

# test programs for the CtoPHC directory :

test_ctophc: test_funky test_state

# test programs for the CtoPHC/Funky directory :

test_funky: ts_arrays pass_poly ts_cpm ts_pieri_solver
# ts_cosup ts_cosupoly ts_cosupsys ts_cosuptoc \
# ts_getsys ts_phc_sys_rw ts_phc_sol_rw ts_pieri_count ts_csols

ts_arrays:
	make -f $(MAKEFILE) ts_arrays

pass_poly:
	make -f $(MAKEFILE) pass_poly

ts_cosup:
	make -f $(MAKEFILE) ts_cosup

ts_cosupoly:
	make -f $(MAKEFILE) ts_cosupoly

ts_cosupsys:
	make -f $(MAKEFILE) ts_cosupsys

ts_cosuptoc:
	make -f $(MAKEFILE) ts_cosuptoc

ts_getsys:
	make -f $(MAKEFILE) ts_getsys

ts_phc_sys_rw:
	make -f $(MAKEFILE) ts_phc_sys_rw

ts_phc_sol_rw:
	make -f $(MAKEFILE) ts_phc_sol_rw

ts_cpm:
	make -f $(MAKEFILE) ts_cpm

ts_pieri_count:
	make -f $(MAKEFILE) ts_pieri_count

ts_pieri_solver:
	make -f $(MAKEFILE) ts_pieri_solver

ts_phc_solver:
	make -f $(MAKEFILE) ts_phc_solver

ts_csols:
	make -f $(MAKEFILE) ts_csols

# test programs in the CtoPHC/State directory :

test_state: ts_filman ts_solcon ts_use_solcon ts_solpool ts_use_solpool \
            ts_syscon ts_use_syscon ts_syspool ts_use_syspool \
            ts_celcon ts_use_celcon ts_pirhom ts_use_c2pieri ts_c2phc \
            ts_c2str ts_use_c2phc ts_use_c2fac ts_unisolve ts_use_giftwrap \
            ts_giftwcon ts_use_sweep ts_runtrack

ts_filman:
	make -f $(MAKEFILE) ts_filman

ts_solcon:
	make -f $(MAKEFILE) ts_solcon

ts_solpool:
	make -f $(MAKEFILE) ts_solpool

ts_use_solcon:
	make -f $(MAKEFILE) ts_use_solcon

ts_use_solpool:
	make -f $(MAKEFILE) ts_use_solpool

ts_syscon:
	make -f $(MAKEFILE) ts_syscon

ts_syspool:
	make -f $(MAKEFILE) ts_syspool

ts_use_syscon:
	make -f $(MAKEFILE) ts_use_syscon

ts_use_syspool:
	make -f $(MAKEFILE) ts_use_syspool

ts_c2phc:
	make -f $(MAKEFILE) ts_c2phc

ts_c2str:
	make -f $(MAKEFILE) ts_c2str

ts_use_c2phc:
	make -f $(MAKEFILE) ts_use_c2phc

ts_pirhom:
	make -f $(MAKEFILE) ts_pirhom

ts_use_c2pieri:
	make -f $(MAKEFILE) ts_use_c2pieri

ts_use_c2fac:
	make -f $(MAKEFILE) ts_use_c2fac

ts_celcon:
	make -f $(MAKEFILE) ts_celcon

ts_use_celcon:
	make -f $(MAKEFILE) ts_use_celcon

ts_unisolve:
	make -f $(MAKEFILE) ts_unisolve

ts_use_giftwrap:
	make -f $(MAKEFILE) ts_use_giftwrap

ts_giftwcon:
	make -f $(MAKEFILE) ts_giftwcon

ts_use_sweep:
	make -f $(MAKEFILE) ts_use_sweep

ts_runtrack:
	make -f $(MAKEFILE) ts_runtrack

# files in feedback :

ts_feedback:
	make -f $(MAKEFILE) ts_feedback

# Test programs for the C interface to GPU acceleration :

test_phctoc: ts_cpu2norm_d ts_cpu2norm_dd ts_cpu2norm_qd ts_cpu2norm \
             ts_gpu2norm_d ts_gpu2norm_dd ts_gpu2norm_qd ts_gpu2norm \
             ts_gpunewton ts_gpunewton_dd ts_gpunewton_qd \
             ts_gpuonepath_d ts_gpuonepath_dd ts_gpuonepath_qd \
             ts_gpumanypaths_d ts_gpumanypaths_dd ts_gpumanypaths_qd \
             ts_gpupath_d ts_gpupath_dd ts_gpupath_qd \
             ts_adepath_d ts_adepath_dd ts_adepath_qd ts_adepath ts_pathpars

ts_cpu2norm_d:
	make -f $(MAKEFILE) ts_cpu2norm_d

ts_cpu2norm_dd:
	make -f $(MAKEFILE) ts_cpu2norm_dd

ts_cpu2norm_qd:
	make -f $(MAKEFILE) ts_cpu2norm_qd

ts_cpu2norm:
	make -f $(MAKEFILE) ts_cpu2norm

ts_gpu2norm_d:
	make -f $(MAKEFILE) ts_gpu2norm_d

ts_gpu2norm_dd:
	make -f $(MAKEFILE) ts_gpu2norm_dd

ts_gpu2norm_qd:
	make -f $(MAKEFILE) ts_gpu2norm_qd

ts_gpu2norm:
	make -f $(MAKEFILE) ts_gpu2norm

ts_gpunewton:
	make -f $(MAKEFILE) ts_gpunewton

ts_gpunewton_dd:
	make -f $(MAKEFILE) ts_gpunewton_dd

ts_gpunewton_qd:
	make -f $(MAKEFILE) ts_gpunewton_qd

ts_gpuonepath_d:
	make -f $(MAKEFILE) ts_gpuonepath_d

ts_gpuonepath_dd:
	make -f $(MAKEFILE) ts_gpuonepath_dd

ts_gpuonepath_qd:
	make -f $(MAKEFILE) ts_gpuonepath_qd

ts_gpumanypaths_d:
	make -f $(MAKEFILE) ts_gpumanypaths_d

ts_gpumanypaths_dd:
	make -f $(MAKEFILE) ts_gpumanypaths_dd

ts_gpumanypaths_qd:
	make -f $(MAKEFILE) ts_gpumanypaths_qd

ts_gpupath_d:
	make -f $(MAKEFILE) ts_gpupath_d

ts_gpupath_dd:
	make -f $(MAKEFILE) ts_gpupath_dd

ts_gpupath_qd:
	make -f $(MAKEFILE) ts_gpupath_qd

ts_adepath_d:
	make -f $(MAKEFILE) ts_adepath_d

ts_adepath_dd:
	make -f $(MAKEFILE) ts_adepath_dd

ts_adepath_qd:
	make -f $(MAKEFILE) ts_adepath_qd

ts_adepath:
	make -f $(MAKEFILE) ts_adepath

ademanypaths.o:
	make -f $(MAKEFILE) ademanypaths.o

lib2path.o:
	make -f $(MAKEFILE) lib2path.o

lib2path_d.o:
	make -f $(MAKEFILE) lib2path_d.o

ts_ademanypaths:
	make -f $(MAKEFILE) ts_ademanypaths

ts_pathpars:
	make -f $(MAKEFILE) ts_pathpars

# Test programs for the Tasking:

test_tasking: ts_mutex ts_jobs ts_mtvol ts_mtvolcon ts_mtverify ts_mtcont \
              ts_mtmva ts_mtddref ts_mtqdref ts_mtsols ts_mtsharp ts_mtperm \
              ts_mtpermtrol ts_mtmbthom ts_mtcelidx ts_corecount ts_mtadcnv \
	      ts_mtserlin ts_mtnewton ts_mthessian ts_mtratapp ts_mtshiftcnv \
	      ts_mtpcscnv ts_mthesscrc ts_mtprdcnv
#             ts_mtmatvec ts_mtpolval ts_mtlinsol

ts_corecount:
	make -f $(MAKEFILE) ts_corecount

ts_mtmatvec:
	make -f $(MAKEFILE) ts_mtmatvec

ts_mtpolval:
	make -f $(MAKEFILE) ts_mtpolval

ts_mtlinsol:
	make -f $(MAKEFILE) ts_mtlinsol

ts_mtlinsolV0:
	make -f $(MAKEFILE) ts_mtlinsolV0

ts_mtlinsolv1:
	make -f $(MAKEFILE) ts_mtlinsolv1

ts_tasking:
	make -f $(MAKEFILE) ts_tasking

ts_mutex:
	make -f $(MAKEFILE) ts_mutex

ts_jobs:
	make -f $(MAKEFILE) ts_jobs

ts_mtsols:
	make -f $(MAKEFILE) ts_mtsols

ts_mtstref:
	make -f $(MAKEFILE) ts_mtstref

ts_mtddref:
	make -f $(MAKEFILE) ts_mtddref

ts_mtqdref:
	make -f $(MAKEFILE) ts_mtqdref

ts_mtcont:
	make -f $(MAKEFILE) ts_mtcont

ts_mtsharp:
	make -f $(MAKEFILE) ts_mtsharp

ts_mtpolco:
	make -f $(MAKEFILE) ts_mtpolco

ts_mtsam:
	make -f $(MAKEFILE) ts_mtsam

ts_mtvol:
	make -f $(MAKEFILE) ts_mtvol

ts_mtvolcon:
	make -f $(MAKEFILE) ts_mtvolcon

ts_mtverify:
	make -f $(MAKEFILE) ts_mtverify

ts_mtmva:
	make -f $(MAKEFILE) ts_mtmva

ts_mtperm:
	make -f $(MAKEFILE) ts_mtperm

ts_mtpermtrol:
	make -f $(MAKEFILE) ts_mtpermtrol

ts_mtmbthom:
	make -f $(MAKEFILE) ts_mtmbthom

ts_mtcelidx:
	make -f $(MAKEFILE) ts_mtcelidx

ts_mtadcnv:
	make -f $(MAKEFILE) ts_mtadcnv

ts_mtserlin:
	make -f $(MAKEFILE) ts_mtserlin

ts_mtnewton:
	make -f $(MAKEFILE) ts_mtnewton

ts_mthessian:
	make -f $(MAKEFILE) ts_mthessian

ts_mthesscrc:
	make -f $(MAKEFILE) ts_mthesscrc

ts_mtratapp:
	make -f $(MAKEFILE) ts_mtratapp

ts_mtshiftcnv:
	make -f $(MAKEFILE) ts_mtshiftcnv

ts_mtpcscnv:
	make -f $(MAKEFILE) ts_mtpcscnv

ts_mtprdcnv:
	make -f $(MAKEFILE) ts_mtprdcnv

# Test programs for the Main interactive drivers :

test_main: ts_mainroco ts_blkroco ts_greet use_phc ts_actopt ts_opthand \
           ts_bbsolve

ts_mainroco:
	make -f $(MAKEFILE) ts_mainroco

ts_blkroco:
	make -f $(MAKEFILE) ts_blkroco

use_phc:
	make -f $(MAKEFILE) use_phc

ts_greet:
	make -f $(MAKEFILE) ts_greet

ts_actopt:
	make -f $(MAKEFILE) ts_actopt

ts_opthand:
	make -f $(MAKEFILE) ts_opthand

ts_bbsolve:
	make -f $(MAKEFILE) ts_bbsolve

# Test programs in the C folder :

test_lib: lib_syscon lib_syscon_cpp lib_solcon lib_solcon_cpp \
          lib_solrep lib_solrep_cpp lib_tabform lib_tabform_cpp \
          lib_read lib_read_cpp lib_newton lib_newton_cpp \
          lib_roco lib_roco_cpp lib_celcon lib_celcon_cpp \
          lib_celphc lib_celphc_cpp lib_scale lib_scale_cpp \
          lib_reduce lib_reduce_cpp lib_multip lib_multip_cpp \
          lib_conpar lib_conpar_cpp lib_poco lib_poco_cpp \
          lib_track lib_track_cpp lib_nxtsol lib_nxtsol_cpp \
          lib_sweep lib_sweep_cpp lib_nbrvar lib_nbrvar_cpp \
          lib_factor lib_factor_cpp lib_member lib_member_cpp \
          lib_pieri lib_pieri_cpp lib_lrhom lib_lrhom_cpp \
          lib_numbtrop lib_numbtrop_cpp lib_series lib_series_cpp \
          lib_padcon lib_padcon_cpp \
          lib_usolve lib_usolve_cpp phc_solve phc_solve_cpp \
          lib_giftwrap lib_giftwrap_cpp ts_list2str \
          lib_adepath_d lib_adepath_dd lib_adepath_qd \
          lib_witset lib_witset_cpp lib_witsols lib_witsols_cpp

ts_list2str:
	make -f $(MAKEFILE) ts_list2str

lib_syscon:
	make -f $(MAKEFILE) lib_syscon

lib_syscon_cpp:
	make -f $(MAKEFILE) lib_syscon_cpp

lib_tabform:
	make -f $(MAKEFILE) lib_tabform

lib_tabform_cpp:
	make -f $(MAKEFILE) lib_tabform_cpp

lib_solcon:
	make -f $(MAKEFILE) lib_solcon

lib_solcon_cpp:
	make -f $(MAKEFILE) lib_solcon_cpp

lib_solrep:
	make -f $(MAKEFILE) lib_solrep

lib_solrep_cpp:
	make -f $(MAKEFILE) lib_solrep_cpp

lib_read:
	make -f $(MAKEFILE) lib_read

lib_read_cpp:
	make -f $(MAKEFILE) lib_read_cpp

lib_newton:
	make -f $(MAKEFILE) lib_newton

lib_newton_cpp:
	make -f $(MAKEFILE) lib_newton_cpp

lib_multip:
	make -f $(MAKEFILE) lib_multip

lib_multip_cpp:
	make -f $(MAKEFILE) lib_multip_cpp

lib_roco:
	make -f $(MAKEFILE) lib_roco

lib_roco_cpp:
	make -f $(MAKEFILE) lib_roco_cpp

lib_celcon:
	make -f $(MAKEFILE) lib_celcon

lib_celcon_cpp:
	make -f $(MAKEFILE) lib_celcon_cpp

lib_celphc:
	make -f $(MAKEFILE) lib_celphc

lib_celphc_cpp:
	make -f $(MAKEFILE) lib_celphc_cpp

lib_scale:
	make -f $(MAKEFILE) lib_scale

lib_scale_cpp:
	make -f $(MAKEFILE) lib_scale_cpp

lib_reduce:
	make -f $(MAKEFILE) lib_reduce

lib_reduce_cpp:
	make -f $(MAKEFILE) lib_reduce_cpp

lib_conpar:
	make -f $(MAKEFILE) lib_conpar

lib_conpar_cpp:
	make -f $(MAKEFILE) lib_conpar_cpp

lib_poco:
	make -f $(MAKEFILE) lib_poco

lib_poco_cpp:
	make -f $(MAKEFILE) lib_poco_cpp

lib_track:
	make -f $(MAKEFILE) lib_track

lib_track_cpp:
	make -f $(MAKEFILE) lib_track_cpp

lib_nxtsol:
	make -f $(MAKEFILE) lib_nxtsol

lib_nxtsol_cpp:
	make -f $(MAKEFILE) lib_nxtsol_cpp

lib_sweep:
	make -f $(MAKEFILE) lib_sweep

lib_sweep_cpp:
	make -f $(MAKEFILE) lib_sweep_cpp

lib_factor:
	make -f $(MAKEFILE) lib_factor

lib_factor_cpp:
	make -f $(MAKEFILE) lib_factor_cpp

lib_member:
	make -f $(MAKEFILE) lib_member

lib_member_cpp:
	make -f $(MAKEFILE) lib_member_cpp

lib_pieri:
	make -f $(MAKEFILE) lib_pieri

lib_pieri_cpp:
	make -f $(MAKEFILE) lib_pieri_cpp

lib_lrhom:
	make -f $(MAKEFILE) lib_lrhom

lib_lrhom_cpp:
	make -f $(MAKEFILE) lib_lrhom_cpp

lib_usolve:
	make -f $(MAKEFILE) lib_usolve

lib_usolve_cpp:
	make -f $(MAKEFILE) lib_usolve_cpp

lib_giftwrap:
	make -f $(MAKEFILE) lib_giftwrap

lib_giftwrap_cpp:
	make -f $(MAKEFILE) lib_giftwrap_cpp

lib_numbtrop:
	make -f $(MAKEFILE) lib_numbtrop

lib_numbtrop_cpp:
	make -f $(MAKEFILE) lib_numbtrop_cpp

lib_series:
	make -f $(MAKEFILE) lib_series

lib_series_cpp:
	make -f $(MAKEFILE) lib_series_cpp

lib_padcon:
	make -f $(MAKEFILE) lib_padcon

lib_padcon_cpp:
	make -f $(MAKEFILE) lib_padcon_cpp

lib_nbrvar:
	make -f $(MAKEFILE) lib_nbrvar

lib_nbrvar_cpp:
	make -f $(MAKEFILE) lib_nbrvar_cpp

lib_adepath_d:
	make -f $(MAKEFILE) lib_adepath_d

lib_adepath_dd:
	make -f $(MAKEFILE) lib_adepath_dd

lib_adepath_qd:
	make -f $(MAKEFILE) lib_adepath_qd

phc_solve:
	make -f $(MAKEFILE) phc_solve

phc_solve_cpp:
	make -f $(MAKEFILE) phc_solve_cpp

lib_witset:
	make -f $(MAKEFILE) lib_witset

lib_witset_cpp:
	make -f $(MAKEFILE) lib_witset_cpp

lib_witsols:
	make -f $(MAKEFILE) lib_witsols

lib_witsols_cpp:
	make -f $(MAKEFILE) lib_witsols_cpp

# Executables with MPI:

test_mpi: adatest mpi2sys mpi2sol mpi2hom mpi2phc_s mpi2phc_d2 mpi2track \
          mpi2phc_dd mpi2phc_qd \
          mpi2cell_s mpi2cell_d mpi2fac mpi2fac_s mpi2fac_d \
          parallel_pieri parallel_subsystem monodromyMPI mpi2padcon

adatest:
	make -f $(MAKEFILE) adatest

mpi2sys:
	make -f $(MAKEFILE) mpi2sys

mpi2sol:
	make -f $(MAKEFILE) mpi2sol

mpi2hom:
	make -f $(MAKEFILE) mpi2hom

mpi2phc_s:
	make -f $(MAKEFILE) mpi2phc_s

mpi2phc_d2:
	make -f $(MAKEFILE) mpi2phc_d2

mpi2phc_dd:
	make -f $(MAKEFILE) mpi2phc_dd

mpi2phc_qd:
	make -f $(MAKEFILE) mpi2phc_qd

mpi2track:
	make -f $(MAKEFILE) mpi2track

parallel_pieri:
	make -f $(MAKEFILE) parallel_pieri

parallel_subsystem:
	make -f $(MAKEFILE) parallel_subsystem

mpi2fac:
	make -f $(MAKEFILE) mpi2fac

mpi2fac_s:
	make -f $(MAKEFILE) mpi2fac_s

mpi2fac_d:
	make -f $(MAKEFILE) mpi2fac_d

monodromyMPI:
	make -f $(MAKEFILE) monodromyMPI

mpi2cell_s:
	make -f $(MAKEFILE) mpi2cell_s

mpi2cell_d:
	make -f $(MAKEFILE) mpi2cell_d

mpi2padcon:
	make -f $(MAKEFILE) mpi2padcon

# cleaning up the object files :

clean:
	make -f $(MAKEFILE) clean
