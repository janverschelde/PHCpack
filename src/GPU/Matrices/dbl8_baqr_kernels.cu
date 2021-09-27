/* The file dbl8_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl8_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#endif
#include "dbl8_baqr_kernels.h"
#include "octo_double_functions.h"
#include "dbl_baqr_flopcounts.h"

using namespace std;

__global__ void dbl8_small_house
 ( double *x0hihihi, double *x0lohihi, double *x0hilohi, double *x0lolohi,
   double *x0hihilo, double *x0lohilo, double *x0hilolo, double *x0lololo,
   double *x1hihihi, double *x1lohihi, double *x1hilohi, double *x1lolohi,
   double *x1hihilo, double *x1lohilo, double *x1hilolo, double *x1lololo,
   int dim, int dimLog2,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
}

__global__ void cmplx8_small_house
 ( double *x0rehihihi, double *x0relohihi,
   double *x0rehilohi, double *x0relolohi,
   double *x0rehihilo, double *x0relohilo,
   double *x0rehilolo, double *x0relololo,
   double *x0imhihihi, double *x0imlohihi,
   double *x0imhilohi, double *x0imlolohi,
   double *x0imhihilo, double *x0imlohilo,
   double *x0imhilolo, double *x0imlololo,
   double *x1rehihihi, double *x1relohihi,
   double *x1rehilohi, double *x1relolohi,
   double *x1rehihilo, double *x1relohilo,
   double *x1rehilolo, double *x1relololo,
   double *x1imhihihi, double *x1imlohihi,
   double *x1imhilohi, double *x1imlolohi,
   double *x1imhihilo, double *x1imlohilo,
   double *x1imhilolo, double *x1imlololo,
   int dim, int dimLog2,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
}

__global__ void dbl8_large_sum_of_squares
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *sumshihihi, double *sumshilohi,
   double *sumslohihi, double *sumslolohi,
   double *sumshihilo, double *sumshilolo,
   double *sumslohilo, double *sumslololo, int dim, int BS, int BSLog2 )
{
}

__global__ void cmplx8_large_sum_of_squares
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int dim, int BS, int BSLog2 )
{
}

__global__ void dbl8_sum_accumulator
 ( double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo,
   int nbsums, int nbsumsLog2,
   double *acchihihi, double *acclohihi,
   double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo,
   double *acchilolo, double *acclololo )
{
}

__global__ void dbl8_normalize
 ( int dim, int szt,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *v0hihihi, double *v0lohihi, double *v0hilohi, double *v0lolohi,
   double *v0hihilo, double *v0lohilo, double *v0hilolo, double *v0lololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo )
{
}

__global__ void cmplx8_normalize
 ( int dim, int szt,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *inv0rehihihi, double *inv0relohihi,
   double *inv0rehilohi, double *inv0relolohi,
   double *inv0rehihilo, double *inv0relohilo,
   double *inv0rehilolo, double *inv0relololo,
   double *inv0imhihihi, double *inv0imlohihi,
   double *inv0imhilohi, double *inv0imlolohi,
   double *inv0imhihilo, double *inv0imlohilo,
   double *inv0imhilolo, double *inv0imlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo )
{
}

__global__ void dbl8_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
}

__global__ void cmplx8_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi, 
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
}

__global__ void dbl8_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *RTdotvhihihi, double *RTdotvlohihi,
   double *RTdotvhilohi, double *RTdotvlolohi,
   double *RTdotvhihilo, double *RTdotvlohilo,
   double *RTdotvhilolo, double *RTdotvlololo )
{
}

__global__ void cmplx8_RHdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo )
{
}

__global__ void dbl8_sum_betaRTdotv
 ( int nrows,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *RTdotvhihihi, double *RTdotvlohihi,
   double *RTdotvhilohi, double *RTdotvlolohi,
   double *RTdotvhihilo, double *RTdotvlohilo,
   double *RTdotvhilolo, double *RTdotvlololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo )
{
}

__global__ void cmplx8_sum_betaRHdotv
 ( int nrows,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *RTdotvrehihihi, double *RTdotvrelohihi,
   double *RTdotvrehilohi, double *RTdotvrelolohi,
   double *RTdotvrehihilo, double *RTdotvrelohilo,
   double *RTdotvrehilolo, double *RTdotvrelololo,
   double *RTdotvimhihihi, double *RTdotvimlohihi,
   double *RTdotvimhilohi, double *RTdotvimlolohi,
   double *RTdotvimhihilo, double *RTdotvimlohilo,
   double *RTdotvimhilolo, double *RTdotvimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
}

__global__ void dbl8_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo )
{
}

__global__ void cmplx8_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
}

__global__ void dbl8_beta_times_V
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo )
{
}

__global__ void cmplx8_beta_times_V
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi,
   double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo,
   double *Wimhilolo, double *Wimlololo )
{
}

__global__ void dbl8_initialize_WYT
 ( int dim, int szt,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
}

__global__ void cmplx8_initialize_WYH
 ( int dim, int szt,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo )
{
}

__global__ void dbl8_update_WYT
 ( int dim, int szt,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
}

__global__ void cmplx8_update_WYH
 ( int dim, int szt,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYHrehihihi, double *WYHrelohihi,
   double *WYHrehilohi, double *WYHrelolohi,
   double *WYHrehihilo, double *WYHrelohilo,
   double *WYHrehilolo, double *WYHrelololo,
   double *WYHimhihihi, double *WYHimlohihi,
   double *WYHimhilohi, double *WYHimlolohi,
   double *WYHimhihilo, double *WYHimlohilo,
   double *WYHimhilolo, double *WYHimlololo )
{
}

__global__ void dbl8_beta_next_W
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
}

__global__ void cmplx8_beta_next_W
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYHrehihihi, double *WYHrelohihi,
   double *WYHrehilohi, double *WYHrelolohi,
   double *WYHrehihilo, double *WYHrelohilo,
   double *WYHrehilolo, double *WYHrelololo,
   double *WYHimhihihi, double *WYHimlohihi,
   double *WYHimhilohi, double *WYHimlolohi,
   double *WYHimhihilo, double *WYHimlohilo,
   double *WYHimhilolo, double *WYHimlololo )
{
}

__global__ void dbl8_small_WYT
 ( int nrows, int szt,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
}

__global__ void cmplx8_small_WYH
 ( int nrows, int szt,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *Yrehihihi, double *Yrelohihi, double *Yrehilohi, double *Yrelolohi,
   double *Yrehihilo, double *Yrelohilo, double *Yrehilolo, double *Yrelololo,
   double *Yimhihihi, double *Yimlohihi, double *Yimhilohi, double *Yimlolohi,
   double *Yimhihilo, double *Yimlohilo, double *Yimhilolo, double *Yimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo )
{
}

__global__ void dbl8_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihihi, double *Qlohihi, double *Qhilohi, double *Qlolohi,
   double *Qhihilo, double *Qlohilo, double *Qhilolo, double *Qlololo,
   double *WYThihihi, double *WYTlohihi, double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo, double *WYThilolo, double *WYTlololo,
   double *QWYThihihi, double *QWYTlohihi,
   double *QWYThilohi, double *QWYTlolohi,
   double *QWYThihilo, double *QWYTlohilo,
   double *QWYThilolo, double *QWYTlololo )
{
}

__global__ void cmplx8_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihihi, double *Qrelohihi, double *Qrehilohi, double *Qrelolohi,
   double *Qrehihilo, double *Qrelohilo, double *Qrehilolo, double *Qrelololo,
   double *Qimhihihi, double *Qimlohihi, double *Qimhilohi, double *Qimlolohi,
   double *Qimhihilo, double *Qimlohilo, double *Qimhilolo, double *Qimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo,
   double *QWYTrehihihi, double *QWYTrelohihi,
   double *QWYTrehilohi, double *QWYTrelolohi,
   double *QWYTrehihilo, double *QWYTrelohilo,
   double *QWYTrehilolo, double *QWYTrelololo,
   double *QWYTimhihihi, double *QWYTimlohihi,
   double *QWYTimhilohi, double *QWYTimlolohi,
   double *QWYTimhihilo, double *QWYTimlohilo,
   double *QWYTimhilolo, double *QWYTimlololo )
{
}

__global__ void dbl8_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWThihihi, double *YWTlohihi, double *YWThilohi, double *YWTlolohi,
   double *YWThihilo, double *YWTlohilo, double *YWThilolo, double *YWTlololo,
   double *Chihihi, double *Clohihi, double *Chilohi, double *Clolohi,
   double *Chihilo, double *Clohilo, double *Chilolo, double *Clololo,
   double *YWTChihihi, double *YWTClohihi,
   double *YWTChilohi, double *YWTClolohi,
   double *YWTChihilo, double *YWTClohilo,
   double *YWTChilolo, double *YWTClololo )
{
}

__global__ void cmplx8_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWTrehihihi, double *YWTrelohihi,
   double *YWTrehilohi, double *YWTrelolohi,
   double *YWTrehihilo, double *YWTrelohilo,
   double *YWTrehilolo, double *YWTrelololo,
   double *YWTimhihihi, double *YWTimlohihi,
   double *YWTimhilohi, double *YWTimlolohi,
   double *YWTimhihilo, double *YWTimlohilo,
   double *YWTimhilolo, double *YWTimlololo,
   double *Crehihihi, double *Crelohihi, double *Crehilohi, double *Crelolohi,
   double *Crehihilo, double *Crelohilo, double *Crehilolo, double *Crelololo,
   double *Cimhihihi, double *Cimlohihi, double *Cimhilohi, double *Cimlolohi,
   double *Cimhihilo, double *Cimlohilo, double *Cimhilolo, double *Cimlololo,
   double *YWTCrehihihi, double *YWTCrelohihi,
   double *YWTCrehilohi, double *YWTCrelolohi,
   double *YWTCrehihilo, double *YWTCrelohilo,
   double *YWTCrehilolo, double *YWTCrelololo,
   double *YWTCimhihihi, double *YWTCimlohihi,
   double *YWTCimhilohi, double *YWTCimlolohi,
   double *YWTCimhihilo, double *YWTCimlohilo,
   double *YWTCimhilolo, double *YWTCimlololo )
{
}

__global__ void dbl8_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihihi, double *Qlohihi, double *Qhilohi, double *Qlolohi,
   double *Qhihilo, double *Qlohilo, double *Qhilolo, double *Qlololo,
   double *QWYThihihi, double *QWYTlohihi,
   double *QWYThilohi, double *QWYTlolohi,
   double *QWYThihilo, double *QWYTlohilo,
   double *QWYThilolo, double *QWYTlololo )
{
}

__global__ void cmplx8_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihihi, double *Qrelohihi, double *Qrehilohi, double *Qrelolohi,
   double *Qrehihilo, double *Qrelohilo, double *Qrehilolo, double *Qrelololo,
   double *Qimhihihi, double *Qimlohihi, double *Qimhilohi, double *Qimlolohi,
   double *Qimhihilo, double *Qimlohilo, double *Qimhilolo, double *Qimlololo,
   double *QWYTrehihihi, double *QWYTrelohihi,
   double *QWYTrehilohi, double *QWYTrelolohi,
   double *QWYTrehihilo, double *QWYTrelohilo,
   double *QWYTrehilolo, double *QWYTrelololo,
   double *QWYTimhihihi, double *QWYTimlohihi,
   double *QWYTimhilohi, double *QWYTimlolohi,
   double *QWYTimhihilo, double *QWYTimlohilo,
   double *QWYTimhilolo, double *QWYTimlololo )
{
}

__global__ void dbl8_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *YWTChihihi, double *YWTClohihi,
   double *YWTChilohi, double *YWTClolohi,
   double *YWTChihilo, double *YWTClohilo,
   double *YWTChilolo, double *YWTClololo )
{
}

__global__ void cmplx8_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *YWTCrehihihi, double *YWTCrelohihi,
   double *YWTCrehilohi, double *YWTCrelolohi,
   double *YWTCrehihilo, double *YWTCrelohilo,
   double *YWTCrehilolo, double *YWTCrelololo,
   double *YWTCimhihihi, double *YWTCimlohihi,
   double *YWTCimhilohi, double *YWTCimlolohi,
   double *YWTCimhihilo, double *YWTCimlohilo,
   double *YWTCimhilolo, double *YWTCimlololo )
{
}

void GPU_dbl8_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   const int nrLog2 = ceil(log2((double) nrows1));
   const int rowidx = colidx*(nrows+1);       // start of number in A_h
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   if(verbose)
   {
      cout << "nrows : " << nrows
           << "  nVrows : " << nVrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  nbt : " << nbt << endl;
      cout << "k : " << k 
           << "  L : " << L
           << "  nrows1 : " << nrows1
           << "  colidx : " << colidx
           << "  rowidx : " << rowidx << endl;
   }
   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vhihihi_h[i] = 0.0;
         vlohihi_h[i] = 0.0;
         vhilohi_h[i] = 0.0;
         vlolohi_h[i] = 0.0;
         vhihilo_h[i] = 0.0;
         vlohilo_h[i] = 0.0;
         vhilolo_h[i] = 0.0;
         vlololo_h[i] = 0.0;
      }
      cudaMemcpy(&Vhihihi_d[L*nVrows],vhihihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohihi_d[L*nVrows],vlohihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilohi_d[L*nVrows],vhilohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlolohi_d[L*nVrows],vlolohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhihilo_d[L*nVrows],vhihilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohilo_d[L*nVrows],vlohilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilolo_d[L*nVrows],vhilolo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlololo_d[L*nVrows],vlololo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      betahihihi_h[L] = 0.0;
      betalohihi_h[L] = 0.0;
      betahilohi_h[L] = 0.0;
      betalolohi_h[L] = 0.0;
      betahihilo_h[L] = 0.0;
      betalohilo_h[L] = 0.0;
      betahilolo_h[L] = 0.0;
      betalololo_h[L] = 0.0;
      vhihihi_h[0] = 1.0;
      vlohihi_h[0] = 0.0;
      vhilohi_h[0] = 0.0;
      vlolohi_h[0] = 0.0;
      vhihilo_h[0] = 1.0;
      vlohilo_h[0] = 0.0;
      vhilolo_h[0] = 0.0;
      vlololo_h[0] = 0.0;
      cudaMemcpy(&betahihihi_d[L],&betahihihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohihi_d[L],&betalohihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilohi_d[L],&betahilohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalolohi_d[L],&betalolohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahihilo_d[L],&betahihilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohilo_d[L],&betalohilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilolo_d[L],&betahilolo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalololo_d[L],&betalololo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhihihi_d[L*nVrows+L],vhihihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohihi_d[L*nVrows+L],vlohihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilohi_d[L*nVrows+L],vhilohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlolohi_d[L*nVrows+L],vlolohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhihilo_d[L*nVrows+L],vhihilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohilo_d[L*nVrows+L],vlohilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilolo_d[L*nVrows+L],vhilolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlololo_d[L*nVrows+L],vlololo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl8_small_house<<<1,nrows1>>>
         (&Ahihihi_d[rowidx],&Alohihi_d[rowidx],
          &Ahilohi_d[rowidx],&Alolohi_d[rowidx],
          &Ahihilo_d[rowidx],&Alohilo_d[rowidx],
          &Ahilolo_d[rowidx],&Alololo_d[rowidx],
          &Ahihihi_d[rowidx+1],&Alohihi_d[rowidx+1],
          &Ahilohi_d[rowidx+1],&Alolohi_d[rowidx+1],
          &Ahihilo_d[rowidx+1],&Alohilo_d[rowidx+1],
          &Ahilolo_d[rowidx+1],&Alololo_d[rowidx+1],
          nrows1,nrLog2,&Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
                        &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
                        &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
                        &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
          &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
          &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_small_house(nrows1,nrLog2,add,mul,div,sqrtfun);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahihihi_h[L],&betahihihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohihi_h[L],&betalohihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilohi_h[L],&betahilohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalolohi_h[L],&betalolohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahihilo_h[L],&betahihilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohilo_h[L],&betalohilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilolo_h[L],&betahilolo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalololo_h[L],&betalololo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihihi_h,&Vhihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohihi_h,&Vlohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilohi_h,&Vhilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlolohi_h,&Vlolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihilo_h,&Vhihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohilo_h,&Vlohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilolo_h,&Vhilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlololo_h,&Vlololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);

      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "          "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "          "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "          "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;

      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : "
              << vhihihi_h[i] << "  " << vlohihi_h[i] << endl
              << "       "
              << vhilohi_h[i] << "  " << vlolohi_h[i] << endl
              << "       "
              << vhihilo_h[i] << "  " << vlohilo_h[i] << endl
              << "       "
              << vhilolo_h[i] << "  " << vlololo_h[i] << endl;
   }
}

void GPU_cmplx8_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *vrehihihi_h, double *vrelohihi_h,
   double *vrehilohi_h, double *vrelolohi_h,
   double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h,
   double *vimhilohi_h, double *vimlolohi_h,
   double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
}

void GPU_dbl8_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *sumshihihi_h, double *sumslohihi_h,
   double *sumshilohi_h, double *sumslolohi_h,
   double *sumshihilo_h, double *sumslohilo_h,
   double *sumshilolo_h, double *sumslololo_h,
   double *sumshihihi_d, double *sumslohihi_d,
   double *sumshilohi_d, double *sumslolohi_d,
   double *sumshihilo_d, double *sumslohilo_d,
   double *sumshilolo_d, double *sumslololo_d,
   double *sigmahihihi_h, double *sigmalohihi_h, 
   double *sigmahilohi_h, double *sigmalolohi_h, 
   double *sigmahihilo_h, double *sigmalohilo_h, 
   double *sigmahilolo_h, double *sigmalololo_h, 
   double *sigmahihihi_d, double *sigmalohihi_d,
   double *sigmahilohi_d, double *sigmalolohi_d,
   double *sigmahihilo_d, double *sigmalohilo_d,
   double *sigmahilolo_d, double *sigmalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   // nrows1 = nrows - colidx - 1 = size of Householder vector
   const int nblocks = ceil(((double) nrows1)/szt); // sufficient threads
   const int nblLog2 = ceil(log2((double) nblocks));
   const int sztLog2 = ceil(log2((double) szt));
   const int rowidx = colidx*(nrows+1);         // start of number in A_h
   const int nVrows = nrows - k*szt;             // dimension of V matrix

   cudaEvent_t start,stop;            // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vhihihi_h[i] = 0.0;
         vlohihi_h[i] = 0.0;
         vhilohi_h[i] = 0.0;
         vlolohi_h[i] = 0.0;
         vhihilo_h[i] = 0.0;
         vlohilo_h[i] = 0.0;
         vhilolo_h[i] = 0.0;
         vlololo_h[i] = 0.0;
      }
   }
   vhihihi_h[L] = 1.0;                    // set one on the diagonal

   cudaMemcpy(&Vhihihi_d[L*nVrows],vhihihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlohihi_d[L*nVrows],vlohihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vhilohi_d[L*nVrows],vhilohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlolohi_d[L*nVrows],vlolohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vhihilo_d[L*nVrows],vhihilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlohilo_d[L*nVrows],vlohilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vhilolo_d[L*nVrows],vhilolo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlololo_d[L*nVrows],vlololo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);

   if(verbose)
   {
      cout << "-> launching " << nblocks << " blocks of "
           << szt << " threads to compute the sum of squares ..." << endl;
      cout << "nrows1 : " << nrows1 << "  rowidx : " << rowidx;
      cout << "  ceil(log2(#blocks)) : " << nblLog2;
      cout << "  ceil(log2(szt)) : " << sztLog2 << endl;
   }
   for(int i=0; i<nblocks; i++)
   {
      sumshihihi_h[i] = 0.0;
      sumslohihi_h[i] = 0.0;
      sumshilohi_h[i] = 0.0;
      sumslolohi_h[i] = 0.0;
      sumshihilo_h[i] = 0.0;
      sumslohilo_h[i] = 0.0;
      sumshilolo_h[i] = 0.0;
      sumslololo_h[i] = 0.0;
   }
   cudaMemcpy(sumshihihi_d,sumshihihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohihi_d,sumslohihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilohi_d,sumshilohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslolohi_d,sumslolohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshihilo_d,sumshihilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohilo_d,sumslohilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilolo_d,sumshilolo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslololo_d,sumslololo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl8_large_sum_of_squares<<<nblocks,szt>>>
      (&Ahihihi_d[rowidx+1],&Alohihi_d[rowidx+1],
       &Ahilohi_d[rowidx+1],&Alolohi_d[rowidx+1],
       &Ahihilo_d[rowidx+1],&Alohilo_d[rowidx+1],
       &Ahilolo_d[rowidx+1],&Alololo_d[rowidx+1],
       sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
       sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
       nrows1,szt,sztLog2);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_large_sum_of_squares(nblocks,szt,sztLog2,add,mul);

   if(verbose)
   {
      cout << "-> launching 1 block of " << nblocks
           << " threads to accumulate the sums ..." << endl;
   }
   cudaEventRecord(start);
   dbl8_sum_accumulator<<<1,nblocks>>>
      (sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
       sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
       nblocks,nblLog2,
       sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
       sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_sum_accumulator(nblocks,nblLog2,add);

   cudaMemcpy(sigmahihihi_h,sigmahihihi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohihi_h,sigmalohihi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilohi_h,sigmahilohi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalolohi_h,sigmalolohi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahihilo_h,sigmahihilo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohilo_h,sigmalohilo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilolo_h,sigmahilolo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalololo_h,sigmalololo_d,sizeof(double),
              cudaMemcpyDeviceToHost);

   bool done = false;

   if((sigmahihihi_h[0] == 0.0) && (sigmalohihi_h[0] == 0.0) &&
      (sigmahilohi_h[0] == 0.0) && (sigmalolohi_h[0] == 0.0) &&
      (sigmahihilo_h[0] == 0.0) && (sigmalohilo_h[0] == 0.0) &&
      (sigmahilolo_h[0] == 0.0) && (sigmalololo_h[0] == 0.0))
   {
      betahihihi_h[L] = 0.0; betalohihi_h[L] = 0.0;
      betahilohi_h[L] = 0.0; betalolohi_h[L] = 0.0;
      betahihilo_h[L] = 0.0; betalohilo_h[L] = 0.0;
      betahilolo_h[L] = 0.0; betalololo_h[L] = 0.0;
      done = true;

      if(verbose)
         cout << "Zero sigma value encountered." << endl;
   }
   else // beta is computed on the host instead of by one GPU thread
   {
      // const double x0hi = Ahi_h[rowidx];
      // const double x0lo = Alo_h[rowidx];
      double acchihihi,acclohihi,acchilohi,acclolohi;
      double acchihilo,acclohilo,acchilolo,acclololo;
      double muhihihi,mulohihi,muhilohi,mulolohi;
      double muhihilo,mulohilo,muhilolo,mulololo;
      double v0hihihi,v0lohihi,v0hilohi,v0lolohi;
      double v0hihilo,v0lohilo,v0hilolo,v0lololo;
      double v0p2hihihi,v0p2lohihi,v0p2hilohi,v0p2lolohi;
      double v0p2hihilo,v0p2lohilo,v0p2hilolo,v0p2lololo;
      double x0hihihi,x0lohihi,x0hilohi,x0lolohi;
      double x0hihilo,x0lohilo,x0hilolo,x0lololo;

      cudaMemcpy(&x0hihihi,&Ahihihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lohihi,&Alohihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0hilohi,&Ahilohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lolohi,&Alolohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0hihilo,&Ahihilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lohilo,&Alohilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0hilolo,&Ahilolo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lololo,&Alololo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      // mu = sqrt((*x0)*(*x0) + sigma[0]);
      odf_sqr(x0hihihi,x0lohihi,x0hilohi,x0lolohi,
              x0hihilo,x0lohilo,x0hilolo,x0lololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odf_inc(&acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo,
              sigmahihihi_h[0],sigmalohihi_h[0],
              sigmahilohi_h[0],sigmalolohi_h[0],
              sigmahihilo_h[0],sigmalohilo_h[0],
              sigmahilolo_h[0],sigmalololo_h[0]);
      odf_sqrt(acchihihi,acclohihi,acchilohi,acclolohi,
               acchihilo,acclohilo,acchilolo,acclololo,
               &muhihihi,&mulohihi,&muhilohi,&mulolohi,
               &muhihilo,&mulohilo,&muhilolo,&mulololo);
      if(x0hihihi <= 0.0)
      {
         // v0 = *x0 - mu;
         odf_sub( x0hihihi, x0lohihi, x0hilohi, x0lolohi,
                  x0hihilo, x0lohilo, x0hilolo, x0lololo,
                  muhihihi, mulohihi, muhilohi, mulolohi,
                  muhihilo, mulohilo, muhilolo, mulololo,
                 &v0hihihi,&v0lohihi,&v0hilohi,&v0lolohi,
                 &v0hihilo,&v0lohilo,&v0hilolo,&v0lololo);
      }
      else
      {
         // v0 = -sigma[0]/(*x0 + mu);
         odf_add(  x0hihihi,  x0lohihi,  x0hilohi,  x0lolohi,
                   x0hihilo,  x0lohilo,  x0hilolo,  x0lololo,
                   muhihihi,  mulohihi,  muhilohi,  mulolohi,
                   muhihilo,  mulohilo,  muhilolo,  mulololo,
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_div(sigmahihihi_h[0],sigmalohihi_h[0],
                 sigmahilohi_h[0],sigmalolohi_h[0],
                 sigmahihilo_h[0],sigmalohilo_h[0],
                 sigmahilolo_h[0],sigmalololo_h[0],
                   acchihihi,       acclohihi,     acchilohi,     acclolohi,
                   acchihilo,       acclohilo,     acchilolo,     acclololo,
                   &v0hihihi,       &v0lohihi,     &v0hilohi,     &v0lolohi,
                   &v0hihilo,       &v0lohilo,     &v0hilolo,     &v0lololo);
         odf_minus(&v0hihihi,&v0lohihi,&v0hilohi,&v0lolohi,
                   &v0hihilo,&v0lohilo,&v0hilolo,&v0lololo);
      }
      // v0p2 = v0*v0;
      odf_sqr(   v0hihihi,   v0lohihi,   v0hilohi,   v0lolohi,
                 v0hihilo,   v0lohilo,   v0hilolo,   v0lololo,
              &v0p2hihihi,&v0p2lohihi,&v0p2hilohi,&v0p2lolohi,
              &v0p2hihilo,&v0p2lohilo,&v0p2hilolo,&v0p2lololo);
      // *beta = 2.0*v0p2/(sigma[0] + v0p2);
      odf_add(sigmahihihi_h[0],sigmalohihi_h[0],
              sigmahilohi_h[0],sigmalolohi_h[0],
              sigmahihilo_h[0],sigmalohilo_h[0],
              sigmahilolo_h[0],sigmalololo_h[0],
               v0p2hihihi,      v0p2lohihi,      v0p2hilohi,      v0p2lolohi,
               v0p2hihilo,      v0p2lohilo,      v0p2hilolo,      v0p2lololo,
               &acchihihi,      &acclohihi,      &acchilohi,      &acclolohi,
               &acchihilo,      &acclohilo,      &acchilolo,      &acclololo);
      odf_div( v0p2hihihi,      v0p2lohihi,      v0p2hilohi,      v0p2lolohi,
               v0p2hihilo,      v0p2lohilo,      v0p2hilolo,      v0p2lololo,
                acchihihi,       acclohihi,       acchilohi,       acclolohi,
                acchihilo,       acclohilo,       acchilolo,       acclololo,
              &betahihihi_h[L],&betalohihi_h[L],
              &betahilohi_h[L],&betalolohi_h[L],
              &betahihilo_h[L],&betalohilo_h[L],
              &betahilolo_h[L],&betalololo_h[L]);
      odf_mlt_d(&betahihihi_h[L],&betalohihi_h[L],
                &betahihihi_h[L],&betalohihi_h[L],
                &betahilolo_h[L],&betalololo_h[L],
                &betahilolo_h[L],&betalololo_h[L],2.0);
      sigmahihihi_h[0] = v0hihihi;
      sigmalohihi_h[0] = v0lohihi;
      sigmahilohi_h[0] = v0hilohi;
      sigmalolohi_h[0] = v0lolohi;
      sigmahihilo_h[0] = v0hihilo;
      sigmalohilo_h[0] = v0lohilo;
      sigmahilolo_h[0] = v0hilolo;
      sigmalololo_h[0] = v0lololo;           // v0 needed for normalization
      // update the flop counts
      *add += 3;
      *mul += 3;
      *div += 2;
      *sqrtfun += 1;
   }
   if(verbose)
   {
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "          "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "          "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "          "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;
   }
   cudaMemcpy(&betahihihi_d[L],&betahihihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohihi_d[L],&betalohihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilohi_d[L],&betahilohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalolohi_d[L],&betalolohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahihilo_d[L],&betahihilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohilo_d[L],&betalohilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilolo_d[L],&betahilolo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalololo_d[L],&betalololo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);

   if(!done)  // normalization needed
   {
      // (sigmahi_h, sigmalo_h) has the values for (v0hi, v0lo).
      cudaMemcpy(sigmahihihi_d,sigmahihihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohihi_d,sigmalohihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilohi_d,sigmahilohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalolohi_d,sigmalolohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahihilo_d,sigmahihilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohilo_d,sigmalohilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilolo_d,sigmahilolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalololo_d,sigmalololo_h,sizeof(double),
                 cudaMemcpyHostToDevice);

      if(verbose)
      {
         cout << "-> launching " << nblocks << " blocks of "
              << szt << " threads to normalize ..." << endl;
         cout << "   nrows1 : " << nrows1
              << "  rowidx : " << rowidx << "  nVrows : " << nVrows << endl;
      }
      cudaEventRecord(start);
      dbl8_normalize<<<nblocks,szt>>>
         (nrows1,szt,
          &Ahihihi_d[rowidx+1],&Alohihi_d[rowidx+1],
          &Ahilohi_d[rowidx+1],&Alolohi_d[rowidx+1],
          &Ahihilo_d[rowidx+1],&Alohilo_d[rowidx+1],
          &Ahilolo_d[rowidx+1],&Alololo_d[rowidx+1],
          sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
          sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d,
          &Vhihihi_d[L*nVrows+L+1],&Vlohihi_d[L*nVrows+L+1],
          &Vhilohi_d[L*nVrows+L+1],&Vlolohi_d[L*nVrows+L+1],
          &Vhihilo_d[L*nVrows+L+1],&Vlohilo_d[L*nVrows+L+1],
          &Vhilolo_d[L*nVrows+L+1],&Vlololo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_normalize(nblocks,szt,div);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahihihi_h[L],&betahihihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohihi_h[L],&betalohihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilohi_h[L],&betahilohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalolohi_h[L],&betalolohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahihilo_h[L],&betahihilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohilo_h[L],&betalohilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilolo_h[L],&betahilolo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalololo_h[L],&betalololo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihihi_h,&Vhihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohihi_h,&Vlohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilohi_h,&Vhilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlolohi_h,&Vlolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihilo_h,&Vhihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohilo_h,&Vlohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilolo_h,&Vhilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlololo_h,&Vlololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);

      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "           "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "           "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "           "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : "
              << vhihihi_h[i] << "  " << vlohihi_h[i] << endl
              << "       "
              << vhilohi_h[i] << "  " << vlolohi_h[i] << endl
              << "       "
              << vhihilo_h[i] << "  " << vlohilo_h[i] << endl
              << "       "
              << vhilolo_h[i] << "  " << vlololo_h[i] << endl;
   }
}

void GPU_cmplx8_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *vrehihihi_h, double *vrelohihi_h,
   double *vrehilohi_h, double *vrelolohi_h,
   double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h,
   double *vimhilohi_h, double *vimlolohi_h,
   double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *sumshihihi_h, double *sumslohihi_h,
   double *sumshilohi_h, double *sumslolohi_h,
   double *sumshihilo_h, double *sumslohilo_h,
   double *sumshilolo_h, double *sumslololo_h,
   double *sumshihihi_d, double *sumslohihi_d,
   double *sumshilohi_d, double *sumslolohi_d,
   double *sumshihilo_d, double *sumslohilo_d,
   double *sumshilolo_d, double *sumslololo_d,
   double *sigmahihihi_h, double *sigmalohihi_h,
   double *sigmahilohi_h, double *sigmalolohi_h,
   double *sigmahihilo_h, double *sigmalohilo_h,
   double *sigmahilolo_h, double *sigmalololo_h,
   double *sigmahihihi_d, double *sigmalohihi_d,
   double *sigmahilohi_d, double *sigmalolohi_d,
   double *sigmahihilo_d, double *sigmalohilo_d,
   double *sigmahilolo_d, double *sigmalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
}

void GPU_dbl8_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   cudaEventRecord(start);           // 2nd argument: ncols -> szt
   // changed second argument ncols into szt
   // to avoid updating the next tile
   dbl8_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,
       Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
       Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
       &Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
       &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
       &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
       &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
       &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
       &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_leftRupdate(nrows,ncols,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Ahihihi_h,Ahihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohihi_h,Alohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilohi_h,Ahilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alolohi_h,Alolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahihilo_h,Ahihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohilo_h,Alohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilolo_h,Ahilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alololo_h,Alololo_d,sznum,cudaMemcpyDeviceToHost);

      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahihihi_h[j*nrows+i] << "  "
                 << Alohihi_h[j*nrows+i] << endl
                 << "            "
                 << Ahilohi_h[j*nrows+i] << "  "
                 << Alolohi_h[j*nrows+i] << endl
                 << "            "
                 << Ahihilo_h[j*nrows+i] << "  "
                 << Alohilo_h[j*nrows+i] << endl
                 << "            "
                 << Ahilolo_h[j*nrows+i] << "  "
                 << Alololo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx8_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
}

void GPU_dbl8_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *RTdotvhihihi_h, double *RTdotvlohihi_h,
   double *RTdotvhilohi_h, double *RTdotvlolohi_h,
   double *RTdotvhihilo_h, double *RTdotvlohilo_h,
   double *RTdotvhilolo_h, double *RTdotvlololo_h,
   double *RTdotvhihihi_d, double *RTdotvlohihi_d,
   double *RTdotvhilohi_d, double *RTdotvlolohi_d,
   double *RTdotvhihilo_d, double *RTdotvlohilo_d,
   double *RTdotvhilolo_d, double *RTdotvlololo_d,
   double *whihihi_h, double *wlohihi_h, double *whilohi_h, double *wlolohi_h,
   double *whihilo_h, double *wlohilo_h, double *whilolo_h, double *wlololo_h,
   double *whihihi_d, double *wlohihi_d, double *whilohi_d, double *wlolohi_d,
   double *whihilo_d, double *wlohilo_d, double *whilolo_d, double *wlololo_d,
   double *RTvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix
   const int nhouse = nrows - colidx;  // length of Householder vector
   // total number of entries in R that will be modified
   const int RToffset = colidx*nrows;
   const int dimRTdotv = endcol - colidx;
   // total number of entries in R that will be modified
   const int sizenum = (nrows - colidx)*dimRTdotv;
   const int nbrblocks = (int) ceil(sizenum/((double) szt));

   // changed second argument ncols into endcol
   // to avoid updating the next tile
   // dbl_medium_betaRTv<<<nbrblocks,szt>>>
   //   (nrows,endcol,szt,colidx,A_d,&V_d[L*nVrows+L],&beta_d[L],w_d);
   // number of threads must be ncols - colidx, not endcol - colidx
   // dbl2_small_betaRTv<<<1,nrows-colidx>>> // nrows ...
   //   (nrows,endcol,szt,colidx,Ahi_d,Alo_d,
   //    &Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],
   //    &betahi_d[L],&betalo_d[L],whi_d,wlo_d);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to compute RTdotv ..." << endl;
      cout << "   nhouse : " << nhouse << "  RToffset : " << RToffset
           << "  dimRTdotv : " << dimRTdotv << endl;
   }

   cudaEventRecord(start);
   dbl8_RTdotv<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RToffset,dimRTdotv,
       Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
       Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
       &Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
       &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
       &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
       &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
       RTdotvhihihi_d,RTdotvlohihi_d,RTdotvhilohi_d,RTdotvlolohi_d,
       RTdotvhihilo_d,RTdotvlohilo_d,RTdotvhilolo_d,RTdotvlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RTvlapms += milliseconds;
   cudaEventRecord(start);
   dbl8_sum_betaRTdotv<<<1,dimRTdotv>>>
      (nhouse,&betahihihi_d[L],&betalohihi_d[L],
              &betahilohi_d[L],&betalolohi_d[L],
              &betahihilo_d[L],&betalohilo_d[L],
              &betahilolo_d[L],&betalololo_d[L],
       RTdotvhihihi_d,RTdotvlohihi_d,RTdotvhilohi_d,RTdotvlolohi_d,
       RTdotvhihilo_d,RTdotvlohilo_d,RTdotvhilolo_d,RTdotvlololo_d,
       whihihi_d,wlohihi_d,whilohi_d,wlolohi_d,
       whihilo_d,wlohilo_d,whilolo_d,wlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RTvlapms += milliseconds;
   flopcount_dbl_RTdotv(nhouse,szt,mul);
   flopcount_dbl_sum_betaRTdotv(nhouse,dimRTdotv,add,mul);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to update " << sizenum << " numbers ..." << endl;
      cout << "   nrows : " << nrows << "  endcol : " << endcol
           << "  szt : " << szt << "  colidx : " << colidx << endl;
   }
   cudaEventRecord(start);
   dbl8_medium_subvbetaRTv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
       Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
       &Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
       &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
       &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
       &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
       &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
       &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L],
       whihihi_d,wlohihi_d,whilohi_d,wlolohi_d,
       whihilo_d,wlohilo_d,whilolo_d,wlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
   flopcount_dbl_medium_subvbetaRTv(nrows,endcol,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);
      const size_t szbRTv = dimRTdotv*sizeof(double);
      const size_t szRTdotv = nVrows*szbRTv;

      cudaMemcpy(RTdotvhihihi_h,RTdotvhihihi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlohihi_h,RTdotvlohihi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvhilohi_h,RTdotvhilohi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlolohi_h,RTdotvlolohi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvhihilo_h,RTdotvhihilo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlohilo_h,RTdotvlohilo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvhilolo_h,RTdotvhilolo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlololo_h,RTdotvlololo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);

      cout << "the matrix R^T dot v : " << endl;
      int ix = 0;
      for(int i=0; i<endcol-colidx; i++)
      {
         for(int j=0; j<nhouse; j++)      // must use nhouse
         {
            cout << "RTdotv[" << i << "][" << j << "] : "
                 << RTdotvhihihi_h[ix] << "  "
                 << RTdotvlohihi_h[ix] << endl
                 << "               "
                 << RTdotvhilohi_h[ix] << "  "
                 << RTdotvlolohi_h[ix] << endl
                 << "               "
                 << RTdotvhihilo_h[ix] << "  "
                 << RTdotvlohilo_h[ix] << endl
                 << "               "
                 << RTdotvhilolo_h[ix] << "  "
                 << RTdotvlololo_h[ix] << endl;
            ix = ix + 1;
         }
      }
      cudaMemcpy(whihihi_h,whihihi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlohihi_h,wlohihi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(whilohi_h,whilohi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlolohi_h,wlolohi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(whihilo_h,whihilo_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlohilo_h,wlohilo_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(whilolo_h,whilolo_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlololo_h,wlololo_d,szbRTv,cudaMemcpyDeviceToHost);

      cout << "the vector w = beta*R^T*v : " << endl;
      for(int i=0; i<endcol-colidx; i++)
         cout << "w[" << i << "] : "
              << whihihi_h[i] << "  " << wlohihi_h[i] << endl
              << "       "
              << whilohi_h[i] << "  " << wlolohi_h[i] << endl
              << "       "
              << whihilo_h[i] << "  " << wlohilo_h[i] << endl
              << "       "
              << whilolo_h[i] << "  " << wlololo_h[i] << endl;

      cudaMemcpy(Ahihihi_h,Ahihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohihi_h,Alohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilohi_h,Ahilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alolohi_h,Alolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahihilo_h,Ahihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohilo_h,Alohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilolo_h,Ahilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alololo_h,Alololo_d,sznum,cudaMemcpyDeviceToHost);

      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahihihi_h[j*nrows+i] << "  "
                 << Alohihi_h[j*nrows+i] << endl
                 << "            "
                 << Ahilohi_h[j*nrows+i] << "  "
                 << Alolohi_h[j*nrows+i] << endl
                 << "            "
                 << Ahihilo_h[j*nrows+i] << "  "
                 << Alohilo_h[j*nrows+i] << endl
                 << "            "
                 << Ahilolo_h[j*nrows+i] << "  "
                 << Alololo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx8_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *RHdotvrehihihi_h, double *RHdotvrelohihi_h,
   double *RHdotvrehilohi_h, double *RHdotvrelolohi_h,
   double *RHdotvrehihilo_h, double *RHdotvrelohilo_h,
   double *RHdotvrehilolo_h, double *RHdotvrelololo_h,
   double *RHdotvimhihihi_h, double *RHdotvimlohihi_h,
   double *RHdotvimhilohi_h, double *RHdotvimlolohi_h,
   double *RHdotvimhihilo_h, double *RHdotvimlohilo_h,
   double *RHdotvimhilolo_h, double *RHdotvimlololo_h,
   double *RHdotvrehihihi_d, double *RHdotvrelohihi_d,
   double *RHdotvrehilohi_d, double *RHdotvrelolohi_d,
   double *RHdotvrehihilo_d, double *RHdotvrelohilo_d,
   double *RHdotvrehilolo_d, double *RHdotvrelololo_d,
   double *RHdotvimhihihi_d, double *RHdotvimlohihi_d,
   double *RHdotvimhilohi_d, double *RHdotvimlolohi_d,
   double *RHdotvimhihilo_d, double *RHdotvimlohilo_d,
   double *RHdotvimhilolo_d, double *RHdotvimlololo_d,
   double *wrehihihi_h, double *wrelohihi_h,
   double *wrehilohi_h, double *wrelolohi_h,
   double *wrehihilo_h, double *wrelohilo_h,
   double *wrehilolo_h, double *wrelololo_h,
   double *wimhihihi_h, double *wimlohihi_h,
   double *wimhilohi_h, double *wimlolohi_h,
   double *wimhihilo_h, double *wimlohilo_h,
   double *wimhilolo_h, double *wimlololo_h,
   double *wrehihihi_d, double *wrelohihi_d,
   double *wrehilohi_d, double *wrelolohi_d,
   double *wrehihilo_d, double *wrelohilo_d,
   double *wrehilolo_d, double *wrelololo_d,
   double *wimhihihi_d, double *wimlohihi_d,
   double *wimhilohi_d, double *wimlolohi_d,
   double *wimhihilo_d, double *wimlohilo_d,
   double *wimhilolo_d, double *wimlololo_d,
   double *RHvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose )
{
}

void GPU_dbl8_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vhihihi_h, double *Vlohihi_h, double *Vhilohi_h, double *Vlolohi_h,
   double *Vhihilo_h, double *Vlohilo_h, double *Vhilolo_h, double *Vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *Whihihi_h, double *Wlohihi_h, double *Whilohi_h, double *Wlolohi_h,
   double *Whihilo_h, double *Wlohilo_h, double *Whilolo_h, double *Wlololo_h,
   double *Whihihi_d, double *Wlohihi_d, double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d, double *Whilolo_d, double *Wlololo_d,
   double *WYThihihi_h, double *WYTlohihi_h,
   double *WYThilohi_h, double *WYTlolohi_h,
   double *WYThihilo_h, double *WYTlohilo_h,
   double *WYThilolo_h, double *WYTlololo_h,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   const int nbrblocks1 = (int) ceil(rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_beta_times_V<<<nbrblocks1,szt>>>
      (rowdim,szt,betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                  betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                     Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
                     Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
                     Whihihi_d,   Wlohihi_d,   Whilohi_d,   Wlolohi_d,
                     Whihilo_d,   Wlohilo_d,   Whilolo_d,   Wlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_beta_times_V(rowdim,mul);

   const int nbrblocks2 = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_initialize_WYT<<<nbrblocks2,szt>>>
      (rowdim,szt,  Vhihihi_d,  Vlohihi_d,  Vhilohi_d,  Vlolohi_d,
                    Vhihilo_d,  Vlohilo_d,  Vhilolo_d,  Vlololo_d,
                    Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
                    Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
                  WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
                  WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_initialize_WYT(rowdim,mul);

   for(int j=1; j<szt; j++)
   {
      cudaEventRecord(start);
      dbl8_beta_next_W<<<nbrblocks1,szt>>>
         (rowdim,szt,
          &betahihihi_d[j],&betalohihi_d[j],&betahilohi_d[j],&betalolohi_d[j],
          &betahihilo_d[j],&betalohilo_d[j],&betahilolo_d[j],&betalololo_d[j],
          &Vhihihi_d[j*rowdim],&Vlohihi_d[j*rowdim],
          &Vhilohi_d[j*rowdim],&Vlolohi_d[j*rowdim],
          &Vhihilo_d[j*rowdim],&Vlohilo_d[j*rowdim],
          &Vhilolo_d[j*rowdim],&Vlololo_d[j*rowdim],
          &Whihihi_d[j*rowdim],&Wlohihi_d[j*rowdim],
          &Whilohi_d[j*rowdim],&Wlolohi_d[j*rowdim],
          &Whihilo_d[j*rowdim],&Wlohilo_d[j*rowdim],
          &Whilolo_d[j*rowdim],&Wlololo_d[j*rowdim],
          WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
          WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_beta_next_W(rowdim,add,mul);

      cudaEventRecord(start);
      dbl8_update_WYT<<<nbrblocks2,szt>>>
         (rowdim,szt,
          &Vhihihi_d[j*rowdim],&Vlohihi_d[j*rowdim],
          &Vhilohi_d[j*rowdim],&Vlolohi_d[j*rowdim],
          &Vhihilo_d[j*rowdim],&Vlohilo_d[j*rowdim],
          &Vhilolo_d[j*rowdim],&Vlololo_d[j*rowdim],
          &Whihihi_d[j*rowdim],&Wlohihi_d[j*rowdim],
          &Whilohi_d[j*rowdim],&Wlolohi_d[j*rowdim],
          &Whihilo_d[j*rowdim],&Wlohilo_d[j*rowdim],
          &Whilolo_d[j*rowdim],&Wlololo_d[j*rowdim],
          WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
          WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_update_WYT(rowdim,add,mul);
   }
   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = rowdim*sizeof(double);
      const size_t szVandW = szt*szhouse;
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(betahihihi_h,betahihihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohihi_h,betalohihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilohi_h,betahilohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalolohi_h,betalolohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahihilo_h,betahihilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohilo_h,betalohilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilolo_h,betahilolo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalololo_h,betalololo_d,szbeta,cudaMemcpyDeviceToHost);

      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : "
              << betahihihi_h[j] << "  " << betalohihi_h[j] << endl
              << "          "
              << betahilohi_h[j] << "  " << betalolohi_h[j] << endl
              << "          "
              << betahihilo_h[j] << "  " << betalohilo_h[j] << endl
              << "          "
              << betahilolo_h[j] << "  " << betalololo_h[j] << endl;

      cudaMemcpy(Vhihihi_h,Vhihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlohihi_h,Vlohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vhilohi_h,Vhilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlolohi_h,Vlolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vhihilo_h,Vhihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlohilo_h,Vlohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vhilolo_h,Vhilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlololo_h,Vlololo_d,szVandW,cudaMemcpyDeviceToHost);

      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "V[" << i << "][" << j << "] : "
                 << Vhihihi_h[ix] << "  " << Vlohihi_h[ix] << endl
                 << "          "
                 << Vhilohi_h[ix] << "  " << Vlolohi_h[ix] << endl
                 << "          "
                 << Vhihilo_h[ix] << "  " << Vlohilo_h[ix] << endl
                 << "          "
                 << Vhilolo_h[ix] << "  " << Vlololo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(Whihihi_h,Whihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlohihi_h,Wlohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Whilohi_h,Whilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlolohi_h,Wlolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Whihilo_h,Whihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlohilo_h,Wlohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Whilolo_h,Whilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlololo_h,Wlololo_d,szVandW,cudaMemcpyDeviceToHost);

      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "W[" << i << "][" << j << "] : "
                 << Whihihi_h[ix] << "  " << Wlohihi_h[ix] << endl
                 << "          "
                 << Whilohi_h[ix] << "  " << Wlolohi_h[ix] << endl
                 << "          "
                 << Whihilo_h[ix] << "  " << Wlohilo_h[ix] << endl
                 << "          "
                 << Whilolo_h[ix] << "  " << Wlololo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(WYThihihi_h,WYThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohihi_h,WYTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilohi_h,WYThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlolohi_h,WYTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThihilo_h,WYThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohilo_h,WYTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilolo_h,WYThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlololo_h,WYTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihihi_h[ix] << "  " << WYTlohihi_h[ix] << endl
                 << "            "
                 << WYThilohi_h[ix] << "  " << WYTlolohi_h[ix] << endl
                 << "            "
                 << WYThihilo_h[ix] << "  " << WYTlohilo_h[ix] << endl
                 << "            "
                 << WYThilolo_h[ix] << "  " << WYTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vrehihihi_h, double *Vrelohihi_h,
   double *Vrehilohi_h, double *Vrelolohi_h,
   double *Vrehihilo_h, double *Vrelohilo_h,
   double *Vrehilolo_h, double *Vrelololo_h,
   double *Vimhihihi_h, double *Vimlohihi_h,
   double *Vimhilohi_h, double *Vimlolohi_h,
   double *Vimhihilo_h, double *Vimlohilo_h,
   double *Vimhilolo_h, double *Vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d, 
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *Wrehihihi_h, double *Wrelohihi_h,
   double *Wrehilohi_h, double *Wrelolohi_h,
   double *Wrehihilo_h, double *Wrelohilo_h,
   double *Wrehilolo_h, double *Wrelololo_h,
   double *Wimhihihi_h, double *Wimlohihi_h,
   double *Wimhilohi_h, double *Wimlolohi_h,
   double *Wimhihilo_h, double *Wimlohilo_h,
   double *Wimhilolo_h, double *Wimlololo_h,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *WYHrehihihi_h, double *WYHrelohihi_h,
   double *WYHrehilohi_h, double *WYHrelolohi_h,
   double *WYHrehihilo_h, double *WYHrelohilo_h,
   double *WYHrehilolo_h, double *WYHrelololo_h,
   double *WYHimhihihi_h, double *WYHimlohihi_h,
   double *WYHimhilohi_h, double *WYHimlolohi_h,
   double *WYHimhihilo_h, double *WYHimlohilo_h,
   double *WYHimhilolo_h, double *WYHimlololo_h,
   double *WYHrehihihi_d, double *WYHrelohihi_d,
   double *WYHrehilohi_d, double *WYHrelolohi_d,
   double *WYHrehihilo_d, double *WYHrelohilo_d,
   double *WYHrehilolo_d, double *WYHrelololo_d,
   double *WYHimhihihi_d, double *WYHimlohihi_d,
   double *WYHimhilohi_d, double *WYHimlolohi_d,
   double *WYHimhihilo_d, double *WYHimlohilo_d,
   double *WYHimhilolo_d, double *WYHimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
}

void GPU_dbl8_small_WYT
 ( int nrows, int szt,
   double *Whihihi_d, double *Wlohihi_d,
   double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d,
   double *Whilolo_d, double *Wlololo_d,
   double *Yhihihi_d, double *Ylohihi_d,
   double *Yhilohi_d, double *Ylolohi_d,
   double *Yhihilo_d, double *Ylohilo_d,
   double *Yhilolo_d, double *Ylololo_d,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *WYThihihi_h, double *WYTlohihi_h,
   double *WYThilohi_h, double *WYTlolohi_h,
   double *WYThihilo_h, double *WYTlohilo_h,
   double *WYThilolo_h, double *WYTlololo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   dbl8_small_WYT<<<nbrblocks,szt>>>
      (nrows,szt,  Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
                   Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
                   Yhihihi_d,  Ylohihi_d,  Yhilohi_d,  Ylolohi_d,
                   Yhihilo_d,  Ylohilo_d,  Yhilolo_d,  Ylololo_d,
                 WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
                 WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   // flopcount_dbl_small_WYT(nrows,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYThihihi_h,WYThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohihi_h,WYTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilohi_h,WYThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlolohi_h,WYTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThihilo_h,WYThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohilo_h,WYTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilolo_h,WYThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlololo_h,WYTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihihi_h[ix] << "  " << WYTlohihi_h[ix] << endl
                 << "            "
                 << WYThilohi_h[ix] << "  " << WYTlolohi_h[ix] << endl
                 << "            "
                 << WYThihilo_h[ix] << "  " << WYTlohilo_h[ix] << endl
                 << "            "
                 << WYThilolo_h[ix] << "  " << WYTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_WYH
 ( int nrows, int szt,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *Yrehihihi_d, double *Yrelohihi_d,
   double *Yrehilohi_d, double *Yrelolohi_d,
   double *Yrehihilo_d, double *Yrelohilo_d,
   double *Yrehilolo_d, double *Yrelololo_d,
   double *Yimhihihi_d, double *Yimlohihi_d,
   double *Yimhilohi_d, double *Yimlolohi_d,
   double *Yimhihilo_d, double *Yimlohilo_d,
   double *Yimhilolo_d, double *Yimlololo_d,
   double *WYTrehihihi_d, double *WYTrelohihi_d,
   double *WYTrehilohi_d, double *WYTrelolohi_d,
   double *WYTrehihilo_d, double *WYTrelohilo_d,
   double *WYTrehilolo_d, double *WYTrelololo_d,
   double *WYTimhihihi_d, double *WYTimlohihi_d,
   double *WYTimhilohi_d, double *WYTimlolohi_d,
   double *WYTimhihilo_d, double *WYTimlohilo_d,
   double *WYTimhilolo_d, double *WYTimlololo_d,
   double *WYTrehihihi_h, double *WYTrelohihi_h,
   double *WYTrehilohi_h, double *WYTrelolohi_h,
   double *WYTrehihilo_h, double *WYTrelohilo_h,
   double *WYTrehilolo_h, double *WYTrelololo_h,
   double *WYTimhihihi_h, double *WYTimlohihi_h,
   double *WYTimhilohi_h, double *WYTimlolohi_h,
   double *WYTimhihilo_h, double *WYTimlohilo_h,
   double *WYTimhilolo_h, double *WYTimlololo_h,
   double *lapms, bool verbose )
{
}

void GPU_dbl8_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhihihi_d, double *Ylohihi_d, double *Yhilohi_d, double *Ylolohi_d,
   double *Yhihilo_d, double *Ylohilo_d, double *Yhilolo_d, double *Ylololo_d,
   double *Whihihi_d, double *Wlohihi_d, double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d, double *Whilolo_d, double *Wlololo_d,
   double *YWThihihi_d, double *YWTlohihi_d,
   double *YWThilohi_d, double *YWTlolohi_d,
   double *YWThihilo_d, double *YWTlohilo_d,
   double *YWThilolo_d, double *YWTlololo_d,
   double *YWThihihi_h, double *YWTlohihi_h,
   double *YWThilohi_h, double *YWTlolohi_h,
   double *YWThihilo_h, double *YWTlohilo_h,
   double *YWThilolo_h, double *YWTlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_small_WYT<<<nbrblocks,szt>>>
      (rowdim,szt,  Yhihihi_d,  Ylohihi_d,  Yhilohi_d,  Ylolohi_d,
                    Yhihilo_d,  Ylohilo_d,  Yhilolo_d,  Ylololo_d,
                    Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
                    Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
                  YWThihihi_d,YWTlohihi_d,YWThilohi_d,YWTlolohi_d,
                  YWThihilo_d,YWTlohilo_d,YWThilolo_d,YWTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_WYT(rowdim,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWThihihi_h,YWThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlohihi_h,YWTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWThilohi_h,YWThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlolohi_h,YWTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWThihilo_h,YWThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlohilo_h,YWTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWThilolo_h,YWThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlololo_h,YWTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWThihihi_h[ix] << "  " << YWTlohihi_h[ix] << endl
                 << "            "
                 << YWThilohi_h[ix] << "  " << YWTlolohi_h[ix] << endl
                 << "            "
                 << YWThihilo_h[ix] << "  " << YWTlohilo_h[ix] << endl
                 << "            "
                 << YWThilolo_h[ix] << "  " << YWTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_YWH
 ( int nrows, int szt, int idx,
   double *Yrehihihi_d, double *Yrelohihi_d,
   double *Yrehilohi_d, double *Yrelolohi_d,
   double *Yrehihilo_d, double *Yrelohilo_d,
   double *Yrehilolo_d, double *Yrelololo_d,
   double *Yimhihihi_d, double *Yimlohihi_d,
   double *Yimhilohi_d, double *Yimlolohi_d,
   double *Yimhihilo_d, double *Yimlohilo_d,
   double *Yimhilolo_d, double *Yimlololo_d,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *YWTrehihihi_d, double *YWTrelohihi_d,
   double *YWTrehilohi_d, double *YWTrelolohi_d,
   double *YWTrehihilo_d, double *YWTrelohilo_d,
   double *YWTrehilolo_d, double *YWTrelololo_d,
   double *YWTimhihihi_d, double *YWTimlohihi_d,
   double *YWTimhilohi_d, double *YWTimlolohi_d,
   double *YWTimhihilo_d, double *YWTimlohilo_d,
   double *YWTimhilolo_d, double *YWTimlololo_d,
   double *YWTrehihihi_h, double *YWTrelohihi_h,
   double *YWTrehilohi_h, double *YWTrelolohi_h,
   double *YWTrehihilo_h, double *YWTrelohilo_h,
   double *YWTrehilolo_h, double *YWTrelololo_h,
   double *YWTimhihihi_h, double *YWTimlohihi_h,
   double *YWTimhilohi_h, double *YWTimlolohi_h,
   double *YWTimhihilo_h, double *YWTimlohilo_h,
   double *YWTimhilolo_h, double *YWTimlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
}

void GPU_dbl8_small_QWYT
 ( int dim, int szt, int idx,
   double *Qhihihi_d, double *Qlohihi_d, double *Qhilohi_d, double *Qlolohi_d,
   double *Qhihilo_d, double *Qlohilo_d, double *Qhilolo_d, double *Qlololo_d,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *QWYThihihi_d, double *QWYTlohihi_d,
   double *QWYThilohi_d, double *QWYTlolohi_d,
   double *QWYThihilo_d, double *QWYTlohilo_d,
   double *QWYThilolo_d, double *QWYTlololo_d,
   double *QWYThihihi_h, double *QWYTlohihi_h,
   double *QWYThilohi_h, double *QWYTlolohi_h,
   double *QWYThihilo_h, double *QWYTlohilo_h,
   double *QWYThilolo_h, double *QWYTlololo_h,
   double *Qhihihi_h, double *Qlohihi_h, double *Qhilohi_h, double *Qlolohi_h,
   double *Qhihilo_h, double *Qlohilo_h, double *Qhilolo_h, double *Qlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qhihihi_h,Qhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohihi_h,Qlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilohi_h,Qhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlolohi_h,Qlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhihilo_h,Qhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohilo_h,Qlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilolo_h,Qhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlololo_h,Qlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhihihi_h[ix] << "  " << Qlohihi_h[ix] << endl
                 << "          "
                 << Qhilohi_h[ix] << "  " << Qlolohi_h[ix] << endl
                 << "          "
                 << Qhihilo_h[ix] << "  " << Qlohilo_h[ix] << endl
                 << "          "
                 << Qhilolo_h[ix] << "  " << Qlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   dbl8_small_QWYT<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
          Qhihihi_d,   Qlohihi_d,   Qhilohi_d,   Qlolohi_d,
          Qhihilo_d,   Qlohilo_d,   Qhilolo_d,   Qlololo_d,
        WYThihihi_d, WYTlohihi_d, WYThilohi_d, WYTlolohi_d,
        WYThihilo_d, WYTlohilo_d, WYThilolo_d, WYTlololo_d,
       QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
       QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_QWYT(dim,rowdim,szt,coloff,add,mul);

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYThihihi_h,QWYThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlohihi_h,QWYTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYThilohi_h,QWYThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlolohi_h,QWYTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYThihilo_h,QWYThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlohilo_h,QWYTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYThilolo_h,QWYThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlololo_h,QWYTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "QWYT[" << i << "][" << j << "] : "
                 << QWYThihihi_h[ix] << "  " << QWYTlohihi_h[ix] << endl
                 << "             "
                 << QWYThilohi_h[ix] << "  " << QWYTlolohi_h[ix] << endl
                 << "             "
                 << QWYThihilo_h[ix] << "  " << QWYTlohilo_h[ix] << endl
                 << "             "
                 << QWYThilolo_h[ix] << "  " << QWYTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_QWYH
 ( int dim, int szt, int idx,
   double *Qrehihihi_d, double *Qrelohihi_d,
   double *Qrehilohi_d, double *Qrelolohi_d,
   double *Qrehihilo_d, double *Qrelohilo_d,
   double *Qrehilolo_d, double *Qrelololo_d,
   double *Qimhihihi_d, double *Qimlohihi_d,
   double *Qimhilohi_d, double *Qimlolohi_d,
   double *Qimhihilo_d, double *Qimlohilo_d,
   double *Qimhilolo_d, double *Qimlololo_d,
   double *WYTrehihihi_d, double *WYTrelohihi_d,
   double *WYTrehilohi_d, double *WYTrelolohi_d,
   double *WYTrehihilo_d, double *WYTrelohilo_d,
   double *WYTrehilolo_d, double *WYTrelololo_d,
   double *WYTimhihihi_d, double *WYTimlohihi_d,
   double *WYTimhilohi_d, double *WYTimlolohi_d,
   double *WYTimhihilo_d, double *WYTimlohilo_d,
   double *WYTimhilolo_d, double *WYTimlololo_d,
   double *QWYTrehihihi_d, double *QWYTrelohihi_d,
   double *QWYTrehilohi_d, double *QWYTrelolohi_d,
   double *QWYTrehihilo_d, double *QWYTrelohilo_d,
   double *QWYTrehilolo_d, double *QWYTrelololo_d,
   double *QWYTimhihihi_d, double *QWYTimlohihi_d,
   double *QWYTimhilohi_d, double *QWYTimlolohi_d,
   double *QWYTimhihilo_d, double *QWYTimlohilo_d,
   double *QWYTimhilolo_d, double *QWYTimlololo_d,
   double *QWYTrehihihi_h, double *QWYTrelohihi_h,
   double *QWYTrehilohi_h, double *QWYTrelolohi_h,
   double *QWYTrehihilo_h, double *QWYTrelohilo_h,
   double *QWYTrehilolo_h, double *QWYTrelololo_h,
   double *QWYTimhihihi_h, double *QWYTimlohihi_h,
   double *QWYTimhilohi_h, double *QWYTimlolohi_h,
   double *QWYTimhihilo_h, double *QWYTimlohilo_h,
   double *QWYTimhilolo_h, double *QWYTimlololo_h,
   double *Qrehihihi_h, double *Qrelohihi_h,
   double *Qrehilohi_h, double *Qrelolohi_h,
   double *Qrehihilo_h, double *Qrelohilo_h,
   double *Qrehilolo_h, double *Qrelololo_h,
   double *Qimhihihi_h, double *Qimlohihi_h,
   double *Qimhilohi_h, double *Qimlolohi_h,
   double *Qimhihilo_h, double *Qimlohilo_h,
   double *Qimhilolo_h, double *Qimlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
}

void GPU_dbl8_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWThihihi_d, double *YWTlohihi_d,
   double *YWThilohi_d, double *YWTlolohi_d,
   double *YWThihilo_d, double *YWTlohilo_d,
   double *YWThilolo_d, double *YWTlololo_d,
   double *Chihihi_d, double *Clohihi_d, double *Chilohi_d, double *Clolohi_d,
   double *Chihilo_d, double *Clohilo_d, double *Chilolo_d, double *Clololo_d,
   double *YWTChihihi_d, double *YWTClohihi_d,
   double *YWTChilohi_d, double *YWTClolohi_d,
   double *YWTChihilo_d, double *YWTClohilo_d,
   double *YWTChilolo_d, double *YWTClololo_d,
   double *YWTChihihi_h, double *YWTClohihi_h,
   double *YWTChilohi_h, double *YWTClolohi_h,
   double *YWTChihilo_h, double *YWTClohilo_h,
   double *YWTChilolo_h, double *YWTClololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   if(verbose)
   {
      cout << "in GPU_dbl2_small_YWTC ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff
           << "  nbrblocks : " << nbrblocks << endl;

      double *Chihihi_h = new double[nrows*ncols];
      double *Clohihi_h = new double[nrows*ncols];
      double *Chilohi_h = new double[nrows*ncols];
      double *Clolohi_h = new double[nrows*ncols];
      double *Chihilo_h = new double[nrows*ncols];
      double *Clohilo_h = new double[nrows*ncols];
      double *Chilolo_h = new double[nrows*ncols];
      double *Clololo_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Chihihi_h,Chihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clohihi_h,Clohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Chilohi_h,Chilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clolohi_h,Clolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Chihilo_h,Chihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clohilo_h,Clohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Chilolo_h,Chilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clololo_h,Clololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
            cout << "C_h[" << i << "][" << j << "] : "
                 << Chihihi_h[j*nrows+i] << "  "
                 << Clohihi_h[j*nrows+i] << endl
                 << "            "
                 << Chilohi_h[j*nrows+i] << "  "
                 << Clolohi_h[j*nrows+i] << endl
                 << "            "
                 << Chihilo_h[j*nrows+i] << "  "
                 << Clohilo_h[j*nrows+i] << endl
                 << "            "
                 << Chilolo_h[j*nrows+i] << "  "
                 << Clololo_h[j*nrows+i] << endl;

      free(Chihihi_h); free(Clohihi_h); free(Chilohi_h); free(Clolohi_h);
      free(Chihilo_h); free(Clohilo_h); free(Chilolo_h); free(Clololo_h);
   }

   cudaEventRecord(start);
   dbl8_small_YWTC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
        YWThihihi_d, YWTlohihi_d, YWThilohi_d, YWTlolohi_d,
        YWThihilo_d, YWTlohilo_d, YWThilolo_d, YWTlololo_d,
          Chihihi_d,   Clohihi_d,   Chilohi_d,   Clolohi_d,
          Chihilo_d,   Clohilo_d,   Chilolo_d,   Clololo_d,
       YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
       YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_YWTC(rowdim,coldim,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTChihihi_h,YWTChihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClohihi_h,YWTClohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTChilohi_h,YWTChilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClolohi_h,YWTClolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTChihilo_h,YWTChihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClohilo_h,YWTClohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTChilolo_h,YWTChilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClololo_h,YWTClololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "YWTC[" << i << "][" << j << "] : "
                 << YWTChihihi_h[j*nrows + i] << "  "
                 << YWTClohihi_h[j*nrows + i] << endl
                 << "             "
                 << YWTChilohi_h[j*nrows + i] << "  "
                 << YWTClolohi_h[j*nrows + i] << endl
                 << "             "
                 << YWTChihilo_h[j*nrows + i] << "  "
                 << YWTClohilo_h[j*nrows + i] << endl
                 << "             "
                 << YWTChilolo_h[j*nrows + i] << "  "
                 << YWTClololo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx8_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehihihi_d, double *YWTrelohihi_d,
   double *YWTrehilohi_d, double *YWTrelolohi_d,
   double *YWTrehihilo_d, double *YWTrelohilo_d,
   double *YWTrehilolo_d, double *YWTrelololo_d,
   double *YWTimhihihi_d, double *YWTimlohihi_d,
   double *YWTimhilohi_d, double *YWTimlolohi_d,
   double *YWTimhihilo_d, double *YWTimlohilo_d,
   double *YWTimhilolo_d, double *YWTimlololo_d,
   double *Crehihihi_d, double *Crelohihi_d,
   double *Crehilohi_d, double *Crelolohi_d,
   double *Crehihilo_d, double *Crelohilo_d,
   double *Crehilolo_d, double *Crelololo_d,
   double *Cimhihihi_d, double *Cimlohihi_d,
   double *Cimhilohi_d, double *Cimlolohi_d,
   double *Cimhihilo_d, double *Cimlohilo_d,
   double *Cimhilolo_d, double *Cimlololo_d,
   double *YWTCrehihihi_d, double *YWTCrelohihi_d,
   double *YWTCrehilohi_d, double *YWTCrelolohi_d,
   double *YWTCrehihilo_d, double *YWTCrelohilo_d,
   double *YWTCrehilolo_d, double *YWTCrelololo_d,
   double *YWTCimhihihi_d, double *YWTCimlohihi_d,
   double *YWTCimhilohi_d, double *YWTCimlolohi_d,
   double *YWTCimhihilo_d, double *YWTCimlohilo_d,
   double *YWTCimhilolo_d, double *YWTCimlololo_d,
   double *YWTCrehihihi_h, double *YWTCrelohihi_h,
   double *YWTCrehilohi_h, double *YWTCrelolohi_h,
   double *YWTCrehihilo_h, double *YWTCrelohilo_h,
   double *YWTCrehilolo_h, double *YWTCrelololo_h,
   double *YWTCimhihihi_h, double *YWTCimlohihi_h,
   double *YWTCimhilohi_h, double *YWTCimlolohi_h,
   double *YWTCimhihilo_h, double *YWTCimlohilo_h,
   double *YWTCimhilolo_h, double *YWTCimlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
}

void GPU_dbl8_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhihihi_d, double *Qlohihi_d, double *Qhilohi_d, double *Qlolohi_d,
   double *Qhihilo_d, double *Qlohilo_d, double *Qhilolo_d, double *Qlololo_d,
   double *QWYThihihi_d, double *QWYTlohihi_d,
   double *QWYThilohi_d, double *QWYTlolohi_d,
   double *QWYThihilo_d, double *QWYTlohilo_d,
   double *QWYThilolo_d, double *QWYTlololo_d,
   double *Qhihihi_h, double *Qlohihi_h, double *Qhilohi_h, double *Qlolohi_h,
   double *Qhihilo_h, double *Qlohilo_h, double *Qhilolo_h, double *Qlololo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
       Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
       Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
       QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
       QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_Qupdate(dim,rowdim,add);

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qhihihi_h,Qhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohihi_h,Qlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilohi_h,Qhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlolohi_h,Qlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhihilo_h,Qhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohilo_h,Qlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilolo_h,Qhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlololo_h,Qlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhihihi_h[ix] << "  " << Qlohihi_h[ix] << endl
                 << "          "
                 << Qhilohi_h[ix] << "  " << Qlolohi_h[ix] << endl
                 << "          "
                 << Qhihilo_h[ix] << "  " << Qlohilo_h[ix] << endl
                 << "          "
                 << Qhilolo_h[ix] << "  " << Qlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qrehihihi_d, double *Qrelohihi_d,
   double *Qrehilohi_d, double *Qrelolohi_d,
   double *Qrehihilo_d, double *Qrelohilo_d,
   double *Qrehilolo_d, double *Qrelololo_d,
   double *Qimhihihi_d, double *Qimlohihi_d,
   double *Qimhilohi_d, double *Qimlolohi_d,
   double *Qimhihilo_d, double *Qimlohilo_d,
   double *Qimhilolo_d, double *Qimlololo_d,
   double *QWYTrehihihi_d, double *QWYTrelohihi_d,
   double *QWYTrehilohi_d, double *QWYTrelolohi_d,
   double *QWYTrehihilo_d, double *QWYTrelohilo_d,
   double *QWYTrehilolo_d, double *QWYTrelololo_d,
   double *QWYTimhihihi_d, double *QWYTimlohihi_d,
   double *QWYTimhilohi_d, double *QWYTimlolohi_d,
   double *QWYTimhihilo_d, double *QWYTimlohilo_d,
   double *QWYTimhilolo_d, double *QWYTimlololo_d,
   double *Qrehihihi_h, double *Qrelohihi_h,
   double *Qrehilohi_h, double *Qrelolohi_h,
   double *Qrehihilo_h, double *Qrelohilo_h,
   double *Qrehilolo_h, double *Qrelololo_h,
   double *Qimhihihi_h, double *Qimlohihi_h,
   double *Qimhilohi_h, double *Qimlolohi_h,
   double *Qimhihilo_h, double *Qimlohilo_h,
   double *Qimhilolo_h, double *Qimlololo_h,
   double *lapms, long long int *add, bool verbose )
{
}

void GPU_dbl8_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rhihihi_d, double *Rlohihi_d, double *Rhilohi_d, double *Rlolohi_d,
   double *Rhihilo_d, double *Rlohilo_d, double *Rhilolo_d, double *Rlololo_d,
   double *YWTChihihi_d, double *YWTClohihi_d,
   double *YWTChilohi_d, double *YWTClolohi_d,
   double *YWTChihilo_d, double *YWTClohilo_d,
   double *YWTChilolo_d, double *YWTClololo_d,
   double *Rhihihi_h, double *Rlohihi_h, double *Rhilohi_h, double *Rlolohi_h,
   double *Rhihilo_h, double *Rlohilo_h, double *Rhilolo_h, double *Rlololo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   cudaEventRecord(start);
   dbl8_small_R_add_YWTC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,
       Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
       Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
       YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
       YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_R_add_YWTC(nrows,coldim,szt,rowoff,coloff,add);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rhihihi_h,Rhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlohihi_h,Rlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rhilohi_h,Rhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlolohi_h,Rlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rhihilo_h,Rhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlohilo_h,Rlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rhilolo_h,Rhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlololo_h,Rlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhihihi_h[j*nrows + i] << "  "
                 << Rlohihi_h[j*nrows + i] << endl
                 << "          "
                 << Rhilohi_h[j*nrows + i] << "  "
                 << Rlolohi_h[j*nrows + i] << endl
                 << "          "
                 << Rhihilo_h[j*nrows + i] << "  "
                 << Rlohilo_h[j*nrows + i] << endl
                 << "          "
                 << Rhilolo_h[j*nrows + i] << "  "
                 << Rlololo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx8_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehihihi_d, double *Rrelohihi_d,
   double *Rrehilohi_d, double *Rrelolohi_d,
   double *Rrehihilo_d, double *Rrelohilo_d,
   double *Rrehilolo_d, double *Rrelololo_d,
   double *Rimhihihi_d, double *Rimlohihi_d,
   double *Rimhilohi_d, double *Rimlolohi_d,
   double *Rimhihilo_d, double *Rimlohilo_d,
   double *Rimhilolo_d, double *Rimlololo_d,
   double *YWTCrehihihi_d, double *YWTCrelohihi_d,
   double *YWTCrehilohi_d, double *YWTCrelolohi_d,
   double *YWTCrehihilo_d, double *YWTCrelohilo_d,
   double *YWTCrehilolo_d, double *YWTCrelololo_d,
   double *YWTCimhihihi_d, double *YWTCimlohihi_d,
   double *YWTCimhilohi_d, double *YWTCimlolohi_d,
   double *YWTCimhihilo_d, double *YWTCimlohilo_d,
   double *YWTCimhilolo_d, double *YWTCimlololo_d,
   double *Rrehihihi_h, double *Rrelohihi_h,
   double *Rrehilohi_h, double *Rrelolohi_h,
   double *Rrehihilo_h, double *Rrelohilo_h,
   double *Rrehilolo_h, double *Rrelololo_h,
   double *Rimhihihi_h, double *Rimlohihi_h,
   double *Rimhilohi_h, double *Rimlolohi_h,
   double *Rimhihilo_h, double *Rimlohilo_h,
   double *Rimhilolo_h, double *Rimlololo_h,
   double *lapms, long long int *add, bool verbose )
{
}

void GPU_dbl8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *houselapms, double *RTvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;          // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Ahihihi_h = new double[dim];    // A on the host
   double *Alohihi_h = new double[dim]; 
   double *Ahilohi_h = new double[dim];
   double *Alolohi_h = new double[dim]; 
   double *Ahihilo_h = new double[dim]; 
   double *Alohilo_h = new double[dim]; 
   double *Ahilolo_h = new double[dim];
   double *Alololo_h = new double[dim]; 
   double *Ahihihi_d;                      // A on the device
   double *Alohihi_d; 
   double *Ahilohi_d; 
   double *Alolohi_d; 
   double *Ahihilo_d;
   double *Alohilo_d; 
   double *Ahilolo_d; 
   double *Alololo_d; 
   double *Qhihihi_h = new double[nrows2]; // Q on the host
   double *Qlohihi_h = new double[nrows2]; 
   double *Qhilohi_h = new double[nrows2]; 
   double *Qlolohi_h = new double[nrows2]; 
   double *Qhihilo_h = new double[nrows2];
   double *Qlohilo_h = new double[nrows2]; 
   double *Qhilolo_h = new double[nrows2]; 
   double *Qlololo_h = new double[nrows2]; 
   double *Qhihihi_d;                      // Q on the device
   double *Qlohihi_d;
   double *Qhilohi_d;
   double *Qlolohi_d;
   double *Qhihilo_d;                      // Q on the device
   double *Qlohilo_d;
   double *Qhilolo_d;
   double *Qlololo_d;
   double *vhihihi_h = new double[nrows];  // Householder vector
   double *vlohihi_h = new double[nrows];
   double *vhilohi_h = new double[nrows];
   double *vlolohi_h = new double[nrows];
   double *vhihilo_h = new double[nrows]; 
   double *vlohilo_h = new double[nrows];
   double *vhilolo_h = new double[nrows];
   double *vlololo_h = new double[nrows];
   double *betahihihi_h = new double[szt]; //  beta on the host
   double *betalohihi_h = new double[szt]; 
   double *betahilohi_h = new double[szt]; 
   double *betalolohi_h = new double[szt]; 
   double *betahihilo_h = new double[szt]; 
   double *betalohilo_h = new double[szt]; 
   double *betahilolo_h = new double[szt]; 
   double *betalololo_h = new double[szt]; 
   double *betahihihi_d;                   // beta on the device
   double *betalohihi_d;
   double *betahilohi_d;
   double *betalolohi_d;
   double *betahihilo_d;
   double *betalohilo_d;
   double *betahilolo_d;
   double *betalololo_d;
   double *Vhihihi_h = new double[nrows*szt]; // V matrix
   double *Vlohihi_h = new double[nrows*szt];
   double *Vhilohi_h = new double[nrows*szt];
   double *Vlolohi_h = new double[nrows*szt];
   double *Vhihilo_h = new double[nrows*szt];
   double *Vlohilo_h = new double[nrows*szt];
   double *Vhilolo_h = new double[nrows*szt];
   double *Vlololo_h = new double[nrows*szt];
   double *Vhihihi_d;                         // V on the device
   double *Vlohihi_d;
   double *Vhilohi_d;
   double *Vlolohi_d;
   double *Vhihilo_d;
   double *Vlohilo_d;
   double *Vhilolo_d;
   double *Vlololo_d;
   double *Whihihi_h = new double[nrows*szt]; // W on the host
   double *Wlohihi_h = new double[nrows*szt];
   double *Whilohi_h = new double[nrows*szt];
   double *Wlolohi_h = new double[nrows*szt];
   double *Whihilo_h = new double[nrows*szt];
   double *Wlohilo_h = new double[nrows*szt];
   double *Whilolo_h = new double[nrows*szt];
   double *Wlololo_h = new double[nrows*szt];
   double *Whihihi_d;                         // W on the device
   double *Wlohihi_d;
   double *Whilohi_d;
   double *Wlolohi_d;
   double *Whihilo_d;
   double *Wlohilo_d;
   double *Whilolo_d;
   double *Wlololo_d;
   double *WYThihihi_h = new double[nrows2];  // W*Y^T 
   double *WYTlohihi_h = new double[nrows2];
   double *WYThilohi_h = new double[nrows2];
   double *WYTlolohi_h = new double[nrows2];
   double *WYThihilo_h = new double[nrows2];
   double *WYTlohilo_h = new double[nrows2];
   double *WYThilolo_h = new double[nrows2];
   double *WYTlololo_h = new double[nrows2];
   double *WYThihihi_d;                       // WYT on the device
   double *WYTlohihi_d;
   double *WYThilohi_d;
   double *WYTlolohi_d;
   double *WYThihilo_d;
   double *WYTlohilo_d;
   double *WYThilolo_d;
   double *WYTlololo_d;
   double *YWThihihi_h = new double[nrows2];  // Y*W^T
   double *YWTlohihi_h = new double[nrows2];
   double *YWThilohi_h = new double[nrows2];
   double *YWTlolohi_h = new double[nrows2];
   double *YWThihilo_h = new double[nrows2];
   double *YWTlohilo_h = new double[nrows2];
   double *YWThilolo_h = new double[nrows2];
   double *YWTlololo_h = new double[nrows2];
   double *YWThihihi_d;                       // YWT on the device
   double *YWTlohihi_d;
   double *YWThilohi_d;
   double *YWTlolohi_d;
   double *YWThihilo_d;
   double *YWTlohilo_d;
   double *YWThilolo_d;
   double *YWTlololo_d;
   double *QWYThihihi_h = new double[nrows2]; // Q*WY^T
   double *QWYTlohihi_h = new double[nrows2];
   double *QWYThilohi_h = new double[nrows2];
   double *QWYTlolohi_h = new double[nrows2];
   double *QWYThihilo_h = new double[nrows2];
   double *QWYTlohilo_h = new double[nrows2];
   double *QWYThilolo_h = new double[nrows2];
   double *QWYTlololo_h = new double[nrows2];
   double *QWYThihihi_d;                      // QWYT on the device
   double *QWYTlohihi_d;
   double *QWYThilohi_d;
   double *QWYTlolohi_d;
   double *QWYThihilo_d;
   double *QWYTlohilo_d;
   double *QWYThilolo_d;
   double *QWYTlololo_d;
   double *YWTChihihi_h = new double[dim];    // YWT*C on the host
   double *YWTClohihi_h = new double[dim];
   double *YWTChilohi_h = new double[dim];
   double *YWTClolohi_h = new double[dim];
   double *YWTChihilo_h = new double[dim];
   double *YWTClohilo_h = new double[dim];
   double *YWTChilolo_h = new double[dim];
   double *YWTClololo_h = new double[dim];
   double *YWTChihihi_d;                      // YWTC on the device
   double *YWTClohihi_d;
   double *YWTChilohi_d;
   double *YWTClolohi_d;
   double *YWTChihilo_d;
   double *YWTClohilo_d;
   double *YWTChilolo_d;
   double *YWTClololo_d;
   double *RTdotvhihihi_h = new double[nrows2]; // R^T dotted with v
   double *RTdotvlohihi_h = new double[nrows2];
   double *RTdotvhilohi_h = new double[nrows2];
   double *RTdotvlolohi_h = new double[nrows2];
   double *RTdotvhihilo_h = new double[nrows2];
   double *RTdotvlohilo_h = new double[nrows2];
   double *RTdotvhilolo_h = new double[nrows2];
   double *RTdotvlololo_h = new double[nrows2];
   double *RTdotvhihihi_d;                      // RTdotv on the device
   double *RTdotvlohihi_d;
   double *RTdotvhilohi_d;
   double *RTdotvlolohi_d;
   double *RTdotvhihilo_d;                      // RTdotv on the device
   double *RTdotvlohilo_d;
   double *RTdotvhilolo_d;
   double *RTdotvlololo_d;
   double *bRTvhihihi_h = new double[nrows];  // beta*R^T*v
   double *bRTvlohihi_h = new double[nrows];
   double *bRTvhilohi_h = new double[nrows];
   double *bRTvlolohi_h = new double[nrows];
   double *bRTvhihilo_h = new double[nrows];
   double *bRTvlohilo_h = new double[nrows];
   double *bRTvhilolo_h = new double[nrows];
   double *bRTvlololo_h = new double[nrows];
   double *bRTvhihihi_d;                      // bRTv on the device
   double *bRTvlohihi_d;
   double *bRTvhilohi_d;
   double *bRTvlolohi_d;
   double *bRTvhihilo_d; 
   double *bRTvlohilo_d;
   double *bRTvhilolo_d;
   double *bRTvlololo_d;
   double *sumshihihi_h = new double[nrows];  // subsums for large house
   double *sumslohihi_h = new double[nrows];
   double *sumshilohi_h = new double[nrows];
   double *sumslolohi_h = new double[nrows];
   double *sumshihilo_h = new double[nrows];
   double *sumslohilo_h = new double[nrows];
   double *sumshilolo_h = new double[nrows];
   double *sumslololo_h = new double[nrows];
   double *sumshihihi_d;                      // sums on the device
   double *sumslohihi_d;
   double *sumshilohi_d;
   double *sumslolohi_d;
   double *sumshihilo_d;
   double *sumslohilo_d;
   double *sumshilolo_d;
   double *sumslololo_d;
   double sigmahihihi_h,sigmalohihi_h,sigmahilohi_h,sigmalolohi_h;
   double sigmahihilo_h,sigmalohilo_h,sigmahilolo_h,sigmalololo_h;
   double *sigmahihihi_d;                     // sigma on the device
   double *sigmalohihi_d;
   double *sigmahilohi_d;
   double *sigmalolohi_d;
   double *sigmahihilo_d;
   double *sigmalohilo_d;
   double *sigmahilolo_d;
   double *sigmalololo_d;

   int ix = 0;                          // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Ahihihi_h[ix]   = Ahihihi[i][j];
         Alohihi_h[ix]   = Alohihi[i][j];
         Ahilohi_h[ix]   = Ahilohi[i][j];
         Alolohi_h[ix]   = Alolohi[i][j];
         Ahihilo_h[ix]   = Ahihilo[i][j];
         Alohilo_h[ix]   = Alohilo[i][j];
         Ahilolo_h[ix]   = Ahilolo[i][j];
         Alololo_h[ix++] = Alololo[i][j];
      }

   ix = 0;                              // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qhihihi_h[ix]   = 1.0;
            Qlohihi_h[ix]   = 0.0;
            Qhilohi_h[ix]   = 0.0;
            Qlolohi_h[ix]   = 0.0;
            Qhihilo_h[ix]   = 1.0;
            Qlohilo_h[ix]   = 0.0;
            Qhilolo_h[ix]   = 0.0;
            Qlololo_h[ix++] = 0.0;
         }
         else
         {
            Qhihihi_h[ix]   = 0.0;
            Qlohihi_h[ix]   = 0.0;
            Qhilohi_h[ix]   = 0.0;
            Qlolohi_h[ix]   = 0.0;
            Qhihilo_h[ix]   = 0.0;
            Qlohilo_h[ix]   = 0.0;
            Qhilolo_h[ix]   = 0.0;
            Qlololo_h[ix++] = 0.0;
         }
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Ahihihi_d,sznum);
   cudaMalloc((void**)&Alohihi_d,sznum);
   cudaMalloc((void**)&Ahilohi_d,sznum);
   cudaMalloc((void**)&Alolohi_d,sznum);
   cudaMalloc((void**)&Ahihilo_d,sznum);
   cudaMalloc((void**)&Alohilo_d,sznum);
   cudaMalloc((void**)&Ahilolo_d,sznum);
   cudaMalloc((void**)&Alololo_d,sznum);
   cudaMemcpy(Ahihihi_d,Ahihihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alohihi_d,Alohihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahilohi_d,Ahilohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alolohi_d,Alolohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahihilo_d,Ahihilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alohilo_d,Alohilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahilolo_d,Ahilolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alololo_d,Alololo_h,sznum,cudaMemcpyHostToDevice);

   const size_t szbeta = szt*sizeof(double);
   cudaMalloc((void**)&betahihihi_d,szbeta);
   cudaMalloc((void**)&betalohihi_d,szbeta);
   cudaMalloc((void**)&betahilohi_d,szbeta);
   cudaMalloc((void**)&betalolohi_d,szbeta);
   cudaMalloc((void**)&betahihilo_d,szbeta);
   cudaMalloc((void**)&betalohilo_d,szbeta);
   cudaMalloc((void**)&betahilolo_d,szbeta);
   cudaMalloc((void**)&betalololo_d,szbeta);

   for(int i=0; i<szt; i++)
   {
      betahihihi_h[i] = 0.0;
      betalohihi_h[i] = 0.0;
      betahilohi_h[i] = 0.0;
      betalolohi_h[i] = 0.0;
      betahihilo_h[i] = 0.0;
      betalohilo_h[i] = 0.0;
      betahilolo_h[i] = 0.0;
      betalololo_h[i] = 0.0;
   }
   cudaMemcpy(betahihihi_d,betahihihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohihi_d,betalohihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilohi_d,betahilohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalolohi_d,betalolohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahihilo_d,betahihilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohilo_d,betalohilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilolo_d,betahilolo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalololo_d,betalololo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);  // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vhihihi_d,szVandW + szpad); // pad only in allocation
   cudaMalloc((void**)&Vlohihi_d,szVandW + szpad);
   cudaMalloc((void**)&Vhilohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vlolohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vhihilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vlohilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vhilolo_d,szVandW + szpad);
   cudaMalloc((void**)&Vlololo_d,szVandW + szpad);

   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vhihihi_h[ix] = 0.0;
      Vlohihi_h[ix] = 0.0; 
      Vhilohi_h[ix] = 0.0; 
      Vlolohi_h[ix] = 0.0; 
      Vhihilo_h[ix] = 0.0;
      Vlohilo_h[ix] = 0.0; 
      Vhilolo_h[ix] = 0.0; 
      Vlololo_h[ix++] = 0.0; 
   }
   Vhihihi_h[--ix] = 1.0; // initialize last vector for square tiles

   cudaMemcpy(Vhihihi_d,Vhihihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlohihi_d,Vlohihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhilohi_d,Vhilohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlolohi_d,Vlolohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhihilo_d,Vhihilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlohilo_d,Vlohilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhilolo_d,Vhilolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlololo_d,Vlololo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Whihihi_d,szVandW + szpad); // pad only in allocation
   cudaMalloc((void**)&Wlohihi_d,szVandW + szpad); 
   cudaMalloc((void**)&Whilohi_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlolohi_d,szVandW + szpad); 
   cudaMalloc((void**)&Whihilo_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlohilo_d,szVandW + szpad); 
   cudaMalloc((void**)&Whilolo_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlololo_d,szVandW + szpad); 

   cudaMalloc((void**)&RTdotvhihihi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlohihi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhilohi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlolohi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhihilo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlohilo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhilolo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlololo_d,szVandW + szpad);
   cudaMalloc((void**)&bRTvhihihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlohihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhilohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlolohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhihilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlohilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhilolo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlololo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshihihi_d,szhouse);
   cudaMalloc((void**)&sumslohihi_d,szhouse);
   cudaMalloc((void**)&sumshilohi_d,szhouse);
   cudaMalloc((void**)&sumslolohi_d,szhouse);
   cudaMalloc((void**)&sumshihilo_d,szhouse);
   cudaMalloc((void**)&sumslohilo_d,szhouse);
   cudaMalloc((void**)&sumshilolo_d,szhouse);
   cudaMalloc((void**)&sumslololo_d,szhouse);
   cudaMalloc((void**)&sigmahihihi_d,sizeof(double));
   cudaMalloc((void**)&sigmalohihi_d,sizeof(double));
   cudaMalloc((void**)&sigmahilohi_d,sizeof(double));
   cudaMalloc((void**)&sigmalolohi_d,sizeof(double));
   cudaMalloc((void**)&sigmahihilo_d,sizeof(double));
   cudaMalloc((void**)&sigmalohilo_d,sizeof(double));
   cudaMalloc((void**)&sigmahilolo_d,sizeof(double));
   cudaMalloc((void**)&sigmalololo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYThihihi_d,szWYT + szpad); // pad for W*Y^T product
   cudaMalloc((void**)&WYTlohihi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThilohi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYTlolohi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThihilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTlohilo_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThilolo_d,szWYT + szpad); 
   cudaMalloc((void**)&WYTlololo_d,szWYT + szpad); 
   cudaMalloc((void**)&Qhihihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qlohihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qhilohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qlolohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qhihilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qlohilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qhilolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qlololo_d,szWYT + szpad);
   cudaMemcpy(Qhihihi_d,Qhihihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlohihi_d,Qlohihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhilohi_d,Qhilohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlolohi_d,Qlolohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhihilo_d,Qhihilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlohilo_d,Qlohilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhilolo_d,Qhilolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlololo_d,Qlololo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYThihihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlohihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThilohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlolohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThihilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlohilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThilolo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlololo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWThihihi_d,szYWT + szpad); // pad for Y*W^T product
   cudaMalloc((void**)&YWTlohihi_d,szYWT + szpad);
   cudaMalloc((void**)&YWThilohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlolohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWThihilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlohilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWThilolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlololo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTChihihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTClohihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTChilohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTClolohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTChihilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTClohilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTChilolo_d,sznum + szpad);
   cudaMalloc((void**)&YWTClololo_d,sznum + szpad);

   *houselapms = 0.0; *RTvlapms = 0.0; *tileRlapms = 0.0; *vb2Wlapms = 0.0;
   *WYTlapms = 0.0; *QWYTlapms = 0.0; *Qaddlapms = 0.0;
   *YWTlapms = 0.0; *YWTClapms = 0.0; *Raddlapms = 0.0;
   *addcnt = 0; *mulcnt = 0; *divcnt = 0; *sqrtcnt = 0;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      int colidx,nrows1;

      for(int L=0; L<szt; L++)  // L runs over the columns in one block
      {
         colidx = k*szt + L;              // index of the current column
         nrows1 = nrows - colidx - 1;     // #rows in Householder vector - 1

         if(verbose)
            cout << "-> current column : " << colidx << endl
                 << "-> #nrows in Householder vector - 1 : "
                 << nrows1 << endl;

         if(nrows1 <= szt)
         {
            GPU_dbl8_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                   Ahihihi_h,   Alohihi_h,   Ahilohi_h,   Alolohi_h,
                   Ahihilo_h,   Alohilo_h,   Ahilolo_h,   Alololo_h,
                   Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                   Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
                   vhihihi_h,   vlohihi_h,   vhilohi_h,   vlolohi_h,
                   vhihilo_h,   vlohilo_h,   vhilolo_h,   vlololo_h,
                   Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
                   Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
                betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
                betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
                betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            GPU_dbl8_small_leftRupdate
               (nrows,ncols,szt,colidx,k,L,
                   Ahihihi_h,   Alohihi_h,   Ahilohi_h,   Alolohi_h,
                   Ahihilo_h,   Alohilo_h,   Ahilolo_h,   Alololo_h,
                   Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                   Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
                   Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
                   Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
                betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
                betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
                betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                tileRlapms,addcnt,mulcnt,verbose);
         }
         else
         {
            GPU_dbl8_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                     Ahihihi_h,     Alohihi_h,     Ahilohi_h,     Alolohi_h,
                     Ahihilo_h,     Alohilo_h,     Ahilolo_h,     Alololo_h,
                     Ahihihi_d,     Alohihi_d,     Ahilohi_d,     Alolohi_d,
                     Ahihilo_d,     Alohilo_d,     Ahilolo_d,     Alololo_d,
                     vhihihi_h,     vlohihi_h,     vhilohi_h,     vlolohi_h,
                     vhihilo_h,     vlohilo_h,     vhilolo_h,     vlololo_h,
                     Vhihihi_d,     Vlohihi_d,     Vhilohi_d,     Vlolohi_d,
                     Vhihilo_d,     Vlohilo_d,     Vhilolo_d,     Vlololo_d,
                  betahihihi_h,  betalohihi_h,  betahilohi_h,  betalolohi_h,
                  betahihilo_h,  betalohilo_h,  betahilolo_h,  betalololo_h,
                  betahihihi_d,  betalohihi_d,  betahilohi_d,  betalolohi_d,
                  betahihilo_d,  betalohilo_d,  betahilolo_d,  betalololo_d,
                  sumshihihi_h,  sumslohihi_h,  sumshilohi_h,  sumslolohi_h,
                  sumshihilo_h,  sumslohilo_h,  sumshilolo_h,  sumslololo_h,
                  sumshihihi_d,  sumslohihi_d,  sumshilohi_d,  sumslolohi_d,
                  sumshihilo_d,  sumslohilo_d,  sumshilolo_d,  sumslololo_d,
                &sigmahihihi_h,&sigmalohihi_h,&sigmahilohi_h,&sigmalolohi_h,
                &sigmahihilo_h,&sigmalohilo_h,&sigmahilolo_h,&sigmalololo_h,
                 sigmahihihi_d, sigmalohihi_d, sigmahilohi_d, sigmalolohi_d,
                 sigmahihilo_d, sigmalohilo_d, sigmahilolo_d, sigmalololo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            GPU_dbl8_medium_leftRupdate
               (nrows,ncols,szt,colidx,k,L,
                     Ahihihi_h,     Alohihi_h,     Ahilohi_h,     Alolohi_h,
                     Ahihilo_h,     Alohilo_h,     Ahilolo_h,     Alololo_h,
                     Ahihihi_d,     Alohihi_d,     Ahilohi_d,     Alolohi_d,
                     Ahihilo_d,     Alohilo_d,     Ahilolo_d,     Alololo_d,
                     Vhihihi_d,     Vlohihi_d,     Vhilohi_d,     Vlolohi_d,
                     Vhihilo_d,     Vlohilo_d,     Vhilolo_d,     Vlololo_d,
                  betahihihi_h,  betalohihi_h,  betahilohi_h,  betalolohi_h,
                  betahihilo_h,  betalohilo_h,  betahilolo_h,  betalololo_h,
                  betahihihi_d,  betalohihi_d,  betahilohi_d,  betalolohi_d,
                  betahihilo_d,  betalohilo_d,  betahilolo_d,  betalololo_d,
                RTdotvhihihi_h,RTdotvlohihi_h,RTdotvhilohi_h,RTdotvlolohi_h,
                RTdotvhihilo_h,RTdotvlohilo_h,RTdotvhilolo_h,RTdotvlololo_h,
                RTdotvhihihi_d,RTdotvlohihi_d,RTdotvhilohi_d,RTdotvlolohi_d,
                RTdotvhihilo_d,RTdotvlohilo_d,RTdotvhilolo_d,RTdotvlololo_d,
                  bRTvhihihi_h,  bRTvlohihi_h,  bRTvhilohi_h,  bRTvlolohi_h,
                  bRTvhihilo_h,  bRTvlohilo_h,  bRTvhilolo_h,  bRTvlololo_h,
                  bRTvhihihi_d,  bRTvlohihi_d,  bRTvhilohi_d,  bRTvlolohi_d,
                  bRTvhihilo_d,  bRTvlohilo_d,  bRTvhilolo_d,  bRTvlololo_d,
                RTvlapms,tileRlapms,addcnt,mulcnt,verbose);
         }
      }
      GPU_dbl8_medium_VB_to_W
         (nrows,szt,szt,k,
             Vhihihi_h,   Vlohihi_h,   Vhilohi_h,   Vlolohi_h,
             Vhihilo_h,   Vlohilo_h,   Vhilolo_h,   Vlololo_h,
             Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
             Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
             Whihihi_h,   Wlohihi_h,   Whilohi_h,   Wlolohi_h,
             Whihilo_h,   Wlohilo_h,   Whilolo_h,   Wlololo_h,
             Whihihi_d,   Wlohihi_d,   Whilohi_d,   Wlolohi_d,
             Whihilo_d,   Wlohilo_d,   Whilolo_d,   Wlololo_d,
           WYThihihi_h, WYTlohihi_h, WYThilohi_h, WYTlolohi_h,
           WYThihilo_h, WYTlohilo_h, WYThilolo_h, WYTlololo_h,
           WYThihihi_d, WYTlohihi_d, WYThilohi_d, WYTlolohi_d,
           WYThihilo_d, WYTlohilo_d, WYThilolo_d, WYTlololo_d,
          betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
          betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
          betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
          betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);
/*
      GPU_dbl2_small_WYT
         (nrows-k*szt,szt,Whi_d,Wlo_d,Vhi_d,Vlo_d,WYThi_d,WYTlo_d,
          WYThi_h,WYTlo_h,WYTlapms,verbose);
 */
      GPU_dbl8_small_QWYT
         (nrows,szt,k,
             Qhihihi_d,   Qlohihi_d,   Qhilohi_d,   Qlolohi_d,
             Qhihilo_d,   Qlohilo_d,   Qhilolo_d,   Qlololo_d,
           WYThihihi_d, WYTlohihi_d, WYThilohi_d, WYTlolohi_d,
           WYThihilo_d, WYTlohilo_d, WYThilolo_d, WYTlololo_d,
          QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
          QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d,
          QWYThihihi_h,QWYTlohihi_h,QWYThilohi_h,QWYTlolohi_h,
          QWYThihilo_h,QWYTlohilo_h,QWYThilolo_h,QWYTlololo_h,
             Qhihihi_h,   Qlohihi_h,   Qhilohi_h,   Qlolohi_h,
             Qhihilo_h,   Qlohilo_h,   Qhilolo_h,   Qlololo_h,
          QWYTlapms,addcnt,mulcnt,verbose);

      GPU_dbl8_small_Qupdate
         (nrows,szt,k,
             Qhihihi_d,   Qlohihi_d,   Qhilohi_d,   Qlolohi_d,
             Qhihilo_d,   Qlohilo_d,   Qhilolo_d,   Qlololo_d,
          QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
          QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d,
             Qhihihi_h,   Qlohihi_h,   Qhilohi_h,   Qlolohi_h,
             Qhihilo_h,   Qlohilo_h,   Qhilolo_h,   Qlololo_h,
          Qaddlapms,addcnt,verbose);

      if(k < nbt-1)                                           // update R
      {
         GPU_dbl8_small_YWT
            (nrows,szt,k,
               Vhihihi_d,  Vlohihi_d,  Vhilohi_d,  Vlolohi_d,
               Vhihilo_d,  Vlohilo_d,  Vhilolo_d,  Vlololo_d,
               Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
               Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
             YWThihihi_d,YWTlohihi_d,YWThilohi_d,YWTlolohi_d,
             YWThihilo_d,YWTlohilo_d,YWThilolo_d,YWTlololo_d,
             YWThihihi_h,YWTlohihi_h,YWThilohi_h,YWTlolohi_h,
             YWThihilo_h,YWTlohilo_h,YWThilolo_h,YWTlololo_h,
             YWTlapms,addcnt,mulcnt,verbose);

         GPU_dbl8_small_YWTC
            (nrows,ncols,szt,k,
              YWThihihi_d, YWTlohihi_d, YWThilohi_d, YWTlolohi_d,
              YWThihilo_d, YWTlohilo_d, YWThilolo_d, YWTlololo_d,
                Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
             YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
             YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d,
             YWTChihihi_h,YWTClohihi_h,YWTChilohi_h,YWTClolohi_h,
             YWTChihilo_h,YWTClohilo_h,YWTChilolo_h,YWTClololo_h,
             YWTClapms,addcnt,mulcnt,verbose);

         GPU_dbl8_small_R_add_YWTC
            (nrows,ncols,szt,k,
                Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
             YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
             YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d,
                Ahihihi_h,   Alohihi_h,   Ahilohi_h,   Alolohi_h,
                Ahihilo_h,   Alohilo_h,   Ahilolo_h,   Alololo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qhihihi_h,Qhihihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlohihi_h,Qlohihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhilohi_h,Qhilohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlolohi_h,Qlolohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhihilo_h,Qhihilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlohilo_h,Qlohilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhilolo_h,Qhilolo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlololo_h,Qlololo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qhihihi[i][j] = Qhihihi_h[ix];
         Qlohihi[i][j] = Qlohihi_h[ix];
         Qhilohi[i][j] = Qhilohi_h[ix];
         Qlolohi[i][j] = Qlolohi_h[ix];
         Qhihilo[i][j] = Qhihilo_h[ix];
         Qlohilo[i][j] = Qlohilo_h[ix];
         Qhilolo[i][j] = Qhilolo_h[ix];
         Qlololo[i][j] = Qlololo_h[ix++];
      }

   cudaMemcpy(Ahihihi_h,Ahihihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alohihi_h,Alohihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahilohi_h,Ahilohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alolohi_h,Alolohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahihilo_h,Ahihilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alohilo_h,Alohilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahilolo_h,Ahilolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alololo_h,Alololo_d,sznum,cudaMemcpyDeviceToHost);

   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rhihihi[i][j] = Ahihihi_h[j*nrows+i];
         Rlohihi[i][j] = Alohihi_h[j*nrows+i];
         Rhilohi[i][j] = Ahilohi_h[j*nrows+i];
         Rlolohi[i][j] = Alolohi_h[j*nrows+i];
         Rhihilo[i][j] = Ahihilo_h[j*nrows+i];
         Rlohilo[i][j] = Alohilo_h[j*nrows+i];
         Rhilolo[i][j] = Ahilolo_h[j*nrows+i];
         Rlololo[i][j] = Alololo_h[j*nrows+i];
      }

   free(Ahihihi_h); free(Alohihi_h); free(Ahilohi_h); free(Alolohi_h);
   free(Ahihilo_h); free(Alohilo_h); free(Ahilolo_h); free(Alololo_h);
   free(Qhihihi_h); free(Qlohihi_h); free(Qhilohi_h); free(Qlolohi_h); 
   free(Qhihilo_h); free(Qlohilo_h); free(Qhilolo_h); free(Qlololo_h); 
   free(vhihihi_h); free(vlohihi_h); free(vhilohi_h); free(vlolohi_h);
   free(vhihilo_h); free(vlohilo_h); free(vhilolo_h); free(vlololo_h);
   free(Vhihihi_h); free(Vlohihi_h); free(Vhilohi_h); free(Vlolohi_h);
   free(Vhihilo_h); free(Vlohilo_h); free(Vhilolo_h); free(Vlololo_h);
   free(Whihihi_h); free(Wlohihi_h); free(Whilohi_h); free(Wlolohi_h);
   free(Whihilo_h); free(Wlohilo_h); free(Whilolo_h); free(Wlololo_h);
   free(sumshihihi_h); free(sumslohihi_h);
   free(sumshilohi_h); free(sumslolohi_h);
   free(sumshihilo_h); free(sumslohilo_h);
   free(sumshilolo_h); free(sumslololo_h);

   free(RTdotvhihihi_h); free(RTdotvlohihi_h);
   free(RTdotvhilohi_h); free(RTdotvlolohi_h);
   free(RTdotvhihilo_h); free(RTdotvlohilo_h);
   free(RTdotvhilolo_h); free(RTdotvlololo_h);
   free(bRTvhihihi_h); free(bRTvlohihi_h);
   free(bRTvhilohi_h); free(bRTvlolohi_h);
   free(bRTvhihilo_h); free(bRTvlohilo_h);
   free(bRTvhilolo_h); free(bRTvlololo_h);
   free(WYThihihi_h); free(QWYThihihi_h);
   free(WYThilohi_h); free(QWYThilohi_h);
   free(WYThihilo_h); free(QWYThihilo_h);
   free(WYThilolo_h); free(QWYThilolo_h);
   free(YWThihihi_h); free(YWTChihihi_h);
   free(YWThilohi_h); free(YWTChilohi_h);
   free(YWThihilo_h); free(YWTChihilo_h);
   free(YWThilolo_h); free(YWTChilolo_h);
   free(WYTlohihi_h); free(QWYTlohihi_h);
   free(WYTlolohi_h); free(QWYTlolohi_h);
   free(WYTlohilo_h); free(QWYTlohilo_h);
   free(WYTlololo_h); free(QWYTlololo_h);
   free(YWTlohihi_h); free(YWTClohihi_h);
   free(YWTlolohi_h); free(YWTClolohi_h);
   free(YWTlohilo_h); free(YWTClohilo_h);
   free(YWTlololo_h); free(YWTClololo_h);
}

void GPU_cmplx8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *houselapms, double *RHvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYHlapms, double *QWYHlapms, double *Qaddlapms,
   double *YWHlapms, double *YWHClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
}
