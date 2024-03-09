// The file hexa_double_gpufun.cu defines the code for the functions
// specified in hexa_double_gpufun.h

#include "double_double_gpufun.h"
#include "quad_double_gpufun.h"
#include "octo_double_gpufun.h"
#include "hexa_double_gpufun.h"

/************************* renormalizations ************************/

__device__ __forceinline__ void hdg_renorm16
 ( double f0, double f1, double f2, double f3, double f4, double f5,
   double f6, double f7, double f8, double f9, double f10, double f11,
   double f12, double f13, double f14, double f15, double f16, double *pr, 
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
{
}

__device__ __forceinline__ void hdg_fast_renorm
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10, double x11,
   double x12, double x13, double x14, double x15, double x16,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
{
}

__device__ __forceinline__ void hdg_renorm_add1
 ( double x0, double x1, double x2, double x3, double x4, double x5,
   double x6, double x7, double x8, double x9, double x10, double x11,
   double x12, double x13, double x14, double x15, double y,
   double *r0, double *r1, double *r2, double *r3, double *r4, double *r5,
   double *r6, double *r7, double *r8, double *r9, double *r10, double *r11,
   double *r12, double *r13, double *r14, double *r15 )
{
}

/*********************** copy and abs *****************************/

__device__ __forceinline__ void hdg_copy 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *b_hihihihi, double *b_lohihihi,
   double *b_hilohihi, double *b_lolohihi,
   double *b_hihilohi, double *b_lohilohi,
   double *b_hilolohi, double *b_lololohi,
   double *b_hihihilo, double *b_lohihilo,
   double *b_hilohilo, double *b_lolohilo,
   double *b_hihilolo, double *b_lohilolo,
   double *b_hilololo, double *b_lolololo )
{
}

__device__ __forceinline__ void hdg_abs 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *b_hihihihi, double *b_lohihihi,
   double *b_hilohihi, double *b_lolohihi,
   double *b_hihilohi, double *b_lohilohi,
   double *b_hilolohi, double *b_lololohi,
   double *b_hihihilo, double *b_lohihilo,
   double *b_hilohilo, double *b_lolohilo,
   double *b_hihilolo, double *b_lohilolo,
   double *b_hilololo, double *b_lolololo )
{
}

/******************* addition and subtraction *********************/

__device__ __forceinline__ void hdg_add 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
}

__device__ __forceinline__ void hdg_inc
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo )
{
}

__device__ __forceinline__ void hdg_inc_d
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b )
{
}

__device__ __forceinline__ void hdg_dec
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo )
{
}

__device__ __forceinline__ void hdg_minus 
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo )
{
}

__device__ __forceinline__ void hdg_sub 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
}

/**************  multiplications and division ***********************/

__device__ __forceinline__ void hdg_mul_pwr2
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
}

__device__ __forceinline__ void hdg_mul
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
}

__device__ __forceinline__ void hdg_sqr
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
}

__device__ __forceinline__ void hdg_mul_hd_d
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
}

__device__ __forceinline__ void hdg_mlt_d
 ( double *a_hihihihi, double *a_lohihihi,
   double *a_hilohihi, double *a_lolohihi,
   double *a_hihilohi, double *a_lohilohi,
   double *a_hilolohi, double *a_lololohi,
   double *a_hihihilo, double *a_lohihilo,
   double *a_hilohilo, double *a_lolohilo,
   double *a_hihilolo, double *a_lohilolo,
   double *a_hilololo, double *a_lolololo,
   double b )
{
}

__device__ __forceinline__ void hdg_div 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double b_hihihihi, double b_lohihihi,
   double b_hilohihi, double b_lolohihi,
   double b_hihilohi, double b_lohilohi,
   double b_hilolohi, double b_lololohi,
   double b_hihihilo, double b_lohihilo,
   double b_hilohilo, double b_lolohilo,
   double b_hihilolo, double b_lohilolo,
   double b_hilololo, double b_lolololo,
   double *c_hihihihi, double *c_lohihihi,
   double *c_hilohihi, double *c_lolohihi,
   double *c_hihilohi, double *c_lohilohi,
   double *c_hilolohi, double *c_lololohi,
   double *c_hihihilo, double *c_lohihilo,
   double *c_hilohilo, double *c_lolohilo,
   double *c_hihilolo, double *c_lohilolo,
   double *c_hilololo, double *c_lolololo )
{
}

/************************* square root *********************************/

__device__ __forceinline__ void hdg_sqrt 
 ( double a_hihihihi, double a_lohihihi,
   double a_hilohihi, double a_lolohihi,
   double a_hihilohi, double a_lohilohi,
   double a_hilolohi, double a_lololohi,
   double a_hihihilo, double a_lohihilo,
   double a_hilohilo, double a_lolohilo,
   double a_hihilolo, double a_lohilolo,
   double a_hilololo, double a_lolololo,
   double *b_hihihihi, double *b_lohihihi,
   double *b_hilohihi, double *b_lolohihi,
   double *b_hihilohi, double *b_lohilohi,
   double *b_hilolohi, double *b_lololohi,
   double *b_hihihilo, double *b_lohihilo,
   double *b_hilohilo, double *b_lolohilo,
   double *b_hihilolo, double *b_lohilolo,
   double *b_hilololo, double *b_lolololo )
{
}
