/*
 * RIAM-COMPACT Natural Terrain Version Solver
 *
 * Copyright (c) 2016 RIAM, Kyushu University.
 * All rights reserved.
 *
 * Copyright (c) 2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 */

#ifndef _RIAMC_F_FUNC_H_
#define _RIAMC_F_FUNC_H_

#ifdef _WIN32

#define init_                 INIT
#define eddy_viscosity_       EDDY_VISCOSITY
#define calc_intermediate_v_  CALC_INTERMEDIATE_V
#define calc_contravariant_v_ CALC_CONTRAVARIANT_V
#define sor_                  SOR
#define calc_new_pressure_    CALC_NEW_PRESSURE
#define calc_new_velocity_    CALC_NEW_VELOCITY
#define display_              DISPLAY
#define display2_             DISPLAY2

#endif // _WIN32

//RIAM_COMPACT;
extern "C" {
extern void init_(int *, int *, int *, int *, int *, int *, int *,
                  const double *XLEN, const double *YLEN, const double *ZLEN, 
                  double *, double *,
                  double* X, double* Y, double* Z,
                  double* U, double* V, double* W,
                  double* P,
                  double* UU, double* VV, double   * WW,
                  double* SGS, double* YJA,
                  double* GX, double* EY, double   * TX, double   * TY,
                  double* TZ, double* C1, double   * C2, double   * C3,
                  double* C5, double* C6, double   * C7, double   * C8,
                  double* C9, double* D1, double   * D2, double   * D3,
                  double* D5, double* D6,
                  const double* LORG , const int* GSIZE, const int*head,
                  const int* nID, int *ierr );

extern void eddy_viscosity_(int *, int *, int *, int *, int *, int *,
                  double *, double *, double *,
                  double* Z, double* U, double* V,
                  double* W,
                  double* GX, double* EY, double* TX, double* TY,
                  double* TZ,
                  double* YJA, double* SGS , double* VAN, 
                  const int* nID );

extern void calc_intermediate_v_(int *, int *, int *, int *, int *, int *,
                  double *, double *, double *,
                  double* U, double* V, double* W,
                  double* UD, double* VD, double* WD,
                  double* UU, double* VV, double* WW,
                  double* YJA, double* SGS,
                  double* GX, double* EY, double* TX, double* TY,
                  double* TZ, double* C1, double* C2, double* C3,
                  double* C5, double* C6, double* C7, double* C8,
                  double* C9, 
                  const int* nID, int *ierr );

extern void calc_contravariant_v_(int *, int *, int *, int *, int *, int *, 
                  double *, 
                  double* U, double* V, double* W, 
                  double* GX, double* EY, double* TX, double* TY, 
                  double* TZ, double* UU, double* VV, 
                  double* WW, 
                  double* YJA, double* RHS,
                  double* CU, double* CV, double* CW, const int* nID);

extern void sor_( int *, int *, int *, int *, int *, int *, int *, int *, int *,  
                  double *, double *, double *, 
                  double* P, 
                  double* C1, double* C2, double* C3, double* C5, 
                  double* C6, double* C7, double* C8, double* C9, 
                  double* RHS, double* ERRP,
                  const int* nID, const int* head);

extern void calc_new_pressure_( int *, int *, int *, int *, int *, int *, 
                  double *, 
                  double* P, double* Z, double* W, 
                  double* C1, double* C2, double* C3, double* C5, 
                  double* C6, double* C7, double* C8, double* C9,
                  const int* nID, const int* GNX, const int*GNY, const int* GNZ,
                  double* PP);


extern void calc_new_velocity_( int *, int *, int *, int *, int *, int *, 
                  int *, int *, int *, int *, int *, int *, 
                  double *, 
                  double* X, double* P, double* U, 
                  double* V, double* W, 
                  double* UD, double* VD, double* WD, 
                  double* UU, double* VV, double* WW, 
                  double* GX, double* EY, double* TX, double* TY, 
                  double* TZ, 
                  double* YJA, 
                  double* D1, double* D2, double* D3, double* D5, 
                  double* D6,
                  const int* nID, const int* head, const int* tail );

extern void display_( int *, int *, int *, int *, int *, int *, int *, int *, int *, 
                  int *, int *, 
                  double *U, double *V, double *W,
                  int *, int *, 
                  double *, double *, double *, int * );

extern void display2_( int *, int *, int *, int *, int *, int *, 
                  double *, double *, 
                  double* U, double* V, double* W, 
                  double* P, 
                  double* UU, double* VV, double* WW,
                  int * );

extern void copy_s2vex_(int *MX2, int *MY2, int *MZ2, int *NX, int *NY, int *NZ,
                  double *S1, double *S2, double *S3, double *V3DEX);

}

#endif // _RIAMC_F_FUNC_H_
