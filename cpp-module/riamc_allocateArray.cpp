/*
 * RIAM-COMPACT Natural Terrain Version Solver
 *
 *
 * Copyright (c) 2015 RIAM, Kyushu University.
 * All rights reserved.
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "../include/riamc.h"

// #################################################################
/**
 * @brief 配列のアロケート 
 */
bool RIAMC::allocateArray()
{

  //仮想セル=1のS3D配列のアロケート
  int vc=1;

  //  X
  if( (X = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //  Y
  if( (Y = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //  Z
  if( (Z = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //  P
  if( (P = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //  U
  if( (U = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //  V
  if( (V = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //  W
  if( (W = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // UU 
  if( (UU = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // VV 
  if( (VV = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // WW 
  if( (WW = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // SGS
  if( (SGS = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //  ERRP
  if( (ERRP = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //仮想セル=2のS3D配列のアロケート
  int vc2=2;
  // UD 
  if( (UD = paraMngr->AllocDoubleS3D( vc2 )) == NULL )
  {
    return false;
  }

  // VD 
  if( (VD = paraMngr->AllocDoubleS3D( vc2 )) == NULL )
  {
    return false;
  }

  // WD 
  if( (WD = paraMngr->AllocDoubleS3D( vc2 )) == NULL )
  {
    return false;
  }

  //仮想セル=1のS3D配列のアロケート
  // GX 
  if( (GX = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // EY 
  if( (EY = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // TX 
  if( (TX = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // TY 
  if( (TY = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // TZ 
  if( (TZ = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  //仮想セル=0のS3D配列のアロケート
  int vc0 = 0;
  // C1
  if( (C1 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // C2
  if( (C2 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // C3
  if( (C3 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // C5
  if( (C5 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // C6
  if( (C6 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // C7
  if( (C7 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // C8
  if( (C8 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // C9
  if( (C9 = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  //仮想セル=1のS3D配列のアロケート
  // D1
  if( (D1 = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // D2
  if( (D2 = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // D3
  if( (D3 = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // D5
  if( (D5 = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // D6
  if( (D6 = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // YJA
  if( (YJA = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // RHS
  if( (RHS = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // ワーク配列 仮想セル=0のSD3配列のアロケート
  // VAN
  if( (VAN = paraMngr->AllocDoubleS3D( vc0 )) == NULL )
  {
    return false;
  }

  // ワーク配列 仮想セル=1のSD3配列のアロケート
  // CU
  if( (CU = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // CV
  if( (CV = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // CW
  if( (CW = paraMngr->AllocDoubleS3D( vc )) == NULL )
  {
    return false;
  }

  // uvw
  if( (uvw = paraMngr->AllocDoubleV3DEx( vc )) == NULL )
  {
    return false;
  }

  return true;

}
