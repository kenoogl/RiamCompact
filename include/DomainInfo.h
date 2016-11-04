#ifndef _FB_DOMAIN_INFO_H_
#define _FB_DOMAIN_INFO_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   DomainInfo.h
 * @brief  FlowBase DomainInfo class Header
 */

//#include "cpm_ParaManager.h"
#include <string>
#include "riamc_Define.h"


class DomainInfo {
  
public:
  //cpm_ParaManager* paraMngr; ///< Cartesian Partition Manager
  

  int procGrp;         ///< プロセスグループ番号
  int myRank;          ///< 自ノードのランク番号
  int numProc;         ///< 全ランク数
  
  int nID[6];          ///< 隣接ブロックのランク番号
  int head[3];         ///< 開始インデクス（グローバルインデクス, Fortran）
  
  int guide;           ///< ガイドセル数
  int G_division[3];   ///< プロセス分割数
  
  double pitch[3];     ///< 格子幅 (Non-dimensional)
  double pitchD[3];    ///< 格子幅 (有次元)
  
  int size[3];            ///< 領域分割数 (Local, Non-dimensional)
  double origin[3];    ///< 領域基点   (Local, Non-dimensional)
  double region[3];    ///< 領域サイズ (Local, Non-dimensional)
  double originD[3];   ///< 領域基点   (Local, 有次元)
  double regionD[3];   ///< 領域サイズ (Local, 有次元)
  
  int G_size[3];          ///< 領域分割数 (Global, Non-dimensional)
  double G_origin[3];  ///< 領域基点   (Global, Non-dimensional)
  double G_region[3];  ///< 領域サイズ (Global, Non-dimensional)
  double G_originD[3]; ///< 領域基点   (Global, 有次元)
  double G_regionD[3]; ///< 領域サイズ (Global, 有次元)

  
  /** コンストラクタ */
  DomainInfo() {
    procGrp = 0;
    myRank  = -1;
    numProc = 0;
    for (int i=0; i<NOFACE; i++) nID[i] = -1;
    
    for (int i=0; i<3; i++)
    {
      head[i]       = 0;
      size[i]       = 0;
      G_size[i]     = 0;
      G_division[i] = 0;
      pitch[i]      = 0.0;
      origin[i]     = 0.0;
      region[i]     = 0.0;
      G_origin[i]   = 0.0;
      G_region[i]   = 0.0;
      pitchD[i]     = 0.0;
      originD[i]    = 0.0;
      regionD[i]    = 0.0;
      G_originD[i]  = 0.0;
      G_regionD[i]  = 0.0;
    }
    
    guide = 1;
    //paraMngr = NULL;
  }
  
  /** デストラクタ */
  virtual ~DomainInfo() {}
  
};

#endif // _FB_DOMAIN_INFO_H_
