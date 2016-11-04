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
#include <string>

// #################################################################
// constructor
RIAMC::RIAMC()
{

  NSTEP = 50000;

  TIME = 0.0e0;

  G_origin[0] = 0.0;
  G_origin[1] = 0.0;
  G_origin[2] = 0.0;

  MX = 501;
  MY = 251;
  MZ = 201;

  size[0] = G_size[0] = MX;
  size[1] = G_size[1] = MY;
  size[2] = G_size[2] = MZ;

  XLEN = 10.0;
  YLEN = 5.0;
  ZLEN = 4.0;

  G_region[0] = XLEN;
  G_region[1] = YLEN;
  G_region[2] = ZLEN;

  ITER = 100;
  OMEGA = 0.9;
  EPS   = 0.001;
  ALPHA = 0.5;

  DT = 0.001;
  IPROG1 = 50000;
  IPROG2 = 500;
  NSOR = 1;

  IA = 101;
  IB = 151;
  JA = 101;
  JB = 151;
  KA = 1;
  KB = 51;

  Box_min[0] = IA;
  Box_min[1] = JA;
  Box_min[2] = KA;
  Box_max[0] = IB;
  Box_max[1] = JB;
  Box_max[2] = KB;

  RE = 10000.0;
  CS = 0.1;

  order_of_PM_key = 0;

  for( int i=0; i<3; i++)
  {
    Box_min[0]=0;
    Box_max[0]=0;
  }

  ICON = 0;

  num_threads = 0;
  Parallelism = 0;
  num_process = 0;
 
}



// #################################################################
/* @brief メモリ消費情報を表示
 * @param [in]     fp    ファイルポインタ
 * @param [in,out] G_mem グローバルメモリサイズ
 * @param [in]     L_mem ローカルメモリサイズ
 * @param [in]     str   表示用文字列
 */
void RIAMC::displayMemoryInfo(FILE* fp, double G_mem, double L_mem, const char* str)
{
/*
    FBUtility::MemoryRequirement(str, G_mem, L_mem, fp);
    fprintf(fp, "\n\n");
*/
}



// #################################################################
/**
 * @brief シミュレーションの1ステップの処理
 */
int RIAMC::MainLoop()
{
  double flop_count = 0.0;  ///< flops計算用

  //自ランクの格子数取得
  const int* lsize = paraMngr->GetLocalNodeSize();
  int MX2 = lsize[0];
  int MY2 = lsize[1];
  int MZ2 = lsize[2];
  int NX  = lsize[0];
  int NY  = lsize[1];
  int NZ  = lsize[2];

  //全体の計算空間での格子数取得
  const int* gsize = paraMngr->GetGlobalNodeSize();
  int NOXYZ = (gsize[0]-2) * (gsize[1]-2) * (gsize[2]-2);

  double REI = 1.0/RE;
  double CS2 = CS * CS;

  int ILAP;
  double RMSP;

  int ICOUNT = 0;
  int T_NSTEP = NSTEP;    //最終計算ステップをダミーにセット
  int vc = 1;             //仮想セル数
  
  // 自ランクの隣接ランク番号を取得
  const int* nID = paraMngr->GetNeighborRankID();

  // 自ランクのhead&tailを取得
  const int* head = paraMngr->GetNodeHeadIndex();
  const int* tail = paraMngr->GetNodeTailIndex();

  double PP;

  for (int ISTEP = 1; ISTEP <= NSTEP; ISTEP++)
  {

    TIME = TIME + DT;
    if( myRank == 0 ) {
    printf("time=%e\n", TIME);fflush(stdout);
    }

    TIMING_start("2-eddy_viscosity");
    flop_count = 0.0 ;
    eddy_viscosity_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &RE, &REI, &CS2, 
                    Z, U, V, W, GX, EY, TX, TY, TZ, YJA, SGS, VAN, nID);
    TIMING_stop("2-eddy_viscosity", flop_count);

    // SGSの袖通信を行う
    if( (paraMngr->BndCommS3D( SGS, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at SGS Comm Band Cell \n");
        return 0;
    }    

    TIMING_start("3-calc_intermediate_v");
    flop_count = 0.0 ;
    int ierr=0;
    calc_intermediate_v_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &ALPHA, &DT, &REI,
                         U, V, W, UD, VD, WD, UU, VV, WW, YJA, SGS, GX, EY, TX, TY,
                         TZ, C1, C2, C3, C5, C6, C7, C8, C9, nID, &ierr);
    TIMING_stop("3-calc_intermediate_v", flop_count);
    if( ierr != 0 ) return 0;
  
    // 計算結果 U,V,Wの袖通信を行う
    if( (paraMngr->BndCommS3D( U, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at U Comm Band Cell \n");
        return 0;
    }    
    if( (paraMngr->BndCommS3D( V, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at V Comm Band Cell \n");
        return 0;
    }    
    if( (paraMngr->BndCommS3D( W, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at W Comm Band Cell \n");
        return 0;
    }    
 
    TIMING_start("4-calc_contravariant_v");
    flop_count = 0.0;
    calc_contravariant_v_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &DT, U, V, W,
                          GX, EY, TX, TY, TZ, UU, VV, WW, YJA, RHS,
                          CU, CV, CW, nID);
    TIMING_stop("4-calc_contravariant_v", flop_count);

    TIMING_start("5-SOR");
    flop_count = 0.0;
    sor_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &ILAP, &ITER, &NOXYZ, &OMEGA, &EPS,
         &RMSP, P, C1, C2, C3, C5, C6, C7, C8, C9, RHS, ERRP, nID, head);
    TIMING_stop("5-SOR", flop_count);
 
    TIMING_start("6-calc_new_pressure");
    flop_count = 0.0;
    calc_new_pressure_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &REI, P, Z, W, C1, C2, 
                       C3, C5, C6, C7, C8, C9, nID, &gsize[0], &gsize[1], &gsize[2],
                       &PP);
    TIMING_stop("6-calc_new_pressure", flop_count);

    // 計算結果 Pの袖通信を行う
    if( (paraMngr->BndCommS3D( P, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at P Comm Band Cell \n");
        return 0;
    } 
    
    TIMING_start("7-calc_new_velocity");
    flop_count = 0.0;
    calc_new_velocity_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &KA, &KB, &JA, &JB, &IA, 
                       &IB, &DT, X, P, U, V, W, UD, VD, WD, UU, VV, WW, GX, EY, 
                       TX, TY, TZ, YJA, D1, D2, D3, D5, D6, 
                       nID, head, tail );
    TIMING_stop("7-calc_new_velocity", flop_count);

    // 計算結果 U,V,W,UU,VV,WWの袖通信を行う
    if( (paraMngr->BndCommS3D( U, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at U Comm Band Cell \n");
        return 0;
    } 
    if( (paraMngr->BndCommS3D( V, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at V Comm Band Cell \n");
        return 0;
    } 
    if( (paraMngr->BndCommS3D( W, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at W Comm Band Cell \n");
        return 0;
    } 
    if( (paraMngr->BndCommS3D( UU, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at UU Comm Band Cell \n");
        return 0;
    } 
    if( (paraMngr->BndCommS3D( VV, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at VV Comm Band Cell \n");
        return 0;
    } 
    if( (paraMngr->BndCommS3D( WW, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
    {
        stamped_printf("\tError at WW Comm Band Cell \n");
        return 0;
    } 
    
    TIMING_start("8-display_&_FILE_I/O");
    flop_count = 0.0;
    display_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &ISTEP, &T_NSTEP, &ILAP, &ICOUNT, 
             &NSOR, U, V, W, &IPROG1, &IPROG2, &TIME, &RE, &RMSP, &myRank);
    TIMING_stop("8-display_&_FILE_I/O", flop_count);

#ifdef _CDM_OUTPUT
    if( ISTEP%IPROG2 == 0 || ISTEP == NSTEP )
    {
       //// Pressure
      if( dfi_p )
      {
         // write
         dfi_p->WriteData( (unsigned)ISTEP
                         , TIME
                         , paraMngr->GetLocalNodeSize()
                         , 1
                         , vc
                         , P
                         , (double*)NULL
                         , false
                         , 0
                         , 0.0
                         );
      }
      //// Velocity
      if( dfi_v )
      {  
         // copy array
         copy_s2vex_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, U, V, W, uvw);
         // write
         dfi_v->WriteData( (unsigned)ISTEP
                         , TIME
                         , paraMngr->GetLocalNodeSize()
                         , 3
                         , vc
                         , uvw
                         , (double*)NULL
                         , false
                         , 0
                         , 0.0
                         );
      }
    }
#endif
    
  }
  
  return 1;
}


// #################################################################
/**
 * @brief シミュレーションの1ステップの処理
 */
int RIAMC::Post()
{

  double flop_count = 0.0;  ///< flops計算用

  const int* lsize = paraMngr->GetLocalNodeSize();

  int MX2 = lsize[0];
  int MY2 = lsize[1];
  int MZ2 = lsize[2];
  int NX  = lsize[0];
  int NY  = lsize[1];
  int NZ  = lsize[2];

  TIMING_start("8-display_&_FILE_I/O");
  flop_count = 0.0;
  display2_(&MX2, &MY2, &MZ2, &NX, &NY, &NZ, &TIME, &RE, U, V, W, P, UU, VV, WW,
            &myRank);
  TIMING_stop("8-display_&_FILE_I/O", flop_count);

  PM.gather();
  PM.print(stdout, "", "HPCS");
  PM.printDetail(stdout);

  return 1;

}
// #################################################################
/**
 * @brief  履歴の表示
 * @param [in] fp   FILE pointer
 * @param [in] step ステップ
 * @param [in] time 無次元時刻
 * @param [in] itr  反復回数
 * @param [in] rms  圧力反復のrms値
 */
void RIAMC::printHistory(FILE* fp, int step, double time, int itr, double rms)
{
  fprintf(fp, "%12d %16.6e %5d %16.6e\n", step, time, itr, rms);
}


// #################################################################
/* @brief 並列化と分割の方法を保持
 * @return 並列モード
 */
string RIAMC::setParallelism()
{
  string para_mode;
  
  num_threads = 1;
  
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
#endif
 
  num_process = 1;
   
  if ( num_threads > 1 )
  {
    Parallelism = OpenMP;
    para_mode = "OpenMP";
  }
  else
  {
    Parallelism = Serial;
    para_mode = "Serial";
  }
  
  return para_mode;
}



// #################################################################
/**
 * @brief タイミング測定区間にラベルを与えるラッパー
 * @param [in] label     ラベル
 * @param [in] type      測定対象タイプ(COMM or CALC)
 * @param [in] exclusive 排他測定フラグ(ディフォルトtrue)
 */
void RIAMC::set_label(const string label, PerfMonitor::Type type, bool exclusive)
{
  // 登録個数のチェック
  order_of_PM_key++;
  
  if ( order_of_PM_key > PM_NUM_MAX )
  {
    fprintf(stdout, "\tThe number of labels for Performance monitor goes over limit.\n");
    exit(0);
  }
  
  // 文字数がTM_LABEL_MAX-1を超えるものはカット
  if ( strlen(label.c_str()) > TM_LABEL_MAX-1 )
  {
    printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
  }
  
  // Performance Monitorへの登録
  PM.setProperties(label, type, exclusive);
}



// #################################################################
/**
 * @brief タイミング測定区間にラベルを与える
 */
void RIAMC::set_timing_label()
{
  // common
  set_label("Allocate_Arrays",         PerfMonitor::CALC);
  // common
  
  
  // Initialization_Section
  set_label("Initialization_Section",  PerfMonitor::CALC, false);
  set_label("1-init",                  PerfMonitor::CALC);
  // Initialization_Section
  
  
  // NS_Section
  set_label("2-eddy_viscosity",        PerfMonitor::CALC);
  set_label("3-calc_intermediate_v",   PerfMonitor::CALC);
  set_label("4-calc_contravariant_v",  PerfMonitor::CALC);
  set_label("5-SOR",                   PerfMonitor::CALC);
  set_label("6-calc_new_pressure",     PerfMonitor::CALC);
  set_label("7-calc_new_velocity",     PerfMonitor::CALC);
  set_label("8-display_&_FILE_I/O",    PerfMonitor::CALC);
  
}
