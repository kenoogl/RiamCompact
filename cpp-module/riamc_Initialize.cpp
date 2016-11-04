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

int RIAMC::Initialize(int argc, char **argv)
{

  int ret = 0;  //リターンコード
  int vc = 1;   //仮想セル

  // 並列管理クラスのインスタンスと初期化
  paraMngr = cpm_ParaManager::get_instance(argc,argv);
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
#else
  num_threads = 1
#endif
  MPI_Comm_size(MPI_COMM_WORLD, &num_process); 

  double flop_count    = 0.0;  ///< flops計算用
  
  // パラメータローダのインスタンス生成
  TextParser tp;
  
  // パラメータのロードと保持
  if( argc>1 ) 
  {
    // 入力ファイルの指定
    string input_file = argv[1];
    int ierror=0;
    
    if ( (ierror = tp.read(input_file)) != TP_NO_ERROR )
    {
      stamped_printf("\tError at reading '%s' file : %d\n", input_file.c_str(), ierror);
      exit(0);
    }

    // 固定パラメータ
    if ( !setParameters(&tp) ) return 0;

  } else {
    stamped_printf("\tError undefined input tp file\n");
    return 0;
  }
  
  // タイミング測定の初期化
  PM.initialize( PM_NUM_MAX );
  myRank = paraMngr->GetMyRankID();
  PM.setRankInfo( myRank );
  PM.setParallelMode(setParallelism(), num_threads, num_process);
  set_timing_label();
  
  // パラメータ表示
  if( myRank == 0 ) {
  printParameters(stdout);
  }
  
  // タイミング測定開始
  TIMING_start("Initialization_Section");

  // 領域分割
  if( (ret = paraMngr->NodeInit(G_division, size, G_origin, G_region)) 
             !=CPM_SUCCESS )
  {
      stamped_printf("\tError at NodeInit errocode %d\n",ret);
      return 0;
  }

  //自ランクのHeadIndexの取得
  const int* head = paraMngr->GetNodeHeadIndex();

  //自ランクのTailIndexの取得
  const int* tail = paraMngr->GetNodeTailIndex();

#ifdef _CDM_OUTPUT
  //CDM initialize
  //CDM用headとtailの作成
  int cdm_head[3] = {head[0]+1, head[1]+1, head[2]+1};
  int cdm_tail[3] = {tail[0]+1, tail[1]+1, tail[2]+1};
 
  //Pressure
  dfi_p = cdm_DFI::WriteInit( MPI_COMM_WORLD
                            , cdm_DFI::Generate_DFI_Name("prs")
                            , "sph"
                            , "prs"
                            , CDM::E_CDM_FMT_SPH
                            , vc
                            , CDM::E_CDM_FLOAT64
                            , 1
                            , "proc.dfi"
                            , paraMngr->GetGlobalNodeSize()
                            , paraMngr->GetPitch()
                            , paraMngr->GetGlobalOrigin()
                            , paraMngr->GetDivNum()
                            , cdm_head
                            , cdm_tail
                            , "hogehoge"
                            , CDM::E_CDM_OFF ); 

  if( dfi_p )
  {
    dfi_p->setVariableName(0, "Pressure");
    dfi_p->WriteProcDfiFile(MPI_COMM_WORLD, false, 0, 0);
  }

  // Velocity
  dfi_v = cdm_DFI::WriteInit( MPI_COMM_WORLD
                            , cdm_DFI::Generate_DFI_Name("vel")
                            , "sph"
                            , "vel"
                            , CDM::E_CDM_FMT_SPH
                            , vc
                            , CDM::E_CDM_FLOAT64
                            , 3
                            , "proc.dfi"
                            , paraMngr->GetGlobalNodeSize()
                            , paraMngr->GetPitch()
                            , paraMngr->GetGlobalOrigin()
                            , paraMngr->GetDivNum()
                            , cdm_head
                            , cdm_tail
                            , "hogehoge"
                            , CDM::E_CDM_OFF );

  if( dfi_v )
  {
    dfi_v->setVariableName(0, "U");
    dfi_v->setVariableName(1, "V");
    dfi_v->setVariableName(2, "W");
  }
#endif  

  // 配列のアロケート
  TIMING_start("Allocate_Arrays");
  if( !allocateArray() ) 
  {
      stamped_printf("\tError at Allocate_Arrays\n");
      return 0;
  }
  TIMING_stop("Allocate_Arrays");

  
  // 初期化
  TIMING_start("1-init");
  flop_count = 0.0;

  //自ランクの格子数取得
  const int* lsize = paraMngr->GetLocalNodeSize();
  int MX2 = lsize[0];
  int MY2 = lsize[1];
  int MZ2 = lsize[2];
  int NX  = lsize[0];
  int NY  = lsize[1];
  int NZ  = lsize[2];


  //自ランクの隣接ランク番号を取得
  const int* nID = paraMngr->GetNeighborRankID();

  //全体の計算空間の領域長さ取得
  const double* GLen = paraMngr->GetGlobalRegion();
  //全体の計算空間の基点取得
  const double* GOrg = paraMngr->GetGlobalOrigin();
  //全体の計算空間の格子数取得
  const int* Gsize = paraMngr->GetGlobalNodeSize();

  int ierr=0;
  init_( &MX2, &MY2, &MZ2, &NX, &NY, &NZ, &ICON,
         &GLen[0], &GLen[1], &GLen[2], &TIME, &RE,
         X,   Y,   Z,   U,  V,  W,  P,    UU,   VV,   WW,   SGS,  YJA,
         GX,  EY,  TX,  TY, TZ, C1, C2,   C3,   C5,   C6,   C7,   C8,
         C9,  D1,  D2,  D3, D5, D6,
         GOrg, Gsize, head, nID, &ierr);
  TIMING_stop("1-init", flop_count);
  if( ierr != 0 ) return 0;


  // D1の袖通信を行う
  if( (paraMngr->BndCommS3D( D1, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at D1 Comm Band Cell \n");
      return 0;
  }

  // D2の袖通信を行う
  if( (paraMngr->BndCommS3D( D2, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at D2 Comm Band Cell \n");
      return 0;
  }

  // D3の袖通信を行う
  if( (paraMngr->BndCommS3D( D3, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at D3 Comm Band Cell \n");
      return 0;
  }

  // D5の袖通信を行う
  if( (paraMngr->BndCommS3D( D5, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at D5 Comm Band Cell \n");
      return 0;
  }

  // D6の袖通信を行う
  if( (paraMngr->BndCommS3D( D6, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at D6 Comm Band Cell \n");
      return 0;
  }

  // YJAの袖通信を行う
  if( (paraMngr->BndCommS3D( YJA, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at YJA Comm Band Cell \n");
      return 0;
  }
 
  // GXの袖通信を行う
  if( (paraMngr->BndCommS3D( GX, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at GX Comm Band Cell \n");
      return 0;
  }
 
  // EYの袖通信を行う
  if( (paraMngr->BndCommS3D( EY, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at EY Comm Band Cell \n");
      return 0;
  }
 
  // TXの袖通信を行う
  if( (paraMngr->BndCommS3D( TX, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at TX Comm Band Cell \n");
      return 0;
  }
 
  // TYの袖通信を行う
  if( (paraMngr->BndCommS3D( TY, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at TY Comm Band Cell \n");
      return 0;
  }
 
  // TZの袖通信を行う
  if( (paraMngr->BndCommS3D( TZ, MX2, MY2, MZ2, vc, vc, 0)) != CPM_SUCCESS )
  {
      stamped_printf("\tError at TZ Comm Band Cell \n");
      return 0;
  }

  TIMING_stop("Initialization_Section");

  return 1;
}



// #################################################################
/**
 * @brief パラメータの表示
 * @param [in] fp   FILE pointer
 */
void RIAMC::printParameters(FILE* fp)
{
  fprintf(fp, "\n------------\n");
  fprintf(fp, "PARAMETERS\n");
  fprintf(fp, "------------\n\n");
  fprintf(fp, "Domain\n");
  fprintf(fp, "\tGlobalOrigin   : (%16.6e, %16.6e, %16.6e)\n", G_origin[0], G_origin[1], G_origin[2]);
  fprintf(fp, "\tGlobalGrid     : (%16d, %16d, %16d)\n", G_size[0], G_size[1], G_size[2]);
  fprintf(fp, "\tGlobalRegion   : (%16.6e, %16.6e, %16.6e)\n", G_region[0], G_region[1], G_region[2]);
  fprintf(fp, "\tGlobalPitch    : (%16.6e, %16.6e, %16.6e)\n", pitch[0], pitch[1], pitch[2]);
  fprintf(fp, "\tGlobalDivision : (%16d, %16d, %16d)\n", G_division[0], G_division[1], G_division[2]);
  fprintf(fp, "\n");
  fprintf(fp, "Parameter\n");
  fprintf(fp, "\tBOX\n");
  fprintf(fp, "\t\tMIN          : (%16d, %16d, %16d)\n",Box_min[0],Box_min[1],Box_min[2]);
  fprintf(fp, "\t\tMAX          : (%16d, %16d, %16d)\n",Box_max[0],Box_max[1],Box_max[2]);
  fprintf(fp, "\n");
  fprintf(fp, "\tReynolds       : %16.6e\n", RE);
  fprintf(fp, "\tAlpha          : %10.4f\n", ALPHA);
  fprintf(fp, "\tIteration Max  : %10d\n", ITER);
  fprintf(fp, "\tOmega          : %16.6e\n", OMEGA);
  fprintf(fp, "\tEpsilon        : %16.6e\n", EPS);
  fprintf(fp, "\tSmagorinSky Constant     : %16.6e\n", CS);
  fprintf(fp, "\n");
  fprintf(fp, "Time Control\n");
  fprintf(fp, "\tRestart        : %s\n", (ICON==0)?"Initial":"Restart");
  fprintf(fp, "\tLast step      : %16d\n", NSTEP);
  fprintf(fp, "\tDelta T        : %16.6e\n", DT);
  fprintf(fp, "\tInterval for Prs History : %10d\n", NSOR);
  fprintf(fp, "\tInterval for Log         : %10d\n", IPROG1);
  fprintf(fp, "\tInterval for Vis         : %10d\n", IPROG2);
  fprintf(fp, "\n\n");
}



// #################################################################
/**
 * @brief パラメータのロード
 * @param [in] tp      TextParser
 * @note 無次元パラメータ
 */
bool RIAMC::setParameters(TextParser* tp)
{
  string label, str;
  double ct;

  // 全計算領域の基点
  label = "/DomainInfo/GlobalOrigin";
  if ( !tp->getInspectedVector(label, G_origin, 3) )
  {
     G_origin[0]=0.0e0;
     G_origin[1]=0.0e0;
     G_origin[2]=0.0e0;
  }

  // 全計算格子数
  bool bsize = true;
  label = "/DomainInfo/GlobalGrid";
  if ( !tp->getInspectedVector(label, G_size, 3) )
  {
     bsize = false;
  } else {
     size[0] = G_size[0];
     size[1] = G_size[1];
     size[2] = G_size[2];
  }

  // 全計算領域の大きさ
  bool bregn = true;
  label = "/DomainInfo/GlobalRegion";
  if ( !tp->getInspectedVector(label, G_region, 3) )
  {
     bregn = false;
  }
 
  // 格子分割幅
  bool bpit = true;
  label = "/DomainInfo/GlobalPitch";
  if ( !tp->getInspectedVector(label, pitch, 3) )
  {
     bpit = false;
  }

  // 領域分割数
  label = "/DomainInfo/GlobalDivision";
  if ( !tp->getInspectedVector(label, G_division, 3) )
  {
     G_division[0] = 0;
     G_division[1] = 0;
     G_division[2] = 0;
  }
  
  // GlobalGrid の自動計算処理
  if( !bsize )
  {
    if( !bregn || !bpit ) 
    {
      fprintf(stdout, "tParsing error : fail to get GlobalRegion or GlobalPitch\n");
      return false;
    }
    size[0]=G_size[0]=G_region[0]/pitch[0]+1;
    size[1]=G_size[1]=G_region[1]/pitch[1]+1;
    size[2]=G_size[2]=G_region[2]/pitch[2]+1;
  }

  // GlobalRegion の自動計算
  if( !bregn )
  {
    if( !bsize || !bpit )
    {
      fprintf(stdout, "tParsing error : fail to get GlobalGrid or GlobalPitch\n");
      return false;
    }
    G_region[0] = (G_size[0]-1)*pitch[0];
    G_region[1] = (G_size[1]-1)*pitch[1];
    G_region[2] = (G_size[2]-1)*pitch[2];
  }

  //GlobalPitch の自動計算
  if( !bpit )
  {
    if( !bsize || !bregn ) 
    {
      fprintf(stdout, "tParsing error : fail to get GlobalGrid or GlobalRegion\n");
      return false;
    } 
    pitch[0] = G_region[0]/(G_size[0]-1);
    pitch[1] = G_region[1]/(G_size[1]-1);
    pitch[2] = G_region[2]/(G_size[2]-1);
  }

  MX = size[0];
  MY = size[1];
  MZ = size[2];

  XLEN = G_region[0];
  YLEN = G_region[1];
  ZLEN = G_region[2];

  // 立方体位置(始点)
  label = "/Parameter/BOX/MIN";
  if ( !(tp->getInspectedVector(label,Box_min,3 )) )
  {
    Box_min[0] = 101;
    Box_min[1] = 101;
    Box_min[2] = 1;

    Box_max[0] = 151;
    Box_max[1] = 151;
    Box_max[2] = 51;
  }

  IA = Box_min[0];
  JA = Box_min[1];
  KA = Box_min[2];


  // 立方体位置(終点)
  label = "/Parameter/BOX/MAX";
  if ( !(tp->getInspectedVector(label,Box_max,3 )) )
  {
    Box_min[0] = 101;
    Box_min[1] = 101;
    Box_min[2] = 1;

    Box_max[0] = 151;
    Box_max[0] = 151;
    Box_max[0] = 51;
  }

  IB = Box_max[0];
  JB = Box_max[1];
  KB = Box_max[2];

  // 物理パラメータ
  label = "/Parameter/Reynolds";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    RE = 1.0e+4;
  } else {
    RE = (double)ct;
  }
  
  // 風上化の数値粘性係数
  label = "/Parameter/Alpha";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    ALPHA = 0.5;
  } else {
    ALPHA = (double)ct;
  }
  
  // 最大反復数
  label = "/Parameter/IterationMax";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    ITER = 100;
  } else {
    ITER = (int)ct;
  }
  
  // SORの加速係数
  label = "/Parameter/Omega";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    OMEGA = 0.9;
  } else {
    OMEGA = (double)ct;
  }

  // 収束閾値
  label = "/Parameter/Epsilon";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    EPS = 0.001;
  } else {
    EPS = (double)ct;
  }

  // SMAGORINSKY CONSTANT
  label = "/Parameter/SmagorinSky_Constant";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    CS = 0.1;
  } else {
    CS = (double)ct;
  }

 
  // リスタート指定
  label = "/TimeControl/Restart";
  if ( !(tp->getInspectedValue(label, str)) )
  {
    ICON = 0; 
  } else {
    if ( !strcasecmp(str.c_str(), "initial") )
    {
      ICON = 0;
    }
    else if ( !strcasecmp(str.c_str(), "restart") )
    {
      ICON = 1;
    }
    else
    {
      fprintf(stdout, "\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      return false;
    }
  }
  
  // 計算ステップ数
  label = "/TimeControl/LastStep";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    NSTEP = 50000;
  } else {
    NSTEP = (unsigned)ct;
  }
  
  // 時間積分幅
  label = "/TimeControl/DeltaT";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    DT = 0.001;
  } else {
    DT = (double)ct;
  }

  // 圧力履歴の出力インターバル
  label = "/TimeControl/Interval_Pressure_History";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    NSOR = 1;
  } else {
    NSOR = (int)ct; 
  }

  // 可視化データの出力インターバル
  label = "/TimeControl/Interval_Vis_Output";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    IPROG2 = 500;
  } else {
    IPROG2 = (int)ct;
  }
  
  // 可視化データの出力インターバル
  label = "/TimeControl/Interval_Log_Output";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    IPROG1 = 5000;
  } else {
    IPROG1 = (int)ct;
  }

  return true;

}
