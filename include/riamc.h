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

#ifndef _RIAMC_H_
#define _RIAMC_H_

#include <stdio.h>
#include <vector>

//#include "riamc_Version.h"
#include "riamc_Define.h"
#include "DomainInfo.h"
#include "riamc_Ffunc.h"

// Text parser
#include "TextParser.h"

// FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "fj_tool/fjcoll.h"
#endif

// Performance Monitor
#include "PerfMonitor.h"

// cpm_ParaManager
#include "cpm_ParaManager.h"

#ifdef _CDM_OUTPUT
// cdm_DFI
#include "cdm_DFI.h"
#endif


#ifndef _WIN32
#include <unistd.h>
#include <strings.h>
#else
#include "sph_win32_util.h"
#endif
#include <sys/types.h>

#if defined(IA32_LINUX) || defined(IA64_LINUX) || defined(SGI_ALTIX)
#include <sys/stat.h>
#endif

#ifdef MacOSX
#include <sys/uio.h>
#endif

using namespace std;
using namespace pm_lib;

class RIAMC : public DomainInfo {
  
private:

  unsigned NSTEP;               ///< セッションの終了ステップ数
  
  int MX;  ///< X方向格子数
  int MY;  ///< Y方向格子数
  int MZ;  ///< Z方向格子数

  double XLEN;    ///< X方向長さ
  double YLEN;    ///< Y方向長さ
  double ZLEN;    ///< Z方向長さ

  double RE;       ///< レイノルズ数
  double ALPHA;    ///< 数値粘性の係数
  double OMEGA;    ///< SORの加速係数
  double EPS;      ///< 収束閾値
  double CS;       ///< SMAGORINSKY CONSTANT
  int ITER;           ///< 最大反復回数
  
  int order_of_PM_key; ///< PMlib用の登録番号カウンタ < PM_NUM_MAX
 
  // 立方体の位置
  int Box_min[3];     ///< 立方体の始点
  int Box_max[3];     ///< 立方体の終点 
  int IA;             ///< 立方体のコーナー位置（Xの始点)
  int IB;             ///< 立方体のコーナー位置（Xの終点)
  int JA;             ///< 立方体のコーナー位置（Yの始点)
  int JB;             ///< 立方体のコーナー位置（Yの終点)
  int KA;             ///< 立方体のコーナー位置（Zの始点)
  int KB;             ///< 立方体のコーナー位置（Zの終点)


  //TimeControl
  int ICON;          ///< 0 - initial / 1 - restart
  double DT;         ///< 時間積分幅（無次元）
  int NSOR;             ///< 画面表示間隔
  int IPROG1;           ///< logファイル出力間隔
  int IPROG2;           ///< visファイル出力間隔
 
  int num_threads;      ///< スレッド数
  int Parallelism;      ///< 並列モード
  int num_process;      ///< プロセス数

  PerfMonitor PM;            ///< 性能モニタクラス

  cpm_ParaManager *paraMngr; ///< CPMlibクラス

#ifdef _CDM_OUTPUT
  cdm_DFI *dfi_p;       ///< Pressure 出力用 cpm_DFI
  cdm_DFI *dfi_v;       ///< Velocity 出力用 cpm_DFI
#endif

public:

  double TIME;

  double* X;             ///< double   X[MZ+2][MY+2][MX+2]
  double* Y;             ///< double   Y[MZ+2][MY+2][MX+2]
  double* Z;             ///< double   Z[MZ+2][MY+2][MX+2]
  double* P;             ///< double   P[MZ+2][MY+2][MX+2]
  double* U;             ///< double   U[MZ+2][MY+2][MX+2]
  double* V;             ///< double   V[MZ+2][MY+2][MX+2]
  double* W;             ///< double   W[MZ+2][MY+2][MX+2]
  double* UD;            ///< double  UD[MZ+2][MY+2][MX+2]
  double* VD;            ///< double  VD[MZ+2][MY+2][MX+2]
  double* WD;            ///< double  WD[MZ+2][MY+2][MX+2]
  double* UU;            ///< double  UU[MZ+2][MY+2][MX+2]
  double* VV;            ///< double  VV[MZ+2][MY+2][MX+2]
  double* WW;            ///< double  WW[MZ+2][MY+2][MX+2]
  
  double* GX;            ///< double  GX[MZ][MY][MX]
  double* EY;            ///< double  EY[MZ][MY][MX]
  double* TX;            ///< double  TX[MZ][MY][MX]
  double* TY;            ///< double  TY[MZ][MY][MX]
  double* TZ;            ///< double  TZ[MZ][MY][MX]
  double* C1;            ///< double  C1[MZ][MY][MX]
  double* C2;            ///< double  C2[MZ][MY][MX]
  double* C3;            ///< double  C3[MZ][MY][MX]
  double* C5;            ///< double  C5[MZ][MY][MX]
  double* C6;            ///< double  C6[MZ][MY][MX]
  double* C7;            ///< double  C7[MZ][MY][MX]
  double* C8;            ///< double  C8[MZ][MY][MX]
  double* C9;            ///< double  C9[MZ][MY][MX]

  double* D1;            ///< double  D1[MZ+2][MY+2][MX+2]
  double* D2;            ///< double  D2[MZ+2][MY+2][MX+2]
  double* D3;            ///< double  D3[MZ+2][MY+2][MX+2]
  double* D5;            ///< double  D5[MZ+2][MY+2][MX+2]
  double* D6;            ///< double  D6[MZ+2][MY+2][MX+2]
  double* SGS;           ///< double SGS[MZ+2][MY+2][MX+2]
  double* YJA;           ///< double YJA[MZ][MY][MX]
  double* RHS;           ///< double RHS[MZ][MY][MX]
  double* VAN;           ///< double VAN[MZ][MY][MX]

  //ワーク配列
  double* CU;            ///< double  UD[MZ+2][MY+2][MX+2]
  double* CV;            ///< double  UD[MZ+2][MY+2][MX+2]
  double* CW;            ///< double  UD[MZ+2][MY+2][MX+2]

  //SORでのエラー値格納配列
  double* ERRP;          ///< double   ERRP[MZ+2][MY+2][MX+2]

  //Vellocity出力用配列
  double* uvw;           ///< double   uvw[3][MZ+2][MY+2][MX+2]

public:
  /** コンストラクタ */
  RIAMC();
  
  /** デストラクタ */
  ~RIAMC() {};
  


  
  
private:
  
  /**
   * @brief タイミング測定開始
   * @param [in] key ラベル
   */
  inline void TIMING_start(const string key)
  {
    // PMlib Intrinsic profiler
    PM.start(key);
    
    const char* s_label = key.c_str();
    
    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( s_label );
#elif defined __FX_FAPP
    fapp_start( s_label, 0, 0);
#endif
  }
  
  
  /**
   * @brief タイミング測定終了
   * @param [in] key             ラベル
   * @param [in] flopPerTask    「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
   * @param [in] iterationCount  実行「タスク」数 (ディフォルト1)
   */
  inline void TIMING_stop(const string key, double flopPerTask=0.0, int iterationCount=1)
  {
    // Venus FX profiler
    const char* s_label = key.c_str();
    
#if defined __K_FPCOLL
    stop_collection( s_label );
#elif defined __FX_FAPP
    fapp_stop( s_label, 0, 0);
#endif
    
    // PMlib Intrinsic profiler
    PM.stop(key, flopPerTask, (unsigned)iterationCount);
  }

  
  // メモリ使用量の表示
  void displayMemoryInfo(FILE* fp, double G_mem, double L_mem, const char* str);
  
  
  // 履歴表示
  void printHistory(FILE* fp, int step, double time, int itr, double rms);
  
  
  // パラメータの表示
  void printParameters(FILE* fp);
  
  
  // 並列モード文字列の取得
  string setParallelism();
  
  
  // 時間積分幅や物理パラメータの設定
  bool setParameters(TextParser* tp);

  // 配列のアロケート
  bool allocateArray();
  
  
  // タイミング測定区間にラベルを与えるラッパー
  void set_label(const string label, PerfMonitor::Type type, bool exclusive=true);
  
  
  // タイミング測定区間にラベルを与える
  void set_timing_label();
  
  
  
public:
  
  // 初期化格子生成など
  int Initialize(int argc, char **argv);
  
  
  // シミュレーションの1ステップの処理
  int MainLoop();
  
  
  // シミュレーションの終了時の処理
  int Post();
  
};


#endif // _RIAMC_H_
