//     //////////////////////////////////////////////////////////////////
//     //                                                              //
//     // RIAM-COMPACT Natural Terrain Version Solver : 3dtopo-huge.f  //
//     //                                                              //
//     // Large-eddy simulation of neutral airflow over topography     //
//     //                                                              //
//     // This program is coded by Takanori UCHIDA, 2015/05/14         //
//     //                                                              //
//     //////////////////////////////////////////////////////////////////
//

#include "include/riamc.h"

using namespace pm_lib;

int main(int argc, char *argv[]) {

        struct timeval T1,T4;
        double TIMELA;

	gettimeofday(&T1, NULL);

// RIAMC class インスタンス
        RIAMC riamc;

// 初期化
        if( !riamc.Initialize(argc, argv) )
        {
           fprintf(stdout, "\n\tSolver initialize error.\n\n");
           return -1;
        }

// 繰り返し計算処理
        if ( !riamc.MainLoop() )
        {
          fprintf(stdout, "\n\tSolver execution error.\n\n");
          return -1;
        }

// 出力処理
        if( !riamc.Post() )
        {
          printf("\n\tSolver post error.\n");
          return -1;
        }

	gettimeofday(&T4, NULL);
	TIMELA = ( T4.tv_sec - T1.tv_sec ) + ( T4.tv_usec - T1.tv_usec) * 1.0E-6 ;
        if( riamc.myRank == 0 ) {
	printf("Elapsed Time was %f seconds\n", TIMELA);
        }
//
	return 0;
}

