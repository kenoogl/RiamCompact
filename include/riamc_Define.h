#ifndef _RIAMC_DEFINE_H_
#define _RIAMC_DEFINE_H_

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

/**
 * @file   riamc_Define.h
 * @brief  RIAMC Definition Header
 */

#include <float.h>
#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#define FB_VERS "1.0.0"

#define NOFACE 6

// PMlibの登録ラベル個数
#define PM_NUM_MAX 200

// ラベルの最大文字数
#define TM_LABEL_MAX 24

/** 並列化モード */
enum Parallel_mode
{
  Serial=1,
  OpenMP,
  FlatMPI,
  Hybrid
};

#define stamped_printf printf("%s (%d):  ",__FILE__, __LINE__), printf
#define stamped_fprintf fprintf(fp, "%s (%d):  ",__FILE__, __LINE__), fprintf


#endif // _RIAMC_DEFINE_H_
