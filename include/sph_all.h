#pragma once
#include "density.h"
#include "force.h"
#include "init.h"
#include "integrator.h"
#include "io.h"
#include "kernel.h"
#include "constants.h"

#ifdef __CUDACC__
#include "cuda_utils.cuh"
#include "density.cuh"
#include "force.cuh"
#include "integrator.cuh"
#endif
