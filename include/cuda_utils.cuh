#pragma once
#include <float.h>

__device__ __forceinline__ double atomicMax_double(double* address, double val)
{
    // 將 double 指標轉型為 unsigned long long int 指標
    // 這是因為 atomicCAS 支援 64-bit 的 unsigned long long
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    // 讀取記憶體中目前的值
    unsigned long long int old = *address_as_ull;
    unsigned long long int assumed;

    do {
        assumed = old;
        // 將記憶體中的 bit 轉回 double 來進行大小比較
        double assumed_double = __longlong_as_double(assumed);
        // 如果傳入的 val 比較大，就準備寫入 val；否則保持原來的值
        double new_val = fmax(val, assumed_double);
        // 如果 new_val 根本沒有變大，代表不需要更新了，提早結束以節省效能
        if (new_val == assumed_double)break;
        // atomicCAS 核心邏輯：
        // 檢查 address_as_ull 裡面的值是否還是 assumed，如果是，就把它換成 new_val 的 bit 表示
        // 並回傳交換前的值到 old 中。
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(new_val));
    // 如果 old != assumed，代表在我們讀取到準備寫入的瞬間，別的 Thread 已經改了這塊記憶體，我們必須重試 (Spin-lock)
    } while (assumed != old);
    return __longlong_as_double(old);
}


__device__ __forceinline__ double atomicMin_double(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}
