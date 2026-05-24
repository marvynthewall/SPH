# include "sph_system.h"

int main(int argc, char *argv[]) {
    printf("====================================\n");
    printf("   SPH 模擬程式啟動中...             \n");
    printf("====================================\n");

    SPHSystem2D sph;
    int nx = 20;
    int ny = 30;
    int N = nx*ny;
    allocate_sph_system(&sph, N);
    // init_uniform_box(&sph, nx, ny);
    init_sod_2d(&sph, nx, ny);

    write_csv(&sph, "output_0000.csv");

    printf("\n初始化完成 ...\n");


    free_sph_system(&sph);
    return 0;
}
