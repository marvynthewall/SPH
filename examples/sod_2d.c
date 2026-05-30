# include "sph_system.h"

int main(int argc, char *argv[]) {
    printf("====================================\n");
    printf("   SPH 模擬程式啟動中...             \n");
    printf("====================================\n");
    
    const char *output_filename = "output_0000.csv";
    double mass = 0.001;  // default mass

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_filename = argv[i + 1];
            i++; 
        }
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            mass = atof(argv[i + 1]); 
            i++; 
        }
    }

    SPHSystem2D sph;
    int x = 1.0;
    int y = 1.0;

    init_sod_2d_3(&sph, x, y, mass);

    write_csv(&sph, output_filename);

    printf("\n初始化完成 ...\n");
    printf("輸出檔案至: %s\n", output_filename);
    printf("設定粒子質量: %f\n", mass);

    free_sph_system(&sph);
    return 0;
}

