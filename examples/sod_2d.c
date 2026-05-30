# include "sph_system.h"

int main(int argc, char *argv[]) {
    printf("====================================\n");
    printf("   SPH 模擬程式啟動中...             \n");
    printf("====================================\n");
    // 1. 設定預設值
    const char *output_filename = "output_0000.csv";
    double mass = 0.001;  // 預設質量

    // 2. 解析命令列參數
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_filename = argv[i + 1];
            i++; // 跳過已讀取的檔名
        }
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            mass = atof(argv[i + 1]); // 將字串轉換為 double
            i++; // 跳過已讀取的數值
        }
    }

    SPHSystem2D sph;
    int x = 1.0;
    int y = 1.0;

    // 如果你的 init_sod_2d_2 需要用到 mass，記得要把 mass 傳進去，
    // 例如：init_sod_2d_2(&sph, x, y, mass);
    // 這裡先維持你原本的呼叫方式：
    init_sod_2d_2(&sph, x, y, mass);

    write_csv(&sph, output_filename);

    printf("\n初始化完成 ...\n");
    printf("輸出檔案至: %s\n", output_filename);
    printf("設定粒子質量: %f\n", mass);

    free_sph_system(&sph);
    return 0;
}

