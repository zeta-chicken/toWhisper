#include "wave.h"
#include "vector.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

double autocorrelation(double*, size_t, size_t);
void   LevinsonDurbin(double*, size_t, double*, size_t);
double genWhite(void);
double genRosenberg(int, double, uint32_t);
double hanning(int, int);
double hamming(int, int);
double blackman(int, int);
void   ltiFilter(double*, size_t, double*, size_t, double*, size_t);
void   printHelp(void);
int    initMain(int, char**);
int    setWindow(char*);

//LPC次数
int order = 0;
//有声音割合(0.0~1.0)
double rate = 0.0;
//プリエンファシスフィルタの係数(0.0~1.0)
double hpf = 0.97;
//デエンファシスフィルタの係数(0.0~1.0)
double lpf = 0.97;
//入出力ファイル
char *inputFile = NULL, *outputFile = NULL, *eFile = NULL;
//frame幅
double frameT = 20.0;
//各種フラグ
int printFlag = 0;
//窓関数タイプ
double (* const windowFunction[])(int, int) = {hamming, hanning, blackman};
enum {
    HAMMING,
    HANNING,
    BLACKMAN,
    OTHER,
} windowType = HAMMING;

int main(int argc, char* argv[])
{
    if (initMain(argc, argv) != 0) return -1;

    WAVE *src = readWave(inputFile);

    if (getChannels(src) != MONO) {
        return -1;
    }

    VECTOR *v1 = moveFrom(src, L);

    //LPC次数の計算
    if (order == 0)
        order = (double)getSamplerate(src) / 44100.0 * 40.0;

    //フレーム幅をサンプル数に
    int frame = (double)getSamplerate(src) * frameT / 1000.0;
    if (frame%2 != 0) frame++;

    //フレーム幅でちょうど割り切れるようにする
    int last = v1->sz;
    int length = v1->sz - (v1->sz%frame) + frame;

    VECTOR *v = initVector(length);

    //足りない分はゼロづめ
    for (int i=1; i<length; i++) {
        if (i<last) {
            v->elem[i] = v1->elem[i] - hpf*v1->elem[i-1];
        } else {
            v->elem[i] = 0.0;
        }
    }

    freeVector(v1);
    v1 = initVector(length);
    VECTOR *v2 = initVector(length);

    if (printFlag) {
        printf("Sampling Rate         : %dHz\n", getSamplerate(src));
        printf("Quantization Bit rate : %dbit\n\n", getBit(src));
        printf("Window Function       : ");
        switch (windowType) {
        case HAMMING:
            printf("hamming\n");
            break;
        case HANNING:
            printf("hanning\n");
            break;
        case BLACKMAN:
            printf("blackman\n");
            break;
        default:
            printf("\n");
            break;
        }
        printf("Frame Length          : %d\n", frame);
        printf("LPC Order             : %d\n\n", order);
        printf("%s --> %s\n", inputFile, outputFile);
        if (eFile != NULL) printf("%s --> %s\n", inputFile, eFile);
        puts("");
    }

    double *x = malloc(sizeof(double)*frame);
    double *y = malloc(sizeof(double)*frame);
    double *E = malloc(sizeof(double)*frame);       //残差信号
    double *a = malloc(sizeof(double)*(order+1));   //LPC係数

    for (int i=0; i<length/frame*2-1; i++) {
        double max = 0.0;
        for (int j=0; j<frame; j++) {
            x[j] = v->elem[j+i*frame/2] * windowFunction[windowType](j, frame);
        }

        //LPC係数の導出
        LevinsonDurbin(x, frame, a, order);

        //残差信号(声帯音源)の導出
        for (int j=0; j<frame; j++) {
            double e = 0.0;
            for (int n=0; n<order+1; n++) {
                if (j >= n) e += a[n]*x[j-n];
            }
            E[j] = e;
            max += e*e;
        }

        //声帯振動
        for (int j=0; j<frame; j++) v2->elem[j+i*frame/2] += E[j]*10.0;

        //ホワイトノイズの生成
        for (int j=0; j<frame; j++) {
            y[j] = genWhite();
        }

        //残差信号とノイズの二乗平均レベルをそろえる
        max = sqrt(3.0*max/(double)frame);
        for (int j=0; j<frame; j++) {
            y[j] = rate*E[j] + (1.0-rate)*max*y[j];
        }

        //音声合成フィルタ
        for (int j=0; j<frame; j++) {
            for (int n=1; n<order+1; n++) {
                if (j >= n) y[j] -= a[n]*y[j-n];
            }
        }

        for (int j=0; j<frame; j++) v1->elem[j+i*frame/2] += y[j];
    }

    //デエンファシス
    for (int i=1; i<length; i++) v1->elem[i] = lpf*v1->elem[i-1] + v1->elem[i];
    for (int i=1; i<length; i++) v2->elem[i] = 0.97*v2->elem[i-1] + v2->elem[i];

    moveSet(src, v1, L);
    writeWave(src, outputFile);

    moveSet(src, v2, L);
    if (eFile != NULL) writeWave(src, eFile);

    free(x);
    free(y);
    free(E);
    free(a);
    freeVector(v);
    freeWave(src);
}

// 自己相関関数
double autocorrelation(double* x, size_t l, size_t N)
{
    double res = 0.0;
    double r = 0.0, t = 0.0;
    for (int i=0; i<N-l; i++) {
        t = res + (x[i]*x[i+l]+r);
        r = (x[i]*x[i+l] + r) - (t - res);
        res = t;
    }
    return res;
}

/*
 * lpcOrder = k の場合
 * フィルタ係数は
 * a[0], a[1] , ......, a[k]
 * となるため，k+1のメモリ領域を確保して置く必要がある．
 */
void LevinsonDurbin(double* x, size_t length, double* a, size_t lpcOrder)
{
    double lambda = 0.0, E = 0.0;
    double* r = (double*)malloc(sizeof(double)*(lpcOrder+1));
    double* V = (double*)malloc(sizeof(double)*(lpcOrder+1));
    double* U = (double*)malloc(sizeof(double)*(lpcOrder+1));
    for (int i=0; i<lpcOrder+1; i++)
        r[i] = autocorrelation(x, i, length);
    for (int i=0; i<=lpcOrder; i++) a[i] = 0.0;
    a[0] = 1.0;
    a[1] = -r[1]/r[0];
    E = r[0] + r[1]*a[1];
    for (int k=1; k<lpcOrder; k++) {
        lambda = 0.0;
        for (int j=0; j<=k; j++) lambda += a[j]*r[k+1-j];
        lambda /= -E;
        for (int j=0; j<=k+1; j++) {
            U[j] = a[j];
            V[j] = a[k+1-j];
        }
        for (int j=0; j<=k+1; j++) a[j] = U[j]+lambda*V[j];
        E = (1.0-lambda*lambda)*E;
    }
    free(r);
    free(V);
    free(U);

    return;
}

double genWhite(void)
{
    return (double)rand()/(double)RAND_MAX*2.0-1.0;
}

double genRosenberg(int i, double f, uint32_t fs)
{
    double T = (double)fs / f;
    double tau1 = T*0.4, tau2 = T*0.16;

    i = i % (int)T;
    if ((double)i < tau1) return 3.0*pow((double)i/tau1, 2.0) - 2.0*pow((double)i/tau1, 3.0);
    if ((double)i < tau1 + tau2) return 1.0 - pow(((double)i-tau1)/tau2, 2.0);

    return 0.0;
}

double hanning(int i, int frame)
{
    return 0.5 - 0.5*cos(2.0*M_PI/(double)frame*(double)i);
}

double hamming(int i, int frame)
{
    return 0.54 - 0.46*cos(2.0*M_PI/(double)frame*(double)i);
}

double blackman(int i, int frame)
{
    return 0.42 - 0.5*cos(2.0*M_PI/(double)frame*(double)i) + 0.08*cos(4.0*M_PI/(double)frame*(double)i);
}

/*
 * 線形時不変フィルタ
 * a[0]は出力信号のフィルタ係数であるため，常に1にすること
 */
void ltiFilter(double* x, size_t len, double* a, size_t al, double* b, size_t bl)
{
    if (x == NULL) return;
    if (a == NULL) return;
    if (b == NULL) bl = 0;
    double *y = malloc(sizeof(double)*len);

    y[0] = x[0];
    for (size_t i=1; i<len; i++) {
        y[i] = 0.0;
        for (size_t j=0; j<bl; j++) {
            if (i >= j) y[i] += b[j]*x[i-j];
        }
        for (size_t j=1; j<al; j++) {
            if (i >= j) y[i] -= a[j]*y[i-j];
        }
    }

    memcpy(x, y, len*sizeof(double));
    free(y);

    return;
}

void printHelp(void)
{
    printf("The translating voice into whisper voice system\n");
    printf("Version 0.91\n");
    printf("Copyright (C) 2017 zeta\n");
    printf("All rights reserved.\n\n");
    printf("usage:\n");
    printf("\ttoWhisper [options] [infile]\n");
    printf("options:\n");
    printf("\t-o : output file name.               [N/A]     [wave file]\n");
    printf("\t-e : vowel sound file name.          [N/A]     [wave file]\n");
    printf("\t-l : de-emphasis filter coefficient. [0.97]    (-1.0, 1.0)\n");
    printf("\t-r : vowel sound rate.               [0.0]     [ 0.0, 1.0]\n");
    printf("\t-w : window function                 [hamming] [function]\n");
    printf("\t-f : frame length(ms).               [20]\n");
    printf("\t-O : LPC order.                      [auto]\n");
    printf("\t-p : print some information about synthesizing.\n");
    printf("\t-h : print this sentence.\n");
    printf("infile:\n");
    printf("\twave file (MONO File only)\n");
    printf("window function name:\n");
    printf("\thamming\n\thanning\n\tblackman\n\n");
    return;
}

int initMain(int argc, char* argv[])
{
    srand((unsigned)time(NULL));

    for (int i=1; i<argc; i++) {
        if (argv[i][0] != '-') {
            if (inputFile == NULL) {
                inputFile = argv[i];
            } else {
                printHelp();
                return -1;
            }
            continue;
        }

        switch (argv[i][1]) {
        case 'o':
            i++;
            outputFile = argv[i];
            break;
        case 'e':
            i++;
            eFile = argv[i];
            break;
        case 'l':
            i++;
            lpf = atof(argv[i]);
            break;
        case 'r':
            i++;
            rate = atof(argv[i]);
            break;
        case 'w':
            i++;
            windowType = setWindow(argv[i]);
            if (windowType == OTHER) {
                printHelp();
                return -1;
            }
            break;
        case 'f':
            i++;
            frameT = atof(argv[i]);
            break;
        case 'O':
            i++;
            order = atoi(argv[i]);
            break;
        case 'p':
            printFlag = 1;
            break;
        case 'h':
            printHelp();
            return -1;
            break;
        default:
            printHelp();
            return -1;
            break;
        }
    }

    if (inputFile == NULL) {
        printHelp();
        return -1;
    }
    if (outputFile == NULL) {
        printHelp();
        return -1;
    }

    return 0;
}

int setWindow(char* str)
{
    if (!strcmp(str, "hamming")) return HAMMING;
    if (!strcmp(str, "hanning")) return HANNING;
    if (!strcmp(str, "blackman")) return BLACKMAN;

    return OTHER;
}
