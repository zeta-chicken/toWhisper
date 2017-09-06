/*
 * wave.h
 * author:zeta
 *
 * waveファイルを扱うヘッダファイル
 *
 * copy系統
 * VECTOR*の指すデータを複製する
 * (アドレスの指すデータ自身の複製)
 * copyFoo(VECTOR* src, VECTOR* des)
 * {
 *     memcpy(des, src, size);
 * }
 * 今後も参照するデータに対して使う
 *
 * move系統
 * VECTOR*の指すデータを移動する
 * (アドレスのみのコピー)
 * moveFoo(VECTOR* src, VECTOR* des)
 * {
 *     des = src;
 *     src = NULL;
 * }
 * 今後参照することのないデータに対して使う
 *
 */

#ifndef _WAVE_H_20170809_
#define _WAVE_H_20170809_

#include "vector.h"
#include <stdint.h>

typedef struct _wave WAVE;
typedef enum _channel {
    MONO   = 1,
    STEREO = 2
} CHANNEL;
typedef enum _chtype {
    L    = 0,
    R    = 1,
} CH_TYPE;

//サンプリング周波数，量子化ビット数，チャンネル
WAVE* initWave(uint32_t, uint16_t, CHANNEL);
int freeWave(WAVE*);
int copySet(WAVE*, const VECTOR*, CH_TYPE);
int moveSet(WAVE*, VECTOR*, CH_TYPE);
VECTOR* copyFrom(const WAVE*, CH_TYPE);
VECTOR* moveFrom(WAVE*, CH_TYPE);
WAVE* readWave(const char*);
int writeWave(const WAVE*, const char*);
uint32_t getSamplerate(const WAVE*);
CHANNEL getChannels(const WAVE*);
uint16_t getBit(const WAVE*);
int setSamplerate(WAVE*, uint32_t);
int setChannels(WAVE*, CHANNEL);
int setBit(WAVE*, uint16_t);

#endif
