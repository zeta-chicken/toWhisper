#include "wave.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))
#define CLAMP(a,b,c) (MIN(MAX((a),(b)),(c)))

struct _wave {
    VECTOR*   data[2];
    uint32_t  fs;
    uint16_t  bit;
    CHANNEL   ch;
};

struct RiffHeader {
    uint8_t riff[4];
    int32_t size;
    uint8_t type[4];
};

struct FormatChunk {
    uint8_t  id[4];
    int32_t  size;
    int16_t  format;
    uint16_t channels;
    uint32_t samplerate;
    uint32_t bytepersec;
    uint16_t blockalign;
    uint16_t bitswidth;
};

struct DataChunk {
    uint8_t id[4];
    int32_t size;
    uint8_t *waveformData;
};

enum HEADER_TYPE {
    RIFF_HEADER,
    FORMAT_CHUNK,
    DATA_CHUNK,
    OTHER,
};

//波形データの相対参照
static size_t refData(const struct _wave* wave, double* elem[2])
{
    if (wave->data[0] == NULL && wave->data[1] == NULL)
        return 0;

    if (wave->ch == STEREO) {
        if (wave->data[0] == NULL || wave->data[1] == NULL)
            return 0;

        elem[0] = wave->data[0]->elem;
        elem[1] = wave->data[1]->elem;
        return MIN(wave->data[0]->sz, wave->data[1]->sz);
    }

    if (wave->ch == MONO) {
        if (wave->data[0] == NULL) {
            elem[0] = wave->data[1]->elem;
            elem[1] = NULL;
            return wave->data[1]->sz;
        }

        if (wave->data[1] == NULL) {
            elem[0] = wave->data[0]->elem;
            elem[1] = NULL;
            return wave->data[0]->sz;
        }

        elem[0] = wave->data[0]->elem;
        elem[1] = wave->data[1]->elem;
        return MIN(wave->data[0]->sz, wave->data[1]->sz);
    }

    return 0;
}

enum HEADER_TYPE checkHeaderType(const uint8_t id[4])
{
    if (!strncmp((char*)id, "RIFF", 4)) return RIFF_HEADER;
    if (!strncmp((char*)id, "fmt ", 4)) return FORMAT_CHUNK;
    if (!strncmp((char*)id, "data", 4)) return DATA_CHUNK;
    return OTHER;
}

static int writeRIFF(const struct RiffHeader* RIFF, FILE* fp)
{
    if (fp == NULL) return -1;
    fwrite(RIFF->riff,  1, 4, fp);
    fwrite(&RIFF->size, 4, 1, fp);
    fwrite(RIFF->type,  1, 4, fp);
    return 0;
}

static int readRIFF(struct RiffHeader* RIFF, FILE* fp)
{
    if (fp == NULL) return -1;
    fread(RIFF->riff,  1, 4, fp);
    fread(&RIFF->size, 4, 1, fp);
    fread(RIFF->type,  1, 4, fp);
    return 0;
}

static int writeFmt(const struct FormatChunk* fmt, FILE* fp)
{
    if (fp == NULL) return -1;
    fwrite(fmt->id,          1, 4, fp);
    fwrite(&fmt->size,       4, 1, fp);
    fwrite(&fmt->format,     2, 1, fp);
    fwrite(&fmt->channels,   2, 1, fp);
    fwrite(&fmt->samplerate, 4, 1, fp);
    fwrite(&fmt->bytepersec, 4, 1, fp);
    fwrite(&fmt->blockalign, 2, 1, fp);
    fwrite(&fmt->bitswidth,  2, 1, fp);
    return 0;
}

static int readFmt(struct FormatChunk* fmt, FILE* fp)
{
    if (fp == NULL) return -1;
    fread(fmt->id,          1, 4, fp);
    fread(&fmt->size,       4, 1, fp);
    fread(&fmt->format,     2, 1, fp);
    fread(&fmt->channels,   2, 1, fp);
    fread(&fmt->samplerate, 4, 1, fp);
    fread(&fmt->bytepersec, 4, 1, fp);
    fread(&fmt->blockalign, 2, 1, fp);
    fread(&fmt->bitswidth,  2, 1, fp);
    return 0;
}

static int writeData(const struct DataChunk* data, FILE* fp)
{
    if (fp == NULL) return -1;
    fwrite(data->id,           1, 4,          fp);
    fwrite(&data->size,        4, 1,          fp);
    fwrite(data->waveformData, 1, data->size, fp);
    return 0;
}

static int readData(struct DataChunk* data, FILE* fp)
{
    if (fp == NULL) return -1;
    fread(data->id,           1, 4,          fp);
    fread(&data->size,        4, 1,          fp);
    data->waveformData = malloc(data->size);
    fread(data->waveformData, 1, data->size, fp);
    return 0;
}

static int double2charWave(uint8_t* des, double* src[2], uint32_t sz, CHANNEL ch, uint16_t bit)
{
    for (int i=0; i<sz; i++) {
        if (bit == 8) {
            switch (ch) {
            case MONO:
                if (src[1] != NULL)
                    des[i] =
                        round(255.0 * CLAMP((src[0][i] + src[1][i])/4.0 + 0.5, 0.0, 1.0));
                else
                    des[i] = round(255.0 * CLAMP(src[0][i]/2.0 + 0.5, 0.0, 1.0));
                break;

            case STEREO:
                des[2*i]   = round(255.0 * CLAMP(src[0][i]/2.0 + 0.5, 0.0, 1.0));
                des[2*i+1] = round(255.0 * CLAMP(src[1][i]/2.0 + 0.5, 0.0, 1.0));
                break;
            }
        }

        if (bit == 16) {
            int16_t heap = 0;
            switch (ch) {
            case MONO:
                if (src[1] != NULL)
                    heap =
                        CLAMP(round(32767.0 * (src[0][i]+src[1][i])/2.0), -32767.0, 32767.0);
                else
                    heap = CLAMP(round(32767.0*src[0][i]), -32767.0, 32767.0);
                memcpy(des+2*i, &heap, 2);
                break;

            case STEREO:
                heap = CLAMP(round(32767.0*src[0][i]), -32767.0, 32767.0);
                memcpy(des+4*i, &heap, 2);
                heap = CLAMP(round(32767.0*src[1][i]), -32767.0, 32767.0);
                memcpy(des+4*i+2, &heap, 2);
                break;
            }
        }
    }

    return 0;
}

static int char2doubleWave(double* src[2], const uint8_t* des, uint32_t sz, CHANNEL ch, uint16_t bit)
{
    for (int i=0; i<sz; i++) {
        if (bit == 8) {
            uint8_t heap = 0;
            switch (ch) {
            case MONO:
                heap = des[i];
                src[0][i] = (double)heap / 127.5 - 1.0;
                break;

            case STEREO:
                heap = des[2*i];
                src[0][i] = (double)heap / 127.5 - 1.0;
                heap = des[2*i+1];
                src[1][i] = (double)heap / 127.5 - 1.0;
                break;
            }
        }

        if (bit == 16) {
            int16_t heap = 0;
            switch (ch) {
            case MONO:
                memcpy(&heap, des+2*i, 2);
                src[0][i] = (double)heap / 32767.0;
                break;

            case STEREO:
                memcpy(&heap, des+4*i, 2);
                src[0][i] = (double)heap / 32767.0;
                memcpy(&heap, des+4*i+2, 2);
                src[1][i] = (double)heap / 32767.0;
                break;
            }
        }
    }
    return 0;
}

struct _wave* initWave(uint32_t fs, uint16_t bit, CHANNEL ch)
{
    struct _wave* wave = malloc(sizeof(struct _wave));
    if (wave == NULL) {
        fprintf(stderr, "In function, initWave: \
                can not allocate wave memory. return null pointer.");
        return wave;
    }

    wave->data[0] = NULL;
    wave->data[1] = NULL;
    wave->fs   = fs;
    wave->bit  = bit;
    wave->ch   = ch;
    return wave;
}

int freeWave(struct _wave* wave)
{
    if (wave == NULL) return 0;
    freeVector(wave->data[0]);
    freeVector(wave->data[1]);
    free(wave);
    wave = NULL;
    return 0;
}

int copySet(struct _wave* wave, const VECTOR* data, CH_TYPE type)
{
    if (data == NULL) {
        fprintf(stderr, "In function, copySet: 2nd argument, VECTOR* is null pointer.\n");
        return -1;
    }

    if (wave == NULL) {
        fprintf(stderr, "In function, copySet: 1st argument WAVE* is null pointer.\n");
        return -1;
    }

    if (type>1) {
        fprintf(stderr, "In function, copySet: 3rd argument CH_TYPE is not allowed.\n");
        return -1;
    }

    freeVector(wave->data[type]);
    wave->data[type] = initVector(data->sz);
    memcpy(wave->data[type]->elem, data->elem, data->sz*sizeof(double));
    return 0;
}

int moveSet(struct _wave* wave, VECTOR* data, CH_TYPE type)
{
    if (data == NULL) {
        fprintf(stderr, "In function, moveSet: 2nd argument, VECTOR* is null pointer.\n");
        return -1;
    }

    if (wave == NULL) {
        fprintf(stderr, "In function, moveSet: 1st argument WAVE* is null pointer.\n");
        return -1;
    }

    if (type>1) {
        fprintf(stderr, "In function, moveSet: 3rd argument CH_TYPE is not allowed.\n");
        return -1;
    }

    freeVector(wave->data[type]);
    wave->data[type] = data;
    data = NULL;
    return 0;
}

VECTOR* copyFrom(const struct _wave* wave, CH_TYPE type)
{
    if (wave == NULL) {
        fprintf(stderr, "In function, copyFrom: 1st argument WAVE* is null pointer.\n");
        return NULL;
    }

    if (type>1) {
        fprintf(stderr, "In function, copyFrom: 2nd argument CH_TYPE is not allowed.\n");
        return NULL;
    }

    VECTOR* v = initVector(wave->data[type]->sz);

    if (v == NULL) {
        fprintf(stderr, "In function, copyFrom: \
                can not allocate vector memory. return null pointer.\n");
        return v;
    }

    memcpy(v->elem, wave->data[type]->elem, v->sz*sizeof(double));
    return v;
}

VECTOR* moveFrom(struct _wave* wave, CH_TYPE type)
{
    if (wave == NULL) {
        fprintf(stderr, "In function, moveFrom: 1st argument WAVE* is null pointer.\n");
        return NULL;
    }

    if (type>1) {
        fprintf(stderr, "In function, moveFrom: 2nd argument CH_TYPE is not allowed.\n");
        return NULL;
    }

    VECTOR* v = wave->data[type];
    wave->data[type] = NULL;
    return v;
}

struct _wave* readWave(const char* file)
{
    FILE* fp = fopen(file, "rb");
    if (fp == NULL) {
        fprintf(stderr, "In function, readWave: can not open %s.\n", file);
        return NULL;
    }

    struct RiffHeader  RIFF = {{0}};
    struct FormatChunk fmt  = {{0}};
    struct DataChunk   data = {{0}};

    uint8_t id[4] = {0};
    while (fread(id, 1, 4, fp) == 4) {
        fseek(fp, -4, SEEK_CUR);

        switch (checkHeaderType(id)) {
        case RIFF_HEADER:
            readRIFF(&RIFF, fp);
            if (strncmp((char*)RIFF.type, "WAVE", 4)) {
                fprintf(stderr, "In function, readWave: %s is not a WAVE format file.\n", file);
                return NULL;
            }
            break;

        case FORMAT_CHUNK:
            readFmt(&fmt, fp);
            break;

        case DATA_CHUNK:
            readData(&data, fp);
            break;

        default:
            fprintf(stderr, "In function, readWave: Sorry, this file can not be used.\n");
            return NULL;
            break;
        }

        if (RIFF.riff[0] * fmt.id[0] * data.id[0] != 0) break;

        id[0] = 0;
        id[1] = 0;
        id[2] = 0;
        id[3] = 0;
    }

    /*
    if (feof(fp) == 0) {
        fprintf(stderr, "In function, readWave: can not read this file.\n");
        return NULL;
    }
    */

    struct _wave* wave = initWave(fmt.samplerate, fmt.bitswidth, fmt.channels);
    if (wave == NULL) {
        fprintf(stderr, "In function, readWave: can not allocate wave memory.\n");
        return NULL;
    }

    VECTOR* v[2] = {initVector(data.size/fmt.blockalign), initVector(data.size/fmt.blockalign)};
    double* elem[2] = {v[0]->elem, v[1]->elem};
    char2doubleWave(elem, data.waveformData, v[0]->sz, wave->ch, wave->bit);
    free(data.waveformData);
    moveSet(wave, v[0], L);
    moveSet(wave, v[1], R);

    fclose(fp);

    return wave;
}

int writeWave(const struct _wave* wave, const char* file)
{
    double* elem[2];
    const size_t sz = refData(wave, elem);
    if (sz == 0) return -1;

    FILE* fp = fopen(file, "wb");
    if (fp == NULL) return -1;

    //RIFFヘッダ
    struct RiffHeader RIFF = {
        .riff = {'R', 'I', 'F', 'F'},
        .size = 36 + wave->bit/8 * wave->ch * sz,
        .type = {'W', 'A', 'V', 'E'}
    };

    //fmtチャンク
    struct FormatChunk fmt = {
        .id         = {'f', 'm', 't', ' '},
        .size       = 16,
        .format     = 1,
        .channels   = wave->ch,
        .samplerate = wave->fs,
        .bytepersec = wave->fs * wave->bit/8 * wave->ch,
        .blockalign = wave->bit/8 * wave->ch,
        .bitswidth  = wave->bit
    };

    //dataチャンク
    struct DataChunk data = {
        .id           = {'d', 'a', 't', 'a'},
        .size         = wave->bit/8 * wave->ch * sz,
        .waveformData = malloc(wave->bit/8 * wave->ch * sz)
    };

    double2charWave(data.waveformData, elem, sz, wave->ch, wave->bit);

    writeRIFF(&RIFF, fp);
    writeFmt(&fmt, fp);
    writeData(&data, fp);

    free(data.waveformData);
    fclose(fp);

    return 0;
}

uint32_t getSamplerate(const struct _wave* wave)
{
    return wave->fs;
}

CHANNEL getChannels(const struct _wave* wave)
{
    return wave->ch;
}

uint16_t getBit(const struct _wave* wave)
{
    return wave->bit;
}

int setSamplerate(struct _wave* wave, uint32_t fs)
{
    wave->fs = fs;
    return 0;
}

int setChannels(struct _wave* wave, CHANNEL ch)
{
    if (ch > R) return -1;
    wave->ch = ch;
    return 0;
}

int setBit(struct _wave* wave, uint16_t bit)
{
    if (bit != 8 || bit != 16) return -1;
    wave->bit = bit;
    return 0;
}
