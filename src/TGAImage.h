#ifndef __TGAIMAGE_H__
#define __TGAIMAGE_H__

#pragma once
#include <cstdint>
#include <fstream>
#include <vector>
#include "Matrix.h"

typedef std::uint8_t BYTE;
typedef std::uint16_t WORD;

#pragma pack(push, 1)
typedef struct TGAHeader {
    BYTE idLength = 0;          /* 00h  Size of Image ID field */
    BYTE colorMapType = 0;      /* 01h  Color map type */
    BYTE imageType = 0;         /* 02h  Image type code */
    WORD colorMapOrigin = 0;    /* 03h  Color map origin */
    WORD colorLength = 0;       /* 05h  Color map length */
    BYTE colorDepth = 0;        /* 07h  Depth of color map entries */
    WORD xorigin = 0;           /* 08h  X origin of image */
    WORD yorigin = 0;           /* 0Ah  Y origin of image */
    WORD width = 0;             /* 0Eh  Height of image */
    WORD height = 0;            /* 10h  Image pixel size */
    BYTE pixDepth = 0;          /* 11h  Image descriptor byte */
    BYTE imageDescriptor = 0;
} TGA_HEADER;
#pragma pack(pop)

struct TGAColor {
    BYTE bgra[4] = {0, 0, 0, 0};
    BYTE bytespp = 4;
    WORD val;

    TGAColor() : bgra(), bytespp(1) {
        for (int i = 0; i < 4; i ++) bgra[i] = 0;
    }

    TGAColor(const BYTE r, const BYTE g, const BYTE b, const BYTE a = 255) {
        bgra[0] = b;
        bgra[1] = g;
        bgra[2] = r;
        bgra[3] = a;
    }

    TGAColor(const vec4 &color) {
        bgra[0] = color.z;
        bgra[1] = color.y;
        bgra[2] = color.x;
        bgra[3] = color.w;
    }

    BYTE& operator[] (const int i) { return bgra[i]; }

    TGAColor& operator= (const vec4 &c) {
        for (int i = 2; i >= 0; i --) this->bgra[i] = (BYTE)c[2-i];
        this->bgra[3] = c.z ? c.z : 255;
        return *this;
    }

    bool operator< (const TGAColor &c) const {
        for (int i = 0; i < 3; i ++) if (bgra[i] > c.bgra[i]) return false;
        return true;
    }

    bool operator> (const TGAColor &c) const {
        for (int i = 0; i < 3; i ++) if (bgra[i] < c.bgra[i]) return false;
        return true;
    }

    TGAColor operator+ (TGAColor &c) const {
        TGAColor res = *this;
        for (int i = 0; i < 3; i ++) res.bgra[i] = std::min(res.bgra[i] + c.bgra[i], 255);
        return res;
    }

    TGAColor operator* (float intensity) const {
        TGAColor res = *this;
        intensity = (intensity > 1.f ? 1.f : (intensity < 0.f ? 0.f : intensity));
        for (int i = 0; i < 3; i ++) res.bgra[i] = bgra[i] * intensity;
        return res;
    }
};

struct TGAImage {
    enum Format { GRAYSCALE = 1, RGB = 3, RGBA = 4 };
    
    TGAImage() = default;
    TGAImage(const int w, const int h, const int pd);
    bool read_tga_file(const std::string filepath);
    bool write_tga_file(const std::string filepath, const bool vfilp = false, const bool rle = true) const;
    void flip_horizontally();
    void flip_vertically();
    TGAColor get(const int x, const int y) const;
    void set(int x, int y, const TGAColor &c);
    int width() const;
    int height() const;

private:
    bool load_rle_data(std::ifstream &in);
    bool unload_rle_data(std::ofstream &out) const;

    int w = 0;
    int h = 0;
    BYTE pd = 0;
    std::vector<BYTE> data = {};
};

#endif
