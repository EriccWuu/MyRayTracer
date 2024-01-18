#include <iostream>
#include <cstring>
#include "TGAImage.h"

TGAImage::TGAImage(const int w, const int h, const int pd): w(w), h(h), pd(pd), data(w*h*pd, 0) {}


bool TGAImage::read_tga_file(const std::string filepath) {
    std::ifstream in;
    in.open(filepath, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "can't open file " << filepath << "\n";
        return false;
    }

    TGAHeader header;
    in.read(reinterpret_cast<char *>(&header), sizeof(header));
    if (!in.good()) {
        std::cerr << "an error occured while reading the header.\n";
        return false;
    }

    w = header.width;
    h = header.height;
    pd = header.pixDepth >> 3;
    if (w <= 0 || h <= 0 || (pd != GRAYSCALE && pd != RGB && pd != RGBA)) {
        std::cerr << "bad pixdepth (or width/ height) value.\n";
        return false;
    }
    size_t nbytes = w*h*pd;
    data = std::vector<BYTE>(nbytes, 0);
    if (header.imageType == 3 || header.imageType == 2) {
        in.read(reinterpret_cast<char *>(data.data()), nbytes);
        if (!in.good()) {
            std::cerr << "an error occured while reading data.\n";
            return false;
        }
    }
    else if (header.imageType == 10 || header.imageType == 11) {
        if (!load_rle_data(in)) {
            std::cerr << "an error occured while reading data.\n";
            return false;
        }
    }
    else {
        std::cerr << "unknown file format " << (int)header.imageType << "\n";
        return false;
    }

    if (!(header.imageDescriptor & 0x20))
        flip_vertically();
    else if (!(header.imageDescriptor & 0x10))
        flip_horizontally();

    std::cerr << w << "x" << h << "/" << (pd<<3) << "\n";
    return true;
}

bool TGAImage::write_tga_file(const std::string filepath, const bool vflip, const bool rle) const {
    constexpr BYTE developer_area_ref[4] = {0, 0, 0, 0};
    constexpr BYTE extension_area_ref[4] = {0, 0, 0, 0};
    constexpr std::uint8_t footer[18] = {'T','R','U','E','V','I','S','I','O','N','-','X','F','I','L','E','.','\0'};

    std::ofstream out;
    out.open(filepath, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "can't open file " << filepath << "\n";
        return false;
    }

    TGAHeader header = {};
    header.pixDepth = pd << 3;
    header.width = w;
    header.height = h;
    header.imageType = (pd == GRAYSCALE ? (rle ? 11 : 3) : (rle ? 10 : 2));
    header.imageDescriptor = vflip ? 0x00 : 0x20;   // 0x20: top-left origin, 0x00: bottom-left origin
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));
    if (!out.good()) {
        std::cerr << "can't dump the tga file.\n";
        return false;
    }

    if (!rle) {
        out.write(reinterpret_cast<const char *>(data.data()), w*h*pd);
        if (!out.good()) {
            std::cerr << "can't unload raw data.\n";
            return false;
        }
    }
    else if (!unload_rle_data(out)) {
        std::cerr << "can't unload rle data.\n";
        return false;
    }

    out.write(reinterpret_cast<const char *>(developer_area_ref), sizeof(developer_area_ref));
    if (!out.good()) {
        std::cerr << "can't dump the tga file.\n";
        return false;
    }
    out.write(reinterpret_cast<const char *>(extension_area_ref), sizeof(extension_area_ref));
    if (!out.good()) {
        std::cerr << "can't dump the tga file.\n";
        return false;
    }
    out.write(reinterpret_cast<const char *>(footer), sizeof(footer));
    if (!out.good()) {
        std::cerr << "can't dump the tga file.\n";
        return false;
    }

    return true;
}

// load rle compressed data
bool TGAImage::load_rle_data(std::ifstream &in) {
    size_t pixelcount = w*h;
    size_t currentpixel = 0;
    size_t currentbyte  = 0;
    TGAColor colorbuffer;
    do {
        std::uint8_t chunkheader = 0;
        chunkheader = in.get();
        if (!in.good()) {
            std::cerr << "an error occured while reading the data\n";
            return false;
        }
        if (chunkheader<128) {
            chunkheader++;
            for (int i=0; i<chunkheader; i++) {
                in.read(reinterpret_cast<char *>(colorbuffer.bgra), pd);
                if (!in.good()) {
                    std::cerr << "an error occured while reading the header\n";
                    return false;
                }
                for (int t=0; t<pd; t++)
                    data[currentbyte++] = colorbuffer.bgra[t];
                currentpixel++;
                if (currentpixel>pixelcount) {
                    std::cerr << "Too many pixels read\n";
                    return false;
                }
            }
        } else {
            chunkheader -= 127;
            in.read(reinterpret_cast<char *>(colorbuffer.bgra), pd);
            if (!in.good()) {
                std::cerr << "an error occured while reading the header\n";
                return false;
            }
            for (int i=0; i<chunkheader; i++) {
                for (int t=0; t<pd; t++)
                    data[currentbyte++] = colorbuffer.bgra[t];
                currentpixel++;
                if (currentpixel>pixelcount) {
                    std::cerr << "Too many pixels read\n";
                    return false;
                }
            }
        }
    } while (currentpixel < pixelcount);
    return true;
}

// unload rle compressed data
bool TGAImage::unload_rle_data(std::ofstream &out) const {
    const std::uint8_t max_chunk_length = 128;
    size_t npixels = w*h;
    size_t curpix = 0;
    while (curpix<npixels) {
        size_t chunkstart = curpix*pd;
        size_t curbyte = curpix*pd;
        std::uint8_t run_length = 1;
        bool raw = true;
        while (curpix+run_length<npixels && run_length<max_chunk_length) {
            bool succ_eq = true;
            for (int t=0; succ_eq && t<pd; t++)
                succ_eq = (data[curbyte+t]==data[curbyte+t+pd]);
            curbyte += pd;
            if (1==run_length)
                raw = !succ_eq;
            if (raw && succ_eq) {
                run_length--;
                break;
            }
            if (!raw && !succ_eq)
                break;
            run_length++;
        }
        curpix += run_length;
        out.put(raw?run_length-1:run_length+127);
        if (!out.good()) {
            std::cerr << "can't dump the tga file\n";
            return false;
        }
        out.write(reinterpret_cast<const char *>(data.data()+chunkstart), (raw?run_length*pd:pd));
        if (!out.good()) {
            std::cerr << "can't dump the tga file\n";
            return false;
        }
    }
    return true;
}

TGAColor TGAImage::get(const int x, const int y) const {
    if (!data.size() || x < 0 || y < 0 || x >= w | y >= h)
        return {};
    TGAColor ret = {0, 0, 0, 0};
    const BYTE *p = data.data() + (x + y * w) * pd;
    for (int i = pd; i --; ret.bgra[i] = p[i]);
    
    return ret;
}

void TGAImage::set(int x, int y, const TGAColor &c) {
    if (!data.size() || x < 0 || y < 0 || x >= w || y >= h) return;
    memcpy(data.data() + (x + y * w) * pd, c.bgra, pd);
}

void TGAImage::flip_horizontally() {
    int half = w >> 1;
    for (int i = 0; i < half; i ++)
        for (int j = 0; j < h; j ++)
            for (int k = 0; k < pd; k ++) 
                std::swap(data[(i+j*w)*pd+k], data[(w-1-i+j*w)*pd+k]);
}

void TGAImage::flip_vertically() {
    int half = h >> 1;
    for (int i = 0; i < w; i ++)
        for (int j = 0; j < half; j ++)
            for (int k = 0; k < pd; k ++) 
                std::swap(data[(i+j*w)*pd+k], data[(i+(h-1-j)*w)*pd+k]);
}

int TGAImage::width() const {
    return w;
}

int TGAImage::height() const {
    return h;
}

