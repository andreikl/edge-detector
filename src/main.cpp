#include <iostream>
#include <vector>
#include <map>

using namespace std;

struct RGB {
    //char r;
    //char g;
    //char b;
    char gr;
};

int printHelp() {
    cout << "edge-detector [-w] [-h]" << endl;
    cout << "-w width, default 640" << endl;
    cout << "-h height, default 480" << endl;
    cout << "-f frames, default 1" << endl;
    exit(0);
}

const string WIDTH("-w");
const string HEIGHT("-h");
const string FRAMES("-f");

int main(int argc, char *argv[]) {
    map <string, string> argmap;
    for (int i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            argmap[argv[i]] = argv[i + 1];
        }
    }

    auto getWidth = [=] { return argmap.count(WIDTH) ? stoi(argmap.at(WIDTH)): 640; };
    auto getHeight = [=] { return argmap.count(HEIGHT) ? stoi(argmap.at(HEIGHT)) : 480; };
    auto getFrames = [=] { return argmap.count(FRAMES) ? stoi(argmap.at(FRAMES)) : 1; };
    const int width = getWidth();
    const int height = getHeight();
    const int frames = getFrames();

    vector<RGB> image;
    image.reserve(width * height);
    //std::cout << "Frame size (3x): " << image.size() << std:endl;

    for (int k = 0; k < frames; k++) {
        image.clear();
        //reads one frame to array
        std::cin.read((char*)image.data(), width * height * sizeof(RGB));

        if (k == frames - 1) {
            //writes ppm image
            std::cout << "P6\n" << width << " " << height << "\n255\n";
            //std::cout.write((char*)image.data(), width * height * sizeof(RGB));
            //std::cout.flush();
            for (int i = width * height - 1; i >= 0; i--) {
                //std::cout << image[j].r;
                //std::cout << image[j].g;
                //std::cout << image[j].b;
                std::cout << image[i].gr;
                std::cout << image[i].gr;
                std::cout << image[i].gr;
            }
        }
    }
}
