#include <iostream>

int main(int argc, char *argv[]) {
    std::vector<int> numbers = { 1, 2, 3, 4 };
    for (int i : numbers) {
        std::cout << i << '\n';
    }
}