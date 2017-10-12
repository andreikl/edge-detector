#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
    auto x = [=] { return argc; };
    
    std::vector<int> numbers = { 1, 2, 3, 4 };
    for (int i : numbers) {
        std::cout << i << '\n';
    }
    std::cout << x();
}