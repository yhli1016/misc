#include <vector>
#include "base.h"


void vector_val(int n, std::vector<Base> &vec) {
    for(int i = 0; i < n; i++){
        cout << "i = " << i << endl;
        vec[i] = Base(i);
    }
}

int main() {
    // first approach
    // std::vector<Base> vec(10);

    // second approach
    std::vector<Base> vec;
    vec = std::vector<Base>(10);

    vector_val(10, vec);
    cout << "Destruction:" << endl;
}
