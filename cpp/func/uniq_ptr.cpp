#include <memory>
#include "base.h"


void vector_val(int n, std::unique_ptr<Base[]> &vec) {
    for(int i = 0; i < n; i++){
        cout << "i = " << i << endl;
        vec[i] = Base(i);
    }
}

int main() {
    // first approach
    std::unique_ptr<Base[]> vec(new Base[10]);

    // second approach
    //std::unique_ptr<Base[]> vec = nullptr;
    //vec = std::unique_ptr<Base[]>(new Base[10]);

    vector_val(10, vec);
    cout << "Destruction:" << endl;
}
