#include "base.h"

Base return_val(int i) {
    Base b(i);
    b.log("return_val");
    return b;
}

// Base & return_ref(int i) {
//     Base b(i);
//     b.log("return_val");
//     return b;
// }

// Base * return_ptr(int i) {
//     Base b(i);
//     b.log("return_val");
//     return &b;
// }

int main() {
    cout << "Return by value:" << endl;
    Base b1 = return_val(2);
    b1.log("access_val");
    cout << endl;

    cout << "Return by reference:" << endl;
    //Base &b2 = return_ref(3);
    //b2.log("access_ref");
    cout << endl;

    cout << "Return by pointer:" << endl;
    //Base *pb3 = return_ptr(4);
    //pb3->log("access_ptr");
    cout << endl;

    cout << "Destruction:" << endl;
}
