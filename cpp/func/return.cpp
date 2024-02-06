#include "base.h"

Base return_val(int i) {
    Base b(i);
    b.log("return_val");
    return b;
}
/*
Base & return_ref(int i) {
    Base b(i);
    b.log("return_ref");
    return b;
}

Base * return_ptr(int i) {
    Base b(i);
    b.log("return_ptr");
    return &b;
}
*/
class TestClass {
    private:
        Base b;
    public:
        TestClass(int i): b(i){};
        ~TestClass() {};
        const Base& return_ref() const {
            b.log("return_ref");
            return b;
        }
        const Base* return_ptr() const {
            b.log("return_ptr");
            return &b;
        }
};

int main() {
    cout << "Return by value:" << endl;
    Base b1 = return_val(2);
    b1.log("access_val");
    cout << endl;

    cout << "Return by reference:" << endl;
    //Base &b2 = return_ref(3);
    TestClass t2(3);
    const Base &b2 = t2.return_ref();
    b2.log("access_ref");
    cout << endl;

    cout << "Return by pointer:" << endl;
    //Base *pb3 = return_ptr(4);
    TestClass t3(4);
    const Base *pb3 = t3.return_ptr();
    pb3->log("access_ptr");
    cout << endl;

    cout << "Destruction:" << endl;
}
