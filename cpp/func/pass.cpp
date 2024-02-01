#include "base.h"

void pass_val(Base b) {
     b.log("pass_val");
}

void pass_ref(Base &b) {
     b.log("pass_ref");
}

void pass_ptr(Base *pb) {
     pb->log("pass_ptr");
}

int main() {
     cout << "Construct:" << endl;
     Base b(2);
     cout << endl;

     cout << "Pass by value:" << endl;
     pass_val(b);
     cout << endl;

     cout << "Pass by reference:" << endl;
     pass_ref(b);
     cout << endl;

     cout << "Pass by pointer:" << endl;
     pass_ptr(&b);
     cout << endl;

     cout << "Destruct:" << endl;
}
