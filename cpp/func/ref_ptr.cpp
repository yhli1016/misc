#include "base.h"


void pass_ref(Base &rb, Base &rc) {
     std::cout << "Before swap:\n";
     rb.log("value_b");
     rc.log("value_c");
     Base tmp = rb;
     rb = rc;
     rc = tmp;
     std::cout << "After swap:\n";
     rb.log("value_b");
     rc.log("value_c");
}

void pass_ptr(Base *pb, Base *pc) {
     std::cout << "Before swap:\n";
     pb->log("value_b");
     pc->log("value_c");
     //(*pb).log("value_b");
     //(*pc).log("valuc_c");
     Base tmp = *pb;
     *pb = *pc;
     *pc = tmp;
     std::cout << "After swap:\n";
     pb->log("value_b");
     pc->log("value_c");
     //(*pb).log("value_b");
     //(*pc).log("value_c");
}

int main() {
     cout << "Construct:" << endl;
     Base b(2), c(3);
     Base &ref_b = b, &ref_c = c;
     Base *ptr_b = &b, *ptr_c = &c;
     cout << endl;

     cout << "Implicit reference:" << endl;
     pass_ref(b, c);
     cout << endl;
     cout << "Explicit reference:" << endl;
     pass_ref(ref_b, ref_c);
     cout << endl;

     cout << "Implicit pointer:" << endl;
     pass_ptr(&b, &c);
     cout << endl;
     cout << "Explicit pointer:" << endl;
     pass_ptr(ptr_b, ptr_c);
     cout << endl;

     cout << "Destruct:" << endl;
}
