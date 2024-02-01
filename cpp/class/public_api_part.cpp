#include <iostream>
using std::cout;
using std::endl;

class Base {
    private:
        int a = 0;
    private:
        void log(const char* msg) const {
            cout << "\033[32m" << this << "\033[0m ";
            cout << "\033[34m" << a << "\033[0m ";
            cout << "\033[34m" << "Base " << msg << "\033[0m ";
            cout << endl;
        }
    public:
        Base(): a(1) {
            log("Default constructor");
        }
        Base(int _a): a(_a) {
            log("Parameter constructor");
        }
        Base(const Base &other) {
            a = other.a;
            log("Copy constructor");
        }
        ~Base() {
            log("Destructor");
        }
        virtual void print(const char* msg) const {
            cout << "\033[32m" << this << "\033[0m ";
            cout << "\033[34m" << a << "\033[0m ";
            cout << "\033[34m" << "Base " << msg << "\033[0m ";
            cout << endl;
        }
};

class Derived: public Base {
    private:
        int b = 0;
    private:
        void log(const char* msg) const {
            cout << "\033[32m" << this << "\033[0m ";
            cout << "\033[31m" << b << "\033[0m ";
            cout << "\033[31m" << "Derived " << msg << "\033[0m ";
            cout << endl;
        }
    public:
        Derived(): Base(), b(1) {
            log("Default constructor");
        }
        Derived(int _a, int _b): Base(_a), b(_b) {
            log("Parameter constructor");
        }
        Derived(const Derived &other): Base(other) {
            b = other.b;
            log("Copy constructor");
        }
        ~Derived() {
            log("Destructor");
        }
        virtual void print(const char* msg) const override {
            cout << "<<--------------" << endl;
            Base::print(msg);
            cout << "\033[32m" << this << "\033[0m ";
            cout << "\033[31m" << b << "\033[0m ";
            cout << "\033[31m" << "Derived " << msg << "\033[0m ";
            cout << endl;
            cout << "-------------->>" << endl;
        }
};

#ifdef TMPL
    // --------- templates for all classes --------
    template <typename T>
    void test_val(T b) {b.print("test_val");}
    template <typename T>
    void test_ref(T &b) {b.print("test_ref");}
    template <typename T>
    void test_ptr(T *pb) {pb->print("test_ptr");}
#else
    // -------- functions for base class --------
    void test_val(Base b) {b.print("test_val");}
    void test_ref(Base &b) {b.print("test_ref");}
    void test_ptr(Base *pb) {pb->print("test_ptr");}
    #ifdef OVERLOAD
        // -------- overloaded functions --------
        void test_val(Derived b) {b.print("test_val");}
        void test_ref(Derived &b) {b.print("test_ref");}
        void test_ptr(Derived *pb) {pb->print("test_ptr");}
    #endif
#endif

// -------- functions for derived class --------
void test_val_deriv(Derived b) {b.print("test_val");}
void test_ref_deriv(Derived &b) {b.print("test_ref");}
void test_ptr_deriv(Derived *pb) {pb->print("test_ptr");}

int main() {
    cout << "Initialization:" << endl;
    Derived d(1, 2);
    Base &ref_d1 = d;
    Derived &ref_d2 = d;
    Base *ptr_d1 = &d;
    Derived *ptr_d2 = &d;
    cout << endl;

    // Direct access
    cout << "Direct access:" << endl;
    d.print("direct_access");
    cout << endl;

    // Functions for derived type
    cout << "Function for derived type:" << endl;
    cout << "Pass by value:" << endl;
    test_val_deriv(d);
    cout << "Implicit reference:" << endl;
    test_ref_deriv(d);
    cout << "Explicit reference:" << endl;
    test_ref_deriv(ref_d2);
    cout << "Implicit pointer:" << endl;
    test_ptr_deriv(&d);
    cout << "Explicit pointer:" << endl;
    test_ptr_deriv(ptr_d2);
    cout << endl;

    // Pass-by-value
    cout << "Pass by value:" << endl;
    test_val(d);
    cout << endl;

    // Pass-by-reference
    cout << "Pass by reference:" << endl;
    cout << "Implicit:" << endl;
    test_ref(d);
    cout << "Explicit base:" << endl;
    test_ref(ref_d1);
    cout << "Explicit derived:" << endl;
    test_ref(ref_d2);
    cout << endl;

    // Pass-by-pointer
    cout << "Pass by pointer:" << endl;
    cout << "Implicit:" << endl;
    test_ptr(&d);
    cout << "Explicit base:" << endl;
    test_ptr(ptr_d1);
    cout << "Explicit derived:" << endl;
    test_ptr(ptr_d2);
    cout << endl;

    cout << "Finalization:" << endl;
}