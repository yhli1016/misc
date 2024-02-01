#ifndef H_BASE_INCLUDED
#define H_BASE_INCLUDED

#include <iostream>
using std::cout;
using std::endl;

class Base {
    private:
        int a = 0;
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
        Base & operator = (const Base &other) {
            a = other.a;
            log("Copy assignment operator");
            return *this;
        }
#ifdef MOVE
        Base(Base &&other) {
            a = other.a;
            log("Move constructor");
        }
        Base & operator = (Base &&other) {
            a = other.a;
            log("Move assignment operator");
            if (this == &other)
                return *this;
        }
#endif
        ~Base() {
            log("Destructor");
        }
        void log(const char* msg) const {
            cout << "\033[32m" << this << "\033[0m ";
            cout << "\033[33m" << a << "\033[0m ";
            cout << "\033[34m" << msg << "\033[0m ";
            cout << endl;
        }
};

#endif