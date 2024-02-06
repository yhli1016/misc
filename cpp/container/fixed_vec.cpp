#include "base.h"

template<typename T>
class Vector {
    private:
        T *ptr = nullptr;
        int m_size = 0;
    private:
        void log(const char* msg) const {
            cout << "\033[32m" << this << "\033[0m ";
            cout << "\033[34m" << msg << "\033[0m ";
            cout << endl;
        }
    public:
        Vector(int size = 0) {
            if (size > 0) {
                ptr = new T[size];
                m_size = size;
            }
        }
        ~Vector() {
            free();
        }
        Vector(const Vector &other) {
            copy(other);
            log("Vector copy constructor");
        }
        Vector & operator = (const Vector &other) {
            copy(other);
            log("Vector copy assignment operator");
            return *this;
        }
        Vector(Vector &&other) {
            move(other);
            log("Vector move constructor");
        }
        Vector & operator = (Vector &&other) {
            move(other);
            log("Vector move assignment operator");
            if (this == &other)
                return *this;
        }
        T & operator [] (int i) {
            return ptr[i];
        }
        void free() {
            if (ptr != nullptr && m_size > 0) {
                delete[] ptr;
                m_size = 0;
            }
        }
        void reset() {
            ptr = nullptr;
            m_size = 0;
        }
        void copy(Vector &other) {
            if (other.ptr != nullptr && other.m_size > 0) {
                free();
                ptr = new T[other.m_size];
                m_size = other.m_size;
                for (int i = 0; i < m_size; ++i) {
                    ptr[i] = other.ptr[i];
                }
            }
        }
        void move(Vector &other) {
            if (other.ptr != nullptr && other.m_size > 0) {
                free();
                ptr = other.ptr;
                m_size = other.m_size;
                other.reset();
            }
        }
};

void vector_val(int n, Vector<Base> &vec) {
    for(int i = 0; i < n; i++){
        cout << "i = " << i << endl;
        vec[i] = Base(i);
    }
}

int main() {
    // first approach
    // Vector<Base> vec(10);

    // second approach
    Vector<Base> vec;
    vec = Vector<Base>(10);

    vector_val(10, vec);
    cout << "Destruction:" << endl;
}
