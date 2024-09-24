#include "base.h"
#include <memory>

Base return_val(int i) {
    Base b(i);
    b.log("return_val");
    return b;
}

class BaseMap {
    private:
#ifdef SMART_PTR
        const std::shared_ptr<Base> base_ = nullptr;
    public:
        BaseMap(const Base& base): base_(std::make_shared<Base>(base)) {
            cout << "Constructore of map" << endl;
        }
#else
        const Base* base_ = nullptr;
    public:
        BaseMap(const Base& base): base_(&base) {
            cout << "Constructore of map" << endl;
        }
#endif
        const Base& return_ref() const {
            base_->log("return_ref");
            return *base_;
        }
        /*const Base* return_ptr() const {
            base_->log("return_ptr");
            return base_;
        }*/
};

int main() {
    cout << "Return by reference:" << endl;
    Base b1 = return_val(1);
    BaseMap map(b1);
    const Base &b2 = map.return_ref();
    b2.log("access_ref");
    cout << endl;
    cout << "Destruction" << endl;
}
