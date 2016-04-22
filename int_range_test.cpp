#include <iostream>

using namespace std;

int main() {
    int a;
    for (a=0; a<=655360; a++) {
        if (a%10000==0){
            cout << a << endl;
        }
    }
    return 0;
}
