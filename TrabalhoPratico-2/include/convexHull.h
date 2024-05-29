#include "./stack.h"
#include <string>

using namespace std;

struct Dot {
    int x;
    int y;
};

class ConvexHull {
    public:
        Stack<Dot> jarvisMarch(Dot array[], int length);

        Stack<Dot> grahamScan(Dot array[], int length, string orderType);

    private:
        int orientation(Dot first, Dot second, Dot third);

        void order(Dot array[], int length, string type, Dot initialDot);        

        void mergeSort(Dot array[], int const begin, int const end, Dot initialDot);

        void merge(Dot array[], int const begin, int const mid, int const end, Dot initialDot);

        void insertionSort(Dot array[], int length, Dot initialDot);

        void bucketSort(Dot array[], int length, Dot initialDot);
};