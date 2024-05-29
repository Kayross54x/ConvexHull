#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

template <typename T>
class Stack {
    public:
        Stack();

        bool isEmpty();

        bool isFull();

        void push(T x);

        T pop();

        T top() ;

        int size();

    private:
        int topo;
        T elementos_[10000];
        double length;
};