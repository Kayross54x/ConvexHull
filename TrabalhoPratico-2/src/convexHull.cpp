#include "../include/convexHull.h"
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cassert>

//Calcula a distância entre dois pontos
int squaredDist(Dot p1, Dot p2) {
    return ((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
}

//Calcula a angulação entre dois pontos
double angulation(Dot first, Dot second) {
    if(second.x - first.x == 0) {
        if(second.y - first.y > 0) {
            return M_PI/2;
        } else {
            return 3*M_PI/2;
        }
    }
    
    return atan2(second.y - first.y, second.x - first.x);
}

//Retorna o ponto anterior ao topo em uma pilha
Dot nextToTop(Stack<Dot> &S) {
    assert(!S.isEmpty() && "A pilha não possui elementos suficientes, está vazia!");

    if(S.size() == 1) {
        cout << "A pilha possui apenas um elemento" << endl;
        return S.top();
    } else {
        Dot p = S.top();
        S.pop();
        Dot res = S.top();
        S.push(p);
        return res;
    }
}

//Retorna o ponto de menor coordenada X
Dot minXPoint(Dot array[], int length) {
    assert((length >= 2) && "A quantidade de pontos é muito pequena!");

    int i = 0;
    Dot min = array[0];

    while(i < length) {
        if(array[i].x < min.x) {
            min = array[i];
        }
        i++;
    }

    return min;
}

//Retorna o ponto de menor coordenada Y e, em caso de empate, a que tiver a menor coordenada X
int minXYPoint(Dot array[], int length) {
    assert((length >= 2) && "A quantidade de pontos é muito pequena!");
    
    int ymin = array[0].y;
    int min = 0;

    for (int i = 1; i < length; i++) {
        int y = array[i].y;

        if ((y < ymin) || (ymin == y && array[i].x < array[min].x)) {
            ymin = array[i].y;
            min = i;
        }
    }

    return min;
}

//Função auxiliar necessária para o funcionamento do algoritmo de ordenação dos buckets da função bucketSort
bool comparePoints(Dot p1, Dot p2, Dot initial) {
    double angle1 = angulation(initial, p1);
    double angle2 = angulation(initial, p2);

    if (angle1 < angle2)
        return true;
    else if (angle1 > angle2)
        return false;
    else
        return squaredDist(p1, initial) < squaredDist(p2, initial);
}

//Função para verificar qual a orientação da reta que se forma pela ligação dos 3 pontos passados por parâmetro
int ConvexHull::orientation(Dot first, Dot second, Dot third) {
    double d = (third.y - second.y) * (second.x - first.x) - (second.y - first.y) * (third.x - second.x);

    if(d > 0) return 1; //sentido anti-horário
    if(d < 0) return -1; //sentido horário
    return 0; //colineares
}

//Função auxiliar para verificar qual o tipo de ordenação deverá ser feita
void ConvexHull::order(Dot array[], int length, string type, Dot initialDot) {
    assert((type == "mergeSort" || type == "insertionSort" || type == "bucketSort") && "Tipo de ordenação inválido!");

    if(type == "mergeSort") {
        mergeSort(array, 0, length - 1, initialDot);
    } else if(type == "insertionSort") {
        insertionSort(array, length, initialDot);
    } else if(type == "bucketSort") {
        bucketSort(array, length, initialDot);
    }
}

//Função auxiliar necessária para o funcionamento do algoritmo de ordenação mergeSort
void ConvexHull::merge(Dot array[], int const left, int const mid, int const right, Dot initialDot) {
    auto const subArrayOne = mid - left + 1;
    auto const subArrayTwo = right - mid;

    //Cria arrays temporários
    auto *leftArray = new Dot[subArrayOne],
         *rightArray = new Dot[subArrayTwo];

    //Copia os pontos para os dois arrays temporários 
    for (auto i = 0; i < subArrayOne; i++)
        leftArray[i] = array[left + i];
    for (auto j = 0; j < subArrayTwo; j++)
        rightArray[j] = array[mid + 1 + j];

    int indexOfSubArrayOne = 0;
    int indexOfSubArrayTwo = 0;

    int indexOfMergedArray = left;
 
    double angle1, angle2;

    //Realizando o merge dos dois arrays temporários de volta para o array original
    while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
        //O critério de ordenação é a inclinação dos pontos em relação ao ângulo inicial
        angle1 = angulation(initialDot, leftArray[indexOfSubArrayOne]);
        angle2 = angulation(initialDot, rightArray[indexOfSubArrayTwo]);

        if (angle1 < angle2) {
            array[indexOfMergedArray].x = leftArray[indexOfSubArrayOne].x;
            array[indexOfMergedArray].y = leftArray[indexOfSubArrayOne].y;
            
            indexOfSubArrayOne++;
        } else if (angle1 > angle2) {
            array[indexOfMergedArray].x = rightArray[indexOfSubArrayTwo].x;
            array[indexOfMergedArray].y = rightArray[indexOfSubArrayTwo].y;

            indexOfSubArrayTwo++;
        } else {  //Em caso de empate na inclinação, compare as distâncias ao ponto inicial e escolha o ponto que tiver a menor
            double dist1 = squaredDist(initialDot, leftArray[indexOfSubArrayOne]);
            double dist2 = squaredDist(initialDot, rightArray[indexOfSubArrayTwo]);

            if (dist1 < dist2) {  // Compare as distâncias usando '<' em vez de '<='
                array[indexOfMergedArray].x = leftArray[indexOfSubArrayOne].x;
                array[indexOfMergedArray].y = leftArray[indexOfSubArrayOne].y;

                indexOfSubArrayOne++;
            } else {
                array[indexOfMergedArray].x = rightArray[indexOfSubArrayTwo].x;
                array[indexOfMergedArray].y = rightArray[indexOfSubArrayTwo].y;

                indexOfSubArrayTwo++;
            }
        }
        indexOfMergedArray++;
    }

    //Copia os elementos restantes do vetor da esquerda se tiver sobrado algum
    while (indexOfSubArrayOne < subArrayOne) {
        array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
        indexOfSubArrayOne++;
        indexOfMergedArray++;
    }

    //Copia os elementos restantes do vetor da direita se tiver sobrado algum
    while (indexOfSubArrayTwo < subArrayTwo) {
        array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
        indexOfSubArrayTwo++;
        indexOfMergedArray++;
    }

    delete[] leftArray;
    delete[] rightArray;
}

void ConvexHull::mergeSort(Dot array[], int const begin, int const end, Dot initialDot) {
    if (begin >= end) {
        return;
    }

    auto mid = begin + (end - begin) / 2;
    mergeSort(array, begin, mid, initialDot);
    mergeSort(array, mid + 1, end, initialDot);
    merge(array, begin, mid, end, initialDot);
}

void ConvexHull::insertionSort(Dot array[], int length, Dot initialDot) {
    int i = 1;

    while(i < length) {
        int j = i;
        double angle1 = angulation(initialDot, array[j-1]);
        double angle2 = angulation(initialDot, array[j]);
        if (angle1 > angle2) {
            while(j > 0) {
                angle1 = angulation(initialDot, array[j-1]);
                angle2 = angulation(initialDot, array[j]);
                if(angle1 > angle2) {
                    Dot temp = array[j - 1];
                    array[j - 1] = array[j];
                    array[j] = temp;
                }
                j--;
            }
        } else if (angle1 == angle2) {
            double dist1 = squaredDist(initialDot,  array[j-1]);
            double dist2 = squaredDist(initialDot,  array[j]);

            if (dist1 > dist2) {
                while(j > 0) {
                    dist1 = squaredDist(initialDot, array[j-1]);
                    dist2 = squaredDist(initialDot, array[j]);
                    if(dist1 > dist2) {
                        Dot temp = array[j - 1];
                        array[j - 1] = array[j];
                        array[j] = temp;
                    }
                    j--;
                }
            }
        }
        i++;
    }
}

void ConvexHull::bucketSort(Dot arr[], int n, Dot initial) {
    // Definindo o número máximo de baldes (aqui, usaremos 10)
    const int numBuckets = 10;
    
    // Criando os baldes
    Dot* buckets[numBuckets];
    for (int i = 0; i < numBuckets; i++) {
        buckets[i] = nullptr;
    }
    
    // Calculando o tamanho de cada balde
    int bucketSizes[numBuckets] = {0};
    for (int i = 0; i < n; i++) {
        int bucketIndex = numBuckets * angulation(initial, arr[i]) / (2 * M_PI);
        bucketSizes[bucketIndex]++;
    }
    
    // Alocando memória para os baldes
    for (int i = 0; i < numBuckets; i++) {
        buckets[i] = new Dot[bucketSizes[i]];
    }
    
    // Colocando os pontos nos baldes de acordo com suas angulações
    int bucketCounters[numBuckets] = {0};
    for (int i = 0; i < n; i++) {
        int bucketIndex = numBuckets * angulation(initial, arr[i]) / (2 * M_PI);
        buckets[bucketIndex][bucketCounters[bucketIndex]] = arr[i];
        bucketCounters[bucketIndex]++;
    }
    
    // Ordenando os pontos em cada balde usando std::sort
    for (int i = 0; i < numBuckets; i++) {
        std::sort(buckets[i], buckets[i] + bucketSizes[i], [initial](Dot p1, Dot p2) { return comparePoints(p1, p2, initial); });
    }
    
    // Juntando os pontos ordenados de todos os baldes
    int index = 0;
    for (int i = 0; i < numBuckets; i++) {
        for (int j = 0; j < bucketSizes[i]; j++) {
            arr[index++] = buckets[i][j];
        }
        delete[] buckets[i];
    }

    //Algoritmo de complexidade O(n)
}

Stack<Dot> ConvexHull::jarvisMarch(Dot array[], int length) {
    Dot on_hull = minXPoint(array, length); //O(n)

    Dot initialPoint = on_hull;
    Stack<Dot> resultPoints;
    while (1) { //O(h) numero de pontos no resultPoints
        resultPoints.push(on_hull);
        Dot nextPoint = array[0];
        for (int j = 0; j < length; j++) { //O(n)
            int o = orientation(on_hull, nextPoint, array[j]);
            if ((nextPoint.x == on_hull.x && nextPoint.y == on_hull.y) || o == 1) {
       
                nextPoint = array[j];
            }
        }

        on_hull = nextPoint;
        if ((on_hull.x == initialPoint.x) && (on_hull.y == initialPoint.y)){
            break;
        }
        
    }

    //T(n) = nh + n = O(nh) ---- O problema dessa solução é que quando o número de pontos no resultPoints é proximo de n, o tempo de execução aumenta muito e fica O(n²)
    return resultPoints;
}

Stack<Dot> ConvexHull::grahamScan(Dot array[], int length, string orderType) {
    Stack<Dot> resultPoints;
    int minPointIndex = minXYPoint(array, length); //O(n)

    //Colocando o ponto minimo na primeira posição
    Dot temp = array[0];
    array[0] = array[minPointIndex];
    array[minPointIndex] = temp;

    Dot p0 = array[0];

    //Ordenando os pontos restantes de acordo com a angulação em relação ao ponto minimo
    order(&array[1], length - 1, orderType, p0);

    //Removendo os pontos colineares caso existam
    int m = 1;
    for (int i = 1; i < length; i++) {
        while (i < length-1 && orientation(p0, array[i], array[i+1]) == 0) {
            i++;
        }
        array[m] = array[i];
        m++;
    }

    if (m < 3) return resultPoints;

    //Criando a pilha e colocando os 3 primeiros pontos
    resultPoints.push(array[0]);
    resultPoints.push(array[1]);
    resultPoints.push(array[2]);
   
    for (int i = 3; i < m; i++) {
        //Continuo removendo o topo da pilha até que o ângulo formado pelos pontos nextToTop, top e points[i] forme uma curva que não seja no sentido anti-horário        // Keep removing top while the angle formed by
        while (orientation(nextToTop(resultPoints), resultPoints.top(), array[i]) != 1) {
            resultPoints.pop();
        }
        resultPoints.push(array[i]);
    }

    //Complexidade de tempo no pior caso: O(N log N)
    //Complexidade de tempo no caso médio: O(N log N)
    //Complexidade de espaço: O(N)
   
   return resultPoints;
}

//-----------------------------------TAD PILHA-----------------------------------
template <typename T>
Stack<T>::Stack() {
    topo = -1;
    length = 0;
}

template <typename T>
bool Stack<T>::isEmpty() {
    return topo == -1;
}

template <typename T>
bool Stack<T>::isFull() {
    return topo == 9999;
}

template <typename T>
void Stack<T>::push(T x) {
    topo++;
    elementos_[topo] = x;
}

template <typename T>
int Stack<T>::size(){
    return topo+1;
}

template <typename T>
T Stack<T>::pop() {
    T valor = elementos_[topo];
    topo--;
    return valor;
}

template <typename T>
T Stack<T>::top() {
    return elementos_[topo];
}