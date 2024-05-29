#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include "./convexHull.cpp"
#include <time.h>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Falha: Arquivo não encontrado na passagem da linha de comando" << std::endl;
        return 1;
    }

    ifstream input_file(argv[1]);

    if (!input_file.is_open()) {
        cerr << "Erro ao abrir o arquivo de entrada!" << endl;
        cerr << "Verifique se o arquivo de entrada informado na linha de comando está no mesmo diretório deste arquivo makefile" << endl;
    } else {
        try {
            ConvexHull solve;
            Dot Dots[100000] = {};
            string line;

            int length = 0;

            while (getline(input_file, line)) {
                istringstream iss(line);

                Dot ponto;
                int i = 0;
                int x = 0;
                int y = 0;
                while(iss >> line) {
                    if(i == 0) {
                        x = std::stoi(line);
                    } else {
                        y = std::stoi(line);
                        break;
                    }
                    i++;
                }
                ponto.x = x;
                ponto.y = y;
                Dots[length] = ponto;
                length++;
            };

            Dot aux0[length] = {};
            Dot aux1[length] = {};
            Dot aux2[length] = {};
            Dot aux3[length] = {};
            Dot aux4[length] = {};
            for(int i = 0; i < length; i++) {
                aux0[i] = Dots[i];
                aux1[i] = Dots[i];
                aux2[i] = Dots[i];
                aux3[i] = Dots[i];
                aux4[i] = Dots[i];
            }

            Stack<Dot> result = solve.grahamScan(aux0, length, "mergeSort");
            Stack<Dot> convexLines;

            std::ofstream arquivo1("line_data.txt");
            std::ofstream arquivo2("scatter_data.txt");

            cout << "FECHO CONVEXO: " << endl;
            while(!result.isEmpty()) {
                Dot point = result.pop();
                convexLines.push(point);
                if(arquivo1.is_open()) arquivo1 << point.x << " " << point.y << endl;
                if(arquivo2.is_open()) arquivo2 << point.x << " " << point.y << endl;
                cout << point.x << " " << point.y << endl;
            }

            cout << endl;

            clock_t start, end;
            double time;

            start = clock();
            solve.grahamScan(aux1, length, "mergeSort");
            end = clock();
            time = ((double) (end - start)) / CLOCKS_PER_SEC;
            cout << "GRAHAM+MERGESORT: " << time << "s" << endl;

            start = clock();
            solve.grahamScan(aux2, length, "insertionSort");
            end = clock();
            time = ((double) (end - start)) / CLOCKS_PER_SEC;
            cout << "GRAHAM+INSERTIONSORT: " << time << "s" << endl;

            start = clock();
            solve.grahamScan(aux3, length, "bucketSort");
            end = clock();
            time = ((double) (end - start)) / CLOCKS_PER_SEC;
            cout << "GRAHAM+LINEAR: " << time << "s" << endl;

            start = clock();
            solve.jarvisMarch(aux4, length);
            end = clock();
            time = ((double) (end - start)) / CLOCKS_PER_SEC;
            cout << "JARVIS: " << time << "s" << endl << endl;

            arquivo1.close();
            arquivo2.close();

            cout << "RETAS: " << endl;
            
            Dot initialPoint = convexLines.pop();
            Dot lastPoint = initialPoint;

            while(convexLines.size() >= 1) {
                Dot nextToTop = convexLines.pop();   

                if(lastPoint.x == nextToTop.x) {
                    cout << "x = " << lastPoint.x << endl;
                } else {
                    double a = (double)(lastPoint.y - nextToTop.y) / (lastPoint.x - nextToTop.x);
                    double b = lastPoint.y - a * lastPoint.x;
                    string condition = b > 0 ? "+ " : "- ";
                    b = b < 0 ? b * -1 : b;

                    cout << a << " x " << condition << b << " = 0" << endl;
                }

                lastPoint = nextToTop;
            }

            if(lastPoint.x == initialPoint.x) {
                cout << "x = " << lastPoint.x << endl;
            } else {
                double a = (double)(lastPoint.y - initialPoint.y) / (lastPoint.x - initialPoint.x);
                double b = lastPoint.y - a * lastPoint.x;
                string condition = b > 0 ? "+ " : "- ";
                b = b < 0 ? b * -1 : b;

                cout << a << " x " << condition << b << " = 0" << endl;
            }
        } catch (char const* e) {
            cout << e << endl;
        }
        
        input_file.close();
    }

    return 0;
}