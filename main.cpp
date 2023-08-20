#include "GKS1D.h"
#include "GKS2D.h"

int main()
{
    clock_t start, finish;
    start = clock();
    // Accuracy_test_1d();
    // Riemann_problem_1d();
    // Accuracy_test_2d();
    Riemann_problem_1d();
    finish = clock();
    cout << endl << "the time cost is:" << double(finish - start) / CLOCKS_PER_SEC << endl;
    return 0;
}
