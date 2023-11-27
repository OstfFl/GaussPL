#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

void readInput(double**& A, double*& b, int& matrix_size, string file);

void display(double**& A, double*& b, int& matrix_size);

void task1(string file);

void task4(string file);

void task2(string file);

void task3(string file);

void generalGaussElimination(double**& A, double*& b, int& matrix_size);

double* gaussXYZ(double**& A, double*& b, int& matrix_size);

void fragmentGaussElimination(double**& A, double*& b, int& matrix_size, int k);

void bubbleSort(vector<int>& columnIndices, double*& res, int& matrix_size);

void displayResults(int& matrix_size, double*& res)
{
    if (res == nullptr)
    {
        cout << "No results in this function!" << endl;
    }
    else
    {
        for (int i = 0; i < matrix_size; i++)
        {
            cout << "x" << i + 1 << "= " << setw(5) << res[i] << endl;
        }
        cout << endl;
        cout << "-----------------------------------------------------------------------------------------" << endl;
    }
}

int main()
{
    double** A;
    double* b;
    int matrix_size;
    readInput(A, b, matrix_size, R"(zad13.csv)");

    if (matrix_size <= 1)
    {
        cout << "----Small matrix!----" << endl;
        return 0;
    }
    else
    {
        cout << "\n\nTASK 1: " << endl;
        task1("zad1.csv"); // You can also input file paths directly here.

        // Uncomment the following lines to execute other tasks.
         cout << "\n\nTASK 2: " << endl;
         task2("zad14.csv");
         cout << "\n\nTASK 3: " << endl;
         task3("zad14.csv");
    	 cout << "\n\nTASK 4: " << endl;
         task4("zad13.csv");
        return 0;
    }
}

void readInput(double**& A, double*& b, int& matrix_size, string file)
{
    ifstream source_file(file);
    if (!source_file.is_open())
    {
        cout << "The file has not been opened!" << endl;
        exit(0);
    }
    else
    {
        source_file >> matrix_size;

        A = new double* [matrix_size];
        for (int i = 0; i < matrix_size; i++)
        {
            A[i] = new double[matrix_size];
        }

        b = new double[matrix_size];

        char semicolon;
        for (unsigned i = 0; i < matrix_size + 1; i++)
            source_file >> semicolon;

        for (unsigned i = 0; i < matrix_size; i++)
        {
            for (unsigned j = 0; j < matrix_size; j++)
            {
                source_file >> A[i][j];
                source_file >> semicolon;
            }
            source_file >> semicolon;
            source_file >> b[i];
        }
        source_file.close();
    }
}

void display(double**& A, double*& b, int& matrix_size)
{
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << " Size = " << matrix_size << endl;
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            cout << setw(13) << A[i][j];
        }
        cout << setw(20) << b[i];
        cout << "\n";
    }
    cout << "-----------------------------------------------------------------------------------------" << endl;
}

void generalGaussElimination(double**& A, double*& b, int& matrix_size)
{
    for (int k = 0; k < matrix_size - 1; k++)
    {
        fragmentGaussElimination(A, b, matrix_size, k);
    }
}

void fragmentGaussElimination(double**& A, double*& b, int& matrix_size, int k)
{
    for (int i = k + 1; i < matrix_size; i++)
    {
        double c = A[i][k] / A[k][k];
        for (int j = k; j < matrix_size; j++)
        {
            A[i][j] = A[i][j] - c * A[k][j];
        }
        b[i] = b[i] - c * b[k];
    }
}

void bubbleSort(vector<int>& columnIndices, double*& res, int& matrix_size)
{
    int i = matrix_size - 1;
    while (i > 0)
    {
        int j = 0;
        int p = 1;
        while (j < i)
        {
            if (columnIndices[j] > columnIndices[j + 1])
            {
                swap(columnIndices[j], columnIndices[j + 1]);
                swap(res[j], res[j + 1]);
                p = 0;
            }
            j = j + 1;
        }
        if (p == 1)
        {
            break;
        }
        i = i - 1;
    }
}

double* gaussXYZ(double**& A, double*& b, int& matrix_size)
{
    auto* res = new double[matrix_size];
    int n = matrix_size;
    if (A[0][0] == 1.19209e-07)
    {
        cout << "A[0][0] is equal to zero!" << endl;
        return 0;
    }
    else
    {
        for (int i = n - 1; i >= 0; i--)
        {
            res[i] = b[i];
            for (int j = i + 1; j < n; j++)
            {
                res[i] = res[i] - A[i][j] * res[j];
            }
            res[i] = res[i] / A[i][i];
        }
        return res;
    }
}

void task1(string file)
{
    double** A;
    double* b;

    int matrix_size;
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << "READING AND DISPLAYING THE ARRAY" << endl;
    readInput(A, b, matrix_size, file);
    display(A, b, matrix_size);
    cout << "-----------------------------------------------------------------------------------------" << endl;

    generalGaussElimination(A, b, matrix_size);

    display(A, b, matrix_size);
    cout << "Results: " << endl;
    double* res = gaussXYZ(A, b, matrix_size);
    displayResults(matrix_size, res);
    delete[] A;
    delete[] b;
}

void task4(string file)
{
    double** A;
    double* b;
    int matrix_size;
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << "READING AND DISPLAYING THE ARRAY" << endl;
    readInput(A, b, matrix_size, file);
    display(A, b, matrix_size);
    cout << "-----------------------------------------------------------------------------------------" << endl;
    int l = 0;
    vector<int> columnIndices;
    columnIndices.reserve(matrix_size);
    for (int i = 0; i < matrix_size; i++)
    {
        columnIndices.push_back(i);
    }
    for (int k = 0; k < matrix_size - 1; k++)
    {
        int needi = l, needj = l;
        double max = INT_MIN;
        for (int i = l; i < matrix_size; i++)
        {
            for (int j = l; j < matrix_size; j++)
            {
                if (abs(A[i][j]) > max)
                {
                    max = fabs(A[i][j]);
                    needi = i;
                    needj = j;
                }
            }
        }
        for (int i = 0; i < matrix_size; i++)
        {
            swap(A[i][l], A[i][needj]);
        }
        swap(b[l], b[needi]);
        for (int i = 0; i < matrix_size; i++)
        {
            swap(A[l][i], A[needi][i]);
        }
        swap(columnIndices[l], columnIndices[needj]);
        cout << "Step" << endl;
        display(A, b, matrix_size);

        fragmentGaussElimination(A, b, matrix_size, k);
        cout << "Step after Gaussian elimination" << endl;
        display(A, b, matrix_size);

        l++;
    }
    cout << "Final array: " << endl;
    display(A, b, matrix_size);
    double* res = gaussXYZ(A, b, matrix_size);
    bubbleSort(columnIndices, res, matrix_size);
    cout << "Results: " << endl;
    displayResults(matrix_size, res);
}

void task3(string file)
{
    double** A;
    double* b;

    int matrix_size;
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << "READING AND DISPLAYING THE ARRAY" << endl;
    readInput(A, b, matrix_size, file);
    display(A, b, matrix_size);
    cout << "-----------------------------------------------------------------------------------------" << endl;
    int l = 0;
    for (int k = 0; k < matrix_size - 1; k++)
    {
        if (A[0][0] <= 1.19209e-07)
        {
            cout << "A[0][0] is equal to zero!" << endl;
            return;
        }
        int needi = l, needj = l;
        for (int i = l; i < matrix_size; i++)
        {
            double maxm = fabs(A[l][i]);
            for (int j = l + 1; j < matrix_size; j++)
            {
                if (fabs(A[j][i]) > maxm)
                {
                    maxm = A[j][i];
                    needj = j;
                    needi = i;
                }
            }
        }
        for (int i = 0; i < matrix_size; i++)
        {
            swap(A[l][i], A[needj][i]);
        }
        swap(b[l], b[needj]);
        fragmentGaussElimination(A, b, matrix_size, k);

        l++;
    }
    cout << "Final array: " << endl;
    display(A, b, matrix_size);
    double* res = gaussXYZ(A, b, matrix_size);
    cout << "Results: " << endl;
    displayResults(matrix_size, res);
}

void task2(string file)
{
    double** A;
    double* b;
    int matrix_size;
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << "READING AND DISPLAYING THE ARRAY" << endl;
    readInput(A, b, matrix_size, file);
    display(A, b, matrix_size);
    cout << "-----------------------------------------------------------------------------------------" << endl;
    int l = 0;
    vector<int> columnIndices;
    columnIndices.reserve(matrix_size);
    for (int i = 0; i < matrix_size; i++)
    {
        columnIndices.push_back(i);
    }
    int needi = l, needj = l;
    for (int k = 0; k < matrix_size - 1; k++)
    {
        double max = INT_MIN;
        for (int i = l; i < matrix_size; i++)
        {
            for (int j = l; j < matrix_size; j++)
            {
                if (abs(A[i][j]) > max)
                {
                    max = fabs(A[i][j]);
                    needi = i;
                    needj = j;
                }
            }
        }
        for (int i = 0; i < matrix_size; i++)
        {
            swap(A[i][i], A[i][needj]);
        }
        swap(columnIndices[l], columnIndices[needj]);

        l++;
    }

    display(A, b, matrix_size);
    cout << "Final array: " << endl;
    display(A, b, matrix_size);
    double* res = gaussXYZ(A, b, matrix_size);
    bubbleSort(columnIndices, res, matrix_size);
    cout << "Results: " << endl;
    displayResults(matrix_size, res);
}
