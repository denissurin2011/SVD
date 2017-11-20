#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace std;

const long double eps = 1e-8;

vector<vector<long double>> u;
vector<vector<long double>> v;

long double sign(long double a) {
    if (a > 0) {
        return 1;
    }
    if (a < 0) {
        return -1;
    }
    return 0;
}

long double sum_col_doubled (vector<vector<long double>> u, int j, int matrix_size) {
    long double res = 0;
    for (int k = 0; k < matrix_size; k++) {
        res += u[k][j] * u[k][j];
    }
    return res;
}

long double scal(vector<vector<long double>> u, int i, int j, int matrix_size) {
    long double res = 0;
    for (int k = 0; k < matrix_size; k++) {
        res += u[k][i] * u[k][j];
    }
    return res;
}

void rotate_u(int i, int j, long double cos, long double sin, int matrix_size) {
    vector <long double> u_old(matrix_size);
    for (int k = 0; k < matrix_size; k++) {
        u_old[k] = u[k][j];
    }
    for (int k = 0; k < matrix_size; k++) {
        u[k][j] = u_old[k] * cos - sin * u[k][i];
    }
    for (int k = 0; k < matrix_size; k++) {
        long double f = u[k][i];
        u[k][i] = u_old[k] * sin + cos * f;
    }
}

void rotate_v(int i, int j, long double cos, long double sin, int matrix_size) {
    vector <long double> v_old(matrix_size);
    for (int k = 0; k < matrix_size; k++) {
        v_old[k] = v[k][j];
    }
    for (int k = 0; k < matrix_size; k++) {
        v[k][j] = v_old[k] * cos - sin * v[k][i];
    }
    for (int k = 0; k < matrix_size; k++) {
        v[k][i] = v_old[k] * sin + cos * v[k][i];
    }
}

void svd(int matrix_size) {
    srand(time(NULL));
    u.resize(matrix_size);
    for (int i = 0; i < matrix_size; i++) {
        u[i].resize(matrix_size);
        for (int j = 0; j < matrix_size; j++) {
            u[i][j] = rand() % (matrix_size * 100);
        }
    }
    cout << "Matrix" << endl;
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            cout << u[i][j] << " ";
        }
        cout << '\n';
    }
    cout << endl;
    v.resize(matrix_size);
    for (int i = 0; i < matrix_size; i++) {
        v[i].resize(matrix_size);
        for (int j = 0; j < matrix_size; j++) {
            if (i == j) {
                v[i][j] = 1;
            } else {
                v[i][j] = 0;
            }
        }
    }
    for (int y = 0; y < 100; y++) {
        long double count = 0;
        for (int i = 0; i < matrix_size - 1; i++) {
            for (int j = i + 1; j < matrix_size; j++) {
                long double a = sum_col_doubled(u, j, matrix_size);
                long double b = sum_col_doubled(u, i, matrix_size);
                long double c = scal(u, i, j, matrix_size);
                long double d = (b - a) / (2 * c);
                long double t = sign(d) / (abs(d) + sqrt(1 + d * d));
                long double cos = 1 / sqrt(1 + t * t);
                long double sin = t / sqrt(1 + t * t);
                count = max(count, abs(c) / sqrt(a * b));
                rotate_u(i, j, cos, sin, matrix_size);
                rotate_v(i, j, cos, sin, matrix_size);
            }
        }
        if (count < eps) {
            break;
        }
    }
    vector<long double> signul(matrix_size, 0);
    for (int i = 0; i < matrix_size; i++) {
        long double res = 0;
        for (int j = 0; j < matrix_size; j++) {
            res += (u[j][i]) * (u[j][i]);
        }
        signul[i] = sqrt(res);
        for (int j = 0; j < matrix_size; j++) {
            u[j][i] /= signul[i];
        }
    }
    cout << "U:" << endl;
    for (int i = 0; i < u.size(); i++) {
        for (int j = 0; j < u[i].size(); j++) {
            cout << u[i][j] << " ";
        }
        cout << '\n';
    }
    cout << endl;
    cout << "E:" << endl;
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            if (i == j) {
                cout << signul[i] << " ";
                continue;
            }
            cout << 0 << " ";
        }
        cout << '\n';
    }
    cout << endl;
    cout << "V:" << endl;
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            cout << v[i][j] << " ";
        }
        cout << '\n';
    }
}

int main() {
    int matrix_size;
    cin >> matrix_size;
    svd(matrix_size); // for demonstration svd generates random matrix.
    return 0;
}