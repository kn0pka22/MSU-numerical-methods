#include "linearSystemSolver.hpp"

void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            std::cout << std::setw(10) << std::setprecision(4) << val << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const std::vector<double>& vec) {
    for (double val : vec) {
        std::cout << std::setw(10) << std::setprecision(4) << val << " ";
    }
    std::cout << std::endl;
}

double matrixNorm(const std::vector<std::vector<double>>& a) {
    int n = a.size(); 
    double max_norm = 0.0; 
    
    for (int i = 0; i < n; i++) {
        double sum = 0.0; 

        for (int j = 0; j < n; j++) {
            sum += fabs(a[j][i]); 
        }
        max_norm = std::max(max_norm, sum);
    }

    return max_norm; 
}


int solve(std::vector<double>& M, std::vector<double>& b, std::vector<double>& x, std::vector<int>& memory) {
    int n = std::sqrt(M.size()); // Предполагаем, что M - это плоская матрица n x n

    // Проверка на корректность размера матрицы
    if (n * n != M.size()) {
        std::cerr << "Error: Matrix size is not valid." << std::endl;
        return -1;
    }

    // Инициализация памяти
    memory.resize(n); // Убедитесь, что memory имеет правильный размер
    for (int j = 0; j < n; j++) {
        memory[j] = j; // Инициализация индексов
    }

    double mean = 1e-3;
    for (int i = 0; i < n * n; ++i) {
        mean += (M[i] > 0) ? M[i] : -M[i]; // Вычисление среднего значения
    }

    double eps = mean * 1e-10 / (n * n); // Параметр для проверки на сингулярность

    for (int step = 1; step <= n; step++) {
        int indMax = step;
        double max = M[e(step, step, n)];

        // Поиск главного элемента
        for (int i = indMax + 1; i <= n; ++i) {
            if ((M[e(step, i, n)] > 0 && max > 0 && M[e(step, i, n)] > max) ||
                (M[e(step, i, n)] > 0 && !(max > 0) && M[e(step, i, n)] > -max) ||
                (!(M[e(step, i, n)] > 0) && max > 0 && -M[e(step, i, n)] > max) ||
                (!(M[e(step, i, n)] > 0) && !(max > 0) && -M[e(step, i, n)] > -max)) {
                indMax = i;
                max = M[e(step, i, n)];
            }
        }

        // Проверка на сингулярность
        if ((max < eps && max > 0) || (max > -eps && !(max > 0))) {
            return -1; // Сингулярная матрица
        }

        // Меняем столбцы местами, если нужно
        if (indMax != step) {
            for (int i = 1; i <= n; ++i) {
                std::swap(M[e(i, step, n)], M[e(i, indMax, n)]);
            }
            std::swap(memory[indMax - 1], memory[step - 1]);
        }

        // Шаг прямого хода метода Гаусса
        for (int i = step + 1; i <= n; ++i) {
            double v = (M[e(i, step, n)] / max);
            for (int j = step; j <= n; ++j) {
                M[e(i, j, n)] -= v * M[e(step, j, n)];
            }
            b[i - 1] -= v * b[step - 1];
        }
    }

    // Обратный ход метода Гаусса
    for (int step = 1; step <= n; ++step) {
        for (int i = 1; i <= n - step; ++i) {
            b[i - 1] -= b[n - step] * M[e(i, n - step + 1, n)] / M[e(n - step + 1, n - step + 1, n)];
            M[e(i, n - step + 1, n)] = 0; // Обнуляем элемент
        }
        b[n - step] /= M[e(n - step + 1, n - step + 1, n)];
        M[e(n - step + 1, n - step + 1, n)] = 1; // Нормализуем элемент
    }

    // Сохранение результатов в векторе решения
    for (int u = 0; u < n; ++u) {
        x[memory[u]] = b[u];
    }

    return 0; // Успешное завершение
}

