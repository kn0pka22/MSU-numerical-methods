#include <iostream>
#include <cstdlib>

int main() {
    const int TOTAL = 10000;
    int asc_count = 0;
    int desc_count = 0;

    std::cout << "Размер unsigned int: " << sizeof(unsigned int) << " байт" << std::endl;
    std::cout << "Размер unsigned long long: " << sizeof(unsigned long long) << " байт" << std::endl;
    std::cout << std::endl;

    unsigned int a = rand();
    unsigned int b = rand();

    for (int i = 0; i < TOTAL - 2; i++) {
        unsigned long long sum = a;
        sum += b;
        unsigned int c = sum & 0xFFFFFFFF; // Оставляем только младшие 32 бита

        // Проверка условий возрастания и убывания
        if (a < b && b < c) {
            asc_count++;
        } else if (c < b && b < a) {
            desc_count++;
        }
        // Сдвиг значений для следующей итерации
        a = b;
        b = c;
    }

    std::cout << "Результаты:" << std::endl;
    std::cout << "Всего троек: " << TOTAL - 2 << std::endl;
    std::cout << "Возрастающих: " << asc_count << std::endl;
    std::cout << "Убывающих: " << desc_count << std::endl;
    std::cout << "Доля возрастающих: " << (asc_count * 100. / (TOTAL - 2)) << "%" << std::endl;
    std::cout << "Доля убывающих: " << (desc_count * 100. / (TOTAL - 2)) << "%" << std::endl;

    return 0;
}