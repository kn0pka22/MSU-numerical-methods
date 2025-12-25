#include <iostream>
#include <cstring>   
#include <cstdlib>   
#include <ctime> 

const int bufferSize = 1024 * 1024 * 1024; 

int main(){
    char* memoryBlock = (char*)malloc(bufferSize * sizeof(char));
    
    if (memoryBlock == nullptr) {
        std::cout << "Memory allocation failed!" << std::endl;
        return 1;
    }

    // Test memset
    clock_t startTime = clock();
    memset(memoryBlock, ' ', bufferSize);
    clock_t endTime = clock();
    double memsetDuration = double(endTime - startTime) / CLOCKS_PER_SEC;
    
    std::cout << "Memset time: " << memsetDuration << " seconds" << std::endl;

    startTime = clock();
    for (int currentIndex = 0; currentIndex < bufferSize; currentIndex++) {
        memoryBlock[currentIndex] = ' ';
    }
    endTime = clock();
    double indexedLoopDuration = double(endTime - startTime) / CLOCKS_PER_SEC;
    
    std::cout << "по 1му заполнение: " << indexedLoopDuration << " seconds" << std::endl;

    // startTime = clock();
    // unsigned long long int a;
    // unsigned long long int* b;
    // b = 

    // unsigned char* x;
    // for ()

    // char* currentPosition = memoryBlock;
    // char* endPosition = memoryBlock + bufferSize;
    // while (currentPosition < endPosition) {
    //     *currentPosition++ = ' ';
    // }
    // endTime = clock();
    // double pointerLoopDuration = double(endTime - startTime) / CLOCKS_PER_SEC;
    
    // std::cout << "Pointer loop time: " << pointerLoopDuration << " seconds" << std::endl;

    startTime = clock();
    for (int currentIndex = 0; currentIndex < bufferSize; currentIndex += 8) {
        memoryBlock[currentIndex] = ' ';
        if (currentIndex + 1 < bufferSize) memoryBlock[currentIndex + 1] = ' ';
        if (currentIndex + 2 < bufferSize) memoryBlock[currentIndex + 2] = ' ';
        if (currentIndex + 3 < bufferSize) memoryBlock[currentIndex + 3] = ' ';
        if (currentIndex + 4 < bufferSize) memoryBlock[currentIndex + 4] = ' ';
        if (currentIndex + 5 < bufferSize) memoryBlock[currentIndex + 5] = ' ';
        if (currentIndex + 6 < bufferSize) memoryBlock[currentIndex + 6] = ' ';
        if (currentIndex + 7 < bufferSize) memoryBlock[currentIndex + 7] = ' ';
    }
    endTime = clock();
    double unrolledLoopDuration = double(endTime - startTime) / CLOCKS_PER_SEC;
    
    std::cout << "по 8 заполнение: " << unrolledLoopDuration << " seconds" << std::endl;

    free(memoryBlock);

    std::cout << "\nPerformance comparison:" << std::endl;
    std::cout << "Memset is " << indexedLoopDuration / memsetDuration 
              << " times faster than indexed" << std::endl;
    // std::cout << "Memset is " << pointerLoopDuration / memsetDuration 
    //           << " times faster than pointer loop" << std::endl;
    std::cout << "Memset is " << unrolledLoopDuration / memsetDuration 
              << " times faster than 8indexed" << std::endl;

    return 0;
}