#ifndef UTILS_H
#define UTILS_H

class Utils
{
private:
    static int startTime;

public:
    // Added here to avoid include <windows.h>
    static void StartTimer();
    static int GetTickCount();

    // Debug output
    static void DbgPrint(const char *format, ...);
    static void PrintTickCount(const char *desc);
    static void PrintTime(const char *desc);

    // Memory
    static int GetMemorySize();
};

#endif
