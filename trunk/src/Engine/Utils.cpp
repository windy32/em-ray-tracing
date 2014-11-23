#define WIN32_LEAN_AND_MEAN

#include <windows.h>
#include <psapi.h>
#include <stdio.h>
#include "Utils.h"

int Utils::startTime;

void Utils::StartTimer()
{
    startTime = GetTickCount();
}

int Utils::GetTickCount()
{
    return (int)::GetTickCount();
}

void Utils::DbgPrint(const char *format, ...)
{
    char buf[1024];
    va_list argList;
    va_start(argList, format);

    vsprintf_s(buf, sizeof(buf), format, argList);
    fprintf(stderr, "%s", buf);

    va_end(argList);
}

void Utils::PrintTickCount(const char *desc)
{
    int tickCount = GetTickCount();
    fprintf(stderr, "%s: %.2lf\n", desc, tickCount / 1000.0);
}

void Utils::PrintTime(const char *desc)
{
    double time = (GetTickCount() - startTime) / 1000.0;
    fprintf(stderr, "%s: %.2lf\n", desc, time);
}

int Utils::GetMemorySize()
{
    HANDLE hProcess = GetCurrentProcess();  
    PROCESS_MEMORY_COUNTERS pmc;  
    GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc));
    return (int)pmc.PagefileUsage;
}
