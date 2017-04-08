#pragma once

#ifndef LOG_H
#define LOG_H

#include "stdio.h"
#include "time.h"

#ifdef _WIN
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64 

#endif // _WIN


#ifdef _WIN

// MSVC defines this in winsock2.h!?
/*typedef struct timeval {
	long tv_sec;
	long tv_usec;
} timeval;
*/

int gettimeofday(struct timeval * tp, struct timezone * tzp)
{
	// Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
	// This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
	// until 00:00:00 January 1, 1970 
	static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

	SYSTEMTIME  system_time;
	FILETIME    file_time;
	uint64_t    time;

	GetSystemTime(&system_time);
	SystemTimeToFileTime(&system_time, &file_time);
	time = ((uint64_t)file_time.dwLowDateTime);
	time += ((uint64_t)file_time.dwHighDateTime) << 32;

	tp->tv_sec = (long)((time - EPOCH) / 10000000L);
	tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
	return 0;
}

#endif // _WIN

/// printf trace to output with current time in milliseconds
void lprintf(const char *text)
{
	char buf[1024];
	struct timeval tp;
	char buf1[1024];
#ifdef _WIN
	struct timezone *tz = NULL;
	gettimeofday(&tp, tz);
	//gettimeofday(&tp, NULL);
#else
	struct timezone tz;
	gettimeofday(&tp, &tz);
#endif
	struct tm tm0;
	time_t lt = time(NULL);
#ifndef _WIN
	localtime_r(&lt, &tm0);
#else
	localtime_s(&tm0, &lt);
#endif _WIN
	strftime(buf, sizeof(buf), "%Y/%m/%d %H:%M:%S:", &tm0);
	size_t sizestr = strlen(buf);
	sprintf_s(buf + sizestr, 1024 - sizestr, "%03u: ", (unsigned)tp.tv_usec / 1000);
	sizestr = strlen(buf);
	sprintf_s(buf + sizestr, 1024 - sizestr, " %s\n", text);
	printf("%s", buf);
}

/*void test1(int &a)
{
	return;
}
*/

#endif // LOG_H 