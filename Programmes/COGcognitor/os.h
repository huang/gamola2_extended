#ifndef OS_H
#define OS_H

#include <unistd.h>

#if defined(__GNUG__)
	
	#define INT64 long long
	#define UINT64 unsigned long long
	#include <dirent.h>
	#define SPRTR '/'
	void myOpenDir(const char* dirpath)
	{
		chdir(dirpath);
	}

	#elif defined(_MSC_VER)
		#define INT64 __int64
		#define UINT64 unsigned __int64
		#include <direct.h>
		#define SPRTR '\\'
		void myOpenDir(const char* dirpath)
		{
			_chdir(dirpath);
		}
			#else
			#error "Platform not supported. Need to update source code"
			#endif 
#endif
