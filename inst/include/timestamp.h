#ifndef VWS_TIMESTAMP_H
#define VWS_TIMESTAMP_H

#include <Rcpp.h>

namespace vws {

/*
* Logger
*
* Get timestamp for current time as character array.
*/
inline std::string timestamp(void)
{
	time_t timer;
	char buffer[64];
	struct tm* tm_info;

	timer = time(NULL);
	tm_info = localtime(&timer);

	strftime(buffer, 64, "%Y-%m-%d %H:%M:%S", tm_info);
	return buffer;
}

}

#endif

