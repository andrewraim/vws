#ifndef VWS_LOGGER_H
#define VWS_LOGGER_H

#include <Rcpp.h>

namespace vws {

inline void logger(const char* fmt, ...)
{
	// Get the current time in the local time zone
	std::time_t raw = std::time(nullptr);
	std::tm* local = std::localtime(&raw);

	// Write formatted time to a string
	char buffer[64];
	strftime(buffer, 64, "%Y-%m-%d %H:%M:%S", local);

	// Insert "..." placeholders into format string for message; see
	// <https://stackoverflow.com/q/1056411>
	char msg[256];
	va_list args;
	va_start(args, fmt);
	vsnprintf(msg, 255, fmt, args);
	va_end(args);

	// Print the formatted message with timestamp
	Rprintf("%s - %s", buffer, msg);
}

}

#endif
