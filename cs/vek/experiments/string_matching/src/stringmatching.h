#ifndef STRINGMATHCING_H
#define STRINGMATCHING_H

extern const char* naive(const char*, int, const char*, int);
int sse_strcmpeq(volatile const char* needle, volatile const char* haystack);

#endif