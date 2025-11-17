//
// Created by Goswami, Sayan on 11/17/25.
//

#ifndef METAGRAPHRF_PRELUDE_H
#define METAGRAPHRF_PRELUDE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <unistd.h>
#include <string>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define RESET "\x1B[0m"

#define HIGH BLU
#define MED  MAG
#define LOW  CYN
#define LOGLVL LOW

static char *time_str(){
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    static char buf[9];
    strftime (buf, 9,"%T",timeinfo);
    return buf;
}

#define log_error(fmt, ...) do { \
fprintf(stderr, "[%s]" RED "ERROR: " fmt RESET " at %s:%i\n",       \
time_str(), ##__VA_ARGS__, __FILE__, __LINE__);     \
exit(1);                                                                      \
} while(0)

#define log_info(fmt, ...) do { \
fprintf(stderr, "[%s]" GRN "INFO: " fmt RESET "\n", time_str(), ##__VA_ARGS__); \
} while(0)

#define log_debug(lvl, fmt, ...) do { if (lvl[3] <= LOGLVL[3]) { \
fprintf(stderr, "[%s]" lvl "INFO: " fmt RESET "\n", time_str(), ##__VA_ARGS__); \
}} while(0)

#define sitrep(fmt, ...) do { \
fprintf(stderr, "\r" MAG "STATUS: " fmt RESET "", ##__VA_ARGS__); \
fflush(stderr); \
} while(0)

#define stderrflush fprintf(stderr, "\n")

#define log_warn(fmt, ...) do { \
fprintf(stderr, "[%s]" YEL "WARN: " fmt RESET " at %s:%i\n",       \
time_str(), ##__VA_ARGS__, __FILE__, __LINE__);     \
} while(0)

#endif //METAGRAPHRF_PRELUDE_H