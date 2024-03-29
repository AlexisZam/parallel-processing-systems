#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MT_CONF "MT_CONF"

void setaffinity_oncpu(unsigned int cpu) {
    cpu_set_t cpu_mask;
    int err;

    CPU_ZERO(&cpu_mask);
    CPU_SET(cpu, &cpu_mask);

    err = sched_setaffinity(0, sizeof(cpu_set_t), &cpu_mask);
    if (err) {
        perror("sched_setaffinity");
        exit(1);
    }
}

static long parse_int(char *s) {
    long ret;
    char *endptr;

    ret = strtol(s, &endptr, 10);
    if (*endptr != '\0') {
        printf("parse error: '%s' is not a number\n", s);
        exit(1);
    }

    return ret;
}

void get_mtconf_options(unsigned int *nr_cpus, unsigned int **cpus) {
    unsigned int i;
    char *s, *e, *token;

    e = getenv(MT_CONF);
    if (!e) {
        printf("%s empty: setting default mt options: 0\n", MT_CONF);
        *nr_cpus = 1;
        *cpus = malloc(sizeof(unsigned int) * (*nr_cpus));
        if (!*cpus) {
            fprintf(stderr, "mt_get_options: malloc failed\n");
            exit(1);
        }

        *cpus[0] = 0;
        return;
    }

    s = malloc(strlen(e) + 1);
    if (!s) {
        fprintf(stderr, "mt_get_options: malloc failed\n");
        exit(1);
    }

    memcpy(s, e, strlen(e) + 1);
    *nr_cpus = 1;
    for (i = 0; i < strlen(s); i++) {
        if (s[i] == ',') {
            *nr_cpus = *nr_cpus + 1;
        }
    }

    i = 0;
    *cpus = malloc(sizeof(unsigned int) * (*nr_cpus));
    if (!(*cpus)) {
        fprintf(stderr, "mt_get_options: malloc failed\n");
        exit(1);
    }

    token = strtok(s, ",");
    do {
        (*cpus)[i++] = (unsigned int)parse_int(token);
    } while ((token = strtok(NULL, ",")));

    free(s);
    return;
}

void mt_conf_print(unsigned int ncpus, unsigned int *cpus) {
    unsigned int i;

    printf("MT_CONF=");
    for (i = 0; i < ncpus; i++) {
        if (i != 0)
            printf(",");
        printf("%u", cpus[i]);
    }
    printf("\n");
}
