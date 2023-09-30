/*
 * Copyright 2020-2022. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었다.
 */
#ifndef _mRSA_H_
#define _mRSA_H_

#include <stdint.h>

#define BASELEN 12
#define PRIME 1
#define COMPOSITE 0
#define MINIMUM_N 0x8000000000000000

void mRSA_generate_key(uint64_t *e, uint64_t *d, uint64_t *n);
int mRSA_cipher(uint64_t *m, uint64_t k, uint64_t n);

#endif
