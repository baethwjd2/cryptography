/*
 * Copyright 2020-2022. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었다.
 */
#ifndef _ECDSA_H_
#define _ECDSA_H_

/*
 * 타원곡선 P-256의 그룹 소수와 차수의 비트 크기로 값을 임의로 변경해서는 안된다.
 */
#define ECDSA_P256 256

/*
 * SHA-2 계열의 해시함수를 구분하기 위한 색인 값이다.
 * SHA512_224와 SHA512_256는 각각 SHA512/224와 SHA512/256를 의미한다.
 */
#define SHA224      0
#define SHA256      1
#define SHA384      2
#define SHA512      3
#define SHA512_224  4
#define SHA512_256  5

/*
 * 오류 코드 목록이다. 오류가 없으면 0을 사용한다.
 */
#define ECDSA_MSG_TOO_LONG  1
#define ECDSA_SIG_INVALID   2
#define ECDSA_SIG_MISMATCH  3

/*
 * 타원곡선 P-256 상의 점을 나타내기 위한 구조체이다.
 */
typedef struct {
    unsigned char x[ECDSA_P256/8];
    unsigned char y[ECDSA_P256/8];
} ecdsa_p256_t;

void ecdsa_p256_init(void);
void ecdsa_p256_clear(void);
void ecdsa_p256_key(void *d, ecdsa_p256_t *Q);
int ecdsa_p256_sign(const void *msg, size_t len, const void *d, void *r, void *s, int sha2_ndx);
int ecdsa_p256_verify(const void *msg, size_t len, const ecdsa_p256_t *Q, const void *r, const void *s, int sha2_ndx);

#endif
