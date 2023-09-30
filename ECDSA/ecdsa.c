/*
 * Copyright 2020-2022. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었다.
 */
#ifdef __linux__
#include <bsd/stdlib.h>
#elif __APPLE__
#include <stdlib.h>
#else
#include <stdlib.h>
#endif
#include "ecdsa.h"
#include "sha2.h"
#include <gmp.h>

mpz_t p, n;
ecdsa_p256_t G;

/*
 * hash 함수의 출력 길이를 바이트 단위로 반환하는 함수
 * Example: int h_len = select_hash(sha2_ndx);
 */
int select_hash(int sha2_ndx)
{
    switch (sha2_ndx)
    {
    case 0:
        return 224 / 8; // SHA224
    case 1:
        return 256 / 8; // SHA256
    case 2:
        return 384 / 8; // SHA384
    case 3:
        return 512 / 8; // SHA512
    case 4:
        return 224 / 8; // SHA512_224
    case 5:
        return 256 / 8; // SHA512_256
    }

    return -1;
}

/*
 * 해당하는 hash 함수를 실행하는 함수
 * Example: execute_hash(sha2_ndx, m, mLen, M);
 */
void execute_hash(int sha2_ndx, const void* m, size_t len, void* digest)
{
    switch (sha2_ndx)
    {
    case 0: // SHA224
        sha224(m, len, digest);
        break;
    case 1: // SHA256
        sha256(m, len, digest);
        break;
    case 2: // SHA384
        sha384(m, len, digest);
        break;
    case 3: // SHA512
        sha512(m, len, digest);
        break;
    case 4: // SHA512_224
        sha512_224(m, len, digest);
        break;
    case 5: // SHA512_256
        sha512_256(m, len, digest);
        break;
    }
}

/*
 * Initialize 256 bit ECDSA parameters
 * 시스템파라미터 p, n, G의 공간을 할당하고 값을 초기화한다.
 */
void ecdsa_p256_init(void)
{
    mpz_t Gx, Gy;

    mpz_inits(Gx, Gy, NULL);

    // 시스템 파라미터 p, n, Gx, Gy값 설정
    mpz_init_set_str(p, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
    mpz_init_set_str(n, "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551", 16);
    mpz_init_set_str(Gx, "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", 16);
    mpz_init_set_str(Gy, "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5", 16);

    // Gx, Gy는 G로 export
    mpz_export(G.x, NULL, 1, ECDSA_P256 / 8, 1, 0, Gx);
    mpz_export(G.y, NULL, 1, ECDSA_P256 / 8, 1, 0, Gy);

    mpz_clears(Gx, Gy, NULL);
}

/*
 * Clear 256 bit ECDSA parameters
 * 할당된 파라미터 공간을 반납한다.
 */
void ecdsa_p256_clear(void)
{
    // 할당된 시스템 파라미터 p, n 메모리 해제
    mpz_clears(p, n, NULL);
}

/*
 * P=Q일 때 elliptic curve addition을 수행하는 함수
 */
void eliptic_equal_add(mpz_t* Qx, mpz_t* Qy) {
    mpz_t common, temp_x, temp_y, tmp;
    mpz_init_set(temp_x, *Qx);
    mpz_init_set(temp_y, *Qy);
    mpz_inits(common, tmp, NULL);

    // common = (3 *Qx ** 2 -3) / 2 * Qy
    mpz_powm_ui(common, temp_x, 2, p);  // Qx^2
    mpz_mul_si(common, common, 3);  // 3*Qx^2
    mpz_sub_ui(common, common, 3);  // 3*Qx^2-3
    mpz_mul_si(tmp, temp_y, 2); // tmp = 2*Qy
    mpz_invert(tmp, tmp, p);    // tmp = 2*Qy^-1
    mpz_mul(common, common, tmp);   // common = (3*Qx^2-3) / tmp

    // Qx_next = ((3 * Qx ** 2 -3) / 2 * Qy) ** 2 - 2 * Qx
    mpz_powm_ui(*Qx, common, 2, p); // common^2
    mpz_submul_ui(*Qx, temp_x, 2);  // common-2Qx

    // Qy_next = ((3 * Qx ** 2 -3) / 2 * Qy) * (Qx - Qx_next) - Qy
    mpz_sub(tmp, temp_x, *Qx);   // Qx-Qx_next
    mpz_mul(*Qy, common, tmp);   // common * (Qx-Qx_next)
    mpz_sub(*Qy, *Qy, temp_y);  // common * (Qx-Qx_next) - Qy 

    // GF(p)에 정의된 곡선이므로 mod p
    mpz_mod(*Qy, *Qy, p);
    mpz_mod(*Qx, *Qx, p);

    // 사용했던 mpz 변수들 메모리 해제
    mpz_clears(common, temp_x, temp_y, tmp, NULL);
}

/*
 *  P!=Q일 때 elliptic curve addition을 수행하는 함수
 */
void eliptic_diff_add(mpz_t* Px, mpz_t* Py, mpz_t* Qx, mpz_t* Qy) {
    mpz_t common, temp_px, temp_py, temp_qx, temp_qy, tmp;
    mpz_init_set(temp_qx, *Qx);
    mpz_init_set(temp_qy, *Qy);
    mpz_init_set(temp_px, *Px);
    mpz_init_set(temp_py, *Py);
    mpz_init_set(common, temp_qy);
    mpz_init_set(tmp, temp_qx);

    // common = {(Qy-Py)/(Qx-Px)}
    mpz_sub(common, common, temp_py);    // Qy - Py
    mpz_sub(tmp, tmp, temp_px);   // tmp = Qx - Px

    mpz_invert(tmp, tmp, p);    // tmp^-1
    mpz_mul(common, common, tmp);   // common = (Qy - Py) * tmp^-1

    // Qx_next = {(Qy-Py)/(Qx-Px)}^2 - Qx - Px
    mpz_powm_ui(*Qx, common, 2, p); // common^2
    mpz_sub(*Qx, *Qx, temp_qx); // common^2 - Qx
    mpz_sub(*Qx, *Qx, temp_px); // common^2 - Qx - Px

    // Qy_next = {(Qy-Py)/(Qx-Px)} * (Px - Qx_next) - Py
    mpz_sub(tmp, *Px, *Qx);  // (Px - Qx_next)
    mpz_mul(*Qy, common, tmp);  // {(Qy-Py)/(Qx-Px)} * (Px - Qx_next)
    mpz_sub(*Qy, *Qy, temp_py); // {(Qy-Py)/(Qx-Px)} * (Px - Qx_next) - Py

    // GF(p)에 정의된 곡선이므로 mod p
    mpz_mod(*Qy, *Qy, p);
    mpz_mod(*Qx, *Qx, p);

    // 사용했던 mpz 변수들 메모리 해제
    mpz_clears(common, temp_px, temp_py, temp_qx, temp_qy, tmp, NULL);
}

/*
 * Q = dT 에서 d, T를 입력받아 결과를 Qx, Qy에 넘겨주는 함수
 */
void eliptic_mul(mpz_t d, mpz_t* Qx, mpz_t* Qy, const ecdsa_p256_t T) {
    mpz_t Px, Py;

    // 변수들 0으로 초기화
    mpz_inits(Px, Py, NULL);

    // 항등원 O가 있을 때 지정할 flag
    int first_addition_flag = 1;

    // T의 x, y 불러오기    
    mpz_import(Px, ECDSA_P256 / 8, 1, 1, 1, 0, T.x);
    mpz_import(Py, ECDSA_P256 / 8, 1, 1, 1, 0, T.y);

    // Q = dT 계산
    for (int i = 0; i < ECDSA_P256; i++) {
        if (mpz_tstbit(d, i)){
            if (first_addition_flag) {   // O + Q = Q 수행
                mpz_set(*Qx, Px);
                mpz_set(*Qy, Py);
                first_addition_flag = 0;
            }
            else {    // Q에 값 더하기
                eliptic_diff_add(&Px, &Py, Qx, Qy);
                mpz_mod(*Qx, *Qx, p);
                mpz_mod(*Qy, *Qy, p);
            }
        }

        // 만약 Q가 항등원이라면 flag 지정
        if (mpz_cmp_si(*Qx, 0) == 0 && mpz_cmp_si(*Qy, 0) == 0) first_addition_flag = 1;

        // 더할 값 P 두배하기
        eliptic_equal_add(&Px, &Py);
        mpz_mod(Px, Px, p);
        mpz_mod(Py, Py, p);
    }

    // 사용했던 mpz 변수들 메모리 해제
    mpz_clears(Px, Py, NULL);
}


/*
 * ecdsa_p256_key() - generates Q = dG
 * 사용자의 개인키와 공개키를 무작위로 생성한다.
 */
void ecdsa_p256_key(void* d, ecdsa_p256_t* Q)
{
    // d: 개인키, n보다 작은 소수
    // Q: 공개키, d와 G의 곱
    mpz_t d_, Qx, Qy;
    // 변수들 0으로 초기화
    mpz_inits(d_, Qx, Qy, NULL);
    gmp_randstate_t state;

    // state에 random값 초기화
    gmp_randinit_default(state);
    gmp_randseed_ui(state, arc4random());

    // 개인키 d를 랜덤하게 생성
    do {
        mpz_urandomm(d_, state, n); // d 는 [1, n-1]사이
    } while (mpz_probab_prime_p(d_, 50));   // d가 소수인지 검증 

    // mpz_t인 d_의 값을 d에 넘겨주기
    mpz_export(d, NULL, 1, ECDSA_P256 / 8, 1, 0, d_);
    eliptic_mul(d_, &Qx, &Qy, G);   // Q = d * G

    // d와 Q를 넘겨줌 
    mpz_export(Q->x, NULL, 1, ECDSA_P256 / 8, 1, 0, Qx);
    mpz_export(Q->y, NULL, 1, ECDSA_P256 / 8, 1, 0, Qy);

    // 사용했던 mpz 변수들과 state 메모리 해제
    mpz_clears(d_, Qx, Qy, NULL);
    gmp_randclear(state);
}

/*
 * ecdsa_p256_sign(msg, len, d, r, s) - ECDSA Signature Generation
 * 길이가 len 바이트인 메시지 m을 개인키 d로 서명한 결과를 r, s에 저장한다.
 * sha2_ndx는 사용할 SHA-2 해시함수 색인 값으로 SHA224, SHA256, SHA384, SHA512,
 * SHA512_224, SHA512_256 중에서 선택한다. r과 s의 길이는 256비트이어야 한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int ecdsa_p256_sign(const void* msg, size_t len, const void* d, void* _r, void* _s, int sha2_ndx)
{
    unsigned int h_len = select_hash(sha2_ndx);
    unsigned char tmp_e[h_len];
    mpz_t s, r, d_tmp, x1, y1, k, e, tmp;
    gmp_randstate_t state;

    // 해시함수가 허용하는 입력 (2^64 bits = 2^61 bytes) 보다 크면 오류 반환
    if (sha2_ndx < 4 && len >= 0x2000000000000000) return ECDSA_MSG_TOO_LONG;

    // 변수들 0으로 초기화
    mpz_inits(s, r, d_tmp, x1, y1, k, e, tmp, NULL);

    // randstate 생성
    gmp_randinit_default(state);
    gmp_randseed_ui(state, arc4random());

    // d 가져오기
    mpz_import(d_tmp, ECDSA_P256 / 8, 1, 1, 1, 0, d);

    // e = H(m)
    execute_hash(sha2_ndx, msg, len, tmp_e);

    // e가 n의 길이(256 bits)보다 길면 자름
    if (h_len > ECDSA_P256 / 8) h_len = ECDSA_P256 / 8;
    mpz_import(e, h_len, 1, 1, 1, 0, tmp_e);

    while (mpz_cmp_si(s, 0) == 0){  // s가 0, 즉 n의 배수가 아닐 때 까지 반복
        while (mpz_cmp_si(r, 0) == 0){ // r이 0, 즉 n의 배수가 아닐 때까지 반복
            // 비밀값 k 랜덤으로 선택
            mpz_urandomm(k, state, n);  //  (0 < 𝑘 < 𝑛)

            // (x1, y1) = kG
            eliptic_mul(k, &x1, &y1, G);

            // r 연산, r = x1 mod n
            mpz_mod(r, x1, n);
        }

        // s 연산, s = k^-1 * (e + r * d) mod n
        mpz_invert(k, k, n);    // k^-1
        mpz_mul(tmp, r, d_tmp); // r * d
        mpz_add(e, e, tmp); // e + r * d
        mpz_mul(k, k, e);  // k^-1 * (e + r * d)
        mpz_mod(s, k, n);
    }

    // _r, _s에 서명 결과 넘겨주기
    mpz_export(_s, NULL, 1, ECDSA_P256 / 8, 1, 0, s);
    mpz_export(_r, NULL, 1, ECDSA_P256 / 8, 1, 0, r);

    // 선언했던 mpz 메모리 해제
    mpz_clears(s, r, d_tmp, x1, y1, k, e, tmp, NULL);
    gmp_randclear(state);

    // 서명에 성공하였으므로 0 반환
    return 0;
}

/*
 * ecdsa_p256_verify(msg, len, Q, r, s) - ECDSA signature veryfication
 * It returns 0 if valid, nonzero otherwise.
 * 길이가 len 바이트인 메시지 m에 대한 서명이 (r,s)가 맞는지 공개키 Q로 검증한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int ecdsa_p256_verify(const void* msg, size_t len, const ecdsa_p256_t* _Q, const void* _r, const void* _s, int sha2_ndx)
{
    unsigned int h_len = select_hash(sha2_ndx);
    unsigned char tmp_e[h_len];
    mpz_t r, s, x1, y1, u1, u2, u1_tmp_x, u1_tmp_y, u2_tmp_x, u2_tmp_y, e;

    // mpz 변수들을 0으로 초기화
    mpz_inits(r, s, x1, y1, u1, u2, u1_tmp_x, u1_tmp_y, u2_tmp_x, u2_tmp_y, e, NULL);

    // r,s 가져오기
    mpz_import(r, ECDSA_P256 / 8, 1, 1, 1, 0, _r);
    mpz_import(s, ECDSA_P256 / 8, 1, 1, 1, 0, _s);

    // 해시함수가 허용하는 입력 (2^64-1bit = 2^61 byte) 보다 크면 오류 반환
    if (sha2_ndx < 4 && len >= 0x2000000000000000) return ECDSA_MSG_TOO_LONG;

    // r, s가 1~n-1 사이에 있지 않으면 잘못된 서명이므로 에러 반환
    if (mpz_cmp_si(r, 1) < 0 || mpz_cmp(r, n) > 0)  return ECDSA_SIG_INVALID;
    if (mpz_cmp_si(s, 1) < 0 || mpz_cmp(s, n) > 0)  return ECDSA_SIG_INVALID;

    // e = H(m)
    execute_hash(sha2_ndx, msg, len, tmp_e);

    // e가 n의 길이(256 bits)보다 길면 자름
    if (h_len > ECDSA_P256 / 8) h_len = ECDSA_P256 / 8;
    mpz_import(e, h_len, 1, 1, 1, 0, tmp_e);

    // u1과 u2 연산, u1 = e*s^-1 mod n, u2 = r*s^-1 mod n
    mpz_invert(s, s, n); // s^-1
    mpz_mul(u1, e, s); // e * s^-1 
    mpz_mul(u2, r, s); // r * s^-1
    mpz_mod(u1, u1, n); // u1 = e * s^-1 mod n
    mpz_mod(u2, u2, n); // u2 = r * s^-1 mod n

    // x1, y1연산, (x1, y1) = u1G + u2Q
    eliptic_mul(u1, &u1_tmp_x, &u1_tmp_y, G);   // u1*G
    eliptic_mul(u2, &u2_tmp_x, &u2_tmp_y, *_Q); // u2*Q

    // u1*G + u2*Q를 계산하는데 u1*G와 u2*Q가 같다면 (x1,y1)은 접선과 만나는 다른 교점
    if (mpz_cmp(u1_tmp_x, u2_tmp_x) == 0 && mpz_cmp(u1_tmp_y, u2_tmp_y) == 0) {
        // y1 == 0일때(x1, y1) == O이므로 오류 반환
        if (mpz_cmp_si(u1_tmp_y, 0) == 0) return ECDSA_SIG_INVALID;
        eliptic_equal_add(&u2_tmp_x, &u2_tmp_y);
    }
    else //다르다면 (x1,y1)은 두 점을 이어 나온 직선과 만나는 새로운 교점
    {
        // x와 y의 값이 모두 같은 경우를 위에서 제외했으므로 아래 if는 x가 같고 y가 다르다는 뜻.
        // x가 같고 y의 부호가 반대인 경우 (x1, y1) = O. 따라서 오류 반환
        if (mpz_cmp(u1_tmp_x, u2_tmp_x) == 0) return ECDSA_SIG_INVALID;
        eliptic_diff_add(&u1_tmp_x, &u1_tmp_y, &u2_tmp_x, &u2_tmp_y);
    }

    // x1, y1 값 설정
    mpz_set(x1, u2_tmp_x);
    mpz_set(y1, u2_tmp_y);

    // r과 x1이 mod n에서 합동이면 서명이 검증되었으므로 0을 반환
    if (mpz_cmp(r, x1) == 0) {
        mpz_clears(r, s, x1, y1, u1, u2, u1_tmp_x, u1_tmp_y, u2_tmp_x, u2_tmp_y, e, NULL);
        return 0;
    }
    else {
        // 선언했던 mpz 변수들 메모리 해제
        mpz_clears(r, s, x1, y1, u1, u2, u1_tmp_x, u1_tmp_y, u2_tmp_x, u2_tmp_y, e, NULL);
        return ECDSA_SIG_MISMATCH;
    }
}
