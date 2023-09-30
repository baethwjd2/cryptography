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
#include <string.h>
#include <gmp.h>
#include "pkcs.h"
#include "sha2.h"

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
mgf 함수
mgfSeed : mgf를 적용해줄 seed의 포인터
seedLen : seed의 길이
mask : mgf 적용 결과값을 저장해줄 변수
maskLen : 원하는 mgf 적용 결과값의 길이
*/
static unsigned char* mgf(int sha2_ndx, void* mgfSeed, size_t seedLen, void* mask, size_t maskLen)
{
    uint32_t i, count, c;
    unsigned char* mgfIn, * m;
    int hLen = select_hash(sha2_ndx);

    if (maskLen > 0x0100000000 * hLen) // maskLen이 2^32 * hLen보다 긴지 체크
        return NULL;

    // octet string mask 생성
    mgfIn = (unsigned char*)malloc(seedLen + 4);
    memset(mgfIn, 0, seedLen + 4);
    memcpy(mgfIn, mgfSeed, seedLen);
    count = maskLen / hLen + (maskLen % hLen ? 1 : 0);
    m = (unsigned char*)malloc(count * hLen);
    memset(m, 0, count * hLen);

    // i를 옥텟 문자열 c의 길이로 변환
    // mgfSeed의 해쉬와 C를 octet string에 병합
    for (i = 0; i < count; i++)
    {
        c = i;
        mgfIn[seedLen + 3] = c & 0x000000ff;
        c >>= 8;
        mgfIn[seedLen + 2] = c & 0x000000ff;
        c >>= 8;
        mgfIn[seedLen + 1] = c & 0x000000ff;
        c >>= 8;
        execute_hash(sha2_ndx, mgfIn, seedLen + 4, m + i * hLen);
    }

    // mask copy 후 memory 해제
    memcpy(mask, m, maskLen);
    free(mgfIn);
    free(m);
    return mask;
}

/*
 * rsa_generate_key() - generates RSA keys e, d and n in octet strings.
 * If mode = 0, then e = 65537 is used. Otherwise e will be randomly selected.
 * Carmichael's totient function Lambda(n) is used.
 */
void rsa_generate_key(void* _e, void* _d, void* _n, int mode)
{
    mpz_t p, q, lambda, e, d, n, gcd;
    gmp_randstate_t state;

    /*
     * Initialize mpz variables
     */
    mpz_inits(p, q, lambda, e, d, n, gcd, NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, arc4random());
    /*
     * Generate prime p and q such that 2^(RSAKEYSIZE-1) <= p*q < 2^RSAKEYSIZE
     */
    do
    {
        do
        {
            mpz_urandomb(p, state, RSAKEYSIZE / 2);
            mpz_setbit(p, 0);
            mpz_setbit(p, RSAKEYSIZE / 2 - 1);
        } while (mpz_probab_prime_p(p, 50) == 0);
        do
        {
            mpz_urandomb(q, state, RSAKEYSIZE / 2);
            mpz_setbit(q, 0);
            mpz_setbit(q, RSAKEYSIZE / 2 - 1);
        } while (mpz_probab_prime_p(q, 50) == 0);
        mpz_mul(n, p, q);
    } while (!mpz_tstbit(n, RSAKEYSIZE - 1));
    /*
     * Generate e and d using Lambda(n)
     */
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_lcm(lambda, p, q);
    if (mode == 0)
        mpz_set_ui(e, 65537);
    else
        do
        {
            mpz_urandomb(e, state, RSAKEYSIZE);
            mpz_gcd(gcd, e, lambda);
        } while (mpz_cmp(e, lambda) >= 0 || mpz_cmp_ui(gcd, 1) != 0);
        mpz_invert(d, e, lambda);
        /*
         * Convert mpz_t values into octet strings
         */
        mpz_export(_e, NULL, 1, RSAKEYSIZE / 8, 1, 0, e);
        mpz_export(_d, NULL, 1, RSAKEYSIZE / 8, 1, 0, d);
        mpz_export(_n, NULL, 1, RSAKEYSIZE / 8, 1, 0, n);
        /*
         * Free the space occupied by mpz variables
         */
        mpz_clears(p, q, lambda, e, d, n, gcd, NULL);
}

/*
 * rsa_cipher() - compute m^k mod n
 * If m >= n then returns PKCS_MSG_OUT_OF_RANGE, otherwise returns 0 for success.
 */
static int rsa_cipher(void* _m, const void* _k, const void* _n)
{
    mpz_t m, k, n;

    /*
     * Initialize mpz variables
     */
    mpz_inits(m, k, n, NULL);
    /*
     * Convert big-endian octets into mpz_t values
     */
    mpz_import(m, RSAKEYSIZE / 8, 1, 1, 1, 0, _m);
    mpz_import(k, RSAKEYSIZE / 8, 1, 1, 1, 0, _k);
    mpz_import(n, RSAKEYSIZE / 8, 1, 1, 1, 0, _n);
    /*
     * Compute m^k mod n
     */
    if (mpz_cmp(m, n) >= 0)
    {
        mpz_clears(m, k, n, NULL);
        return PKCS_MSG_OUT_OF_RANGE;
    }
    mpz_powm(m, m, k, n);
    /*
     * Convert mpz_t m into the octet string _m
     */
    mpz_export(_m, NULL, 1, RSAKEYSIZE / 8, 1, 0, m);
    /*
     * Free the space occupied by mpz variables
     */
    mpz_clears(m, k, n, NULL);
    return 0;
}

/*
 * rsaes_oaep_encrypt() - RSA encrytion with the EME-OAEP encoding method
 * 길이가 len 바이트인 메시지 m을 공개키 (e,n)으로 암호화한 결과를 c에 저장한다.
 * label은 데이터를 식별하기 위한 라벨 문자열로 NULL을 입력하여 생략할 수 있다.
 * sha2_ndx는 사용할 SHA-2 해시함수 색인 값으로 SHA224, SHA256, SHA384, SHA512,
 * SHA512_224, SHA512_256 중에서 선택한다. c의 크기는 RSAKEYSIZE와 같아야 한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int rsaes_oaep_encrypt(const void* m, size_t mLen, const void* label, const void* e, const void* n, void* c, int sha2_ndx)
{
    // 변수 생성 
    int mHash_byte = select_hash(sha2_ndx);
    int hLen = mHash_byte;
    int SHASIZE = mHash_byte * 8;
    size_t RSA_byte = RSAKEYSIZE / 8;
    size_t DB_byte = RSAKEYSIZE / 8 - mHash_byte - 1;
    if (strlen(label) > hLen) // label이 hash할 수 있는 범위보다 길 경우 오류 반환
        return PKCS_LABEL_TOO_LONG;
    if (mLen > RSA_byte - 2 * hLen - 2) // message의 길이가 너무 길 경우 오류 반환
        return PKCS_MSG_TOO_LONG;
    unsigned char Hlabel[SHASIZE];
    unsigned char DB[DB_byte];
    unsigned char EM[RSAKEYSIZE], maskedDB[DB_byte];
    unsigned char seed[mHash_byte];
    unsigned char maskedSeed[mHash_byte];
    size_t PS_byte = DB_byte - mHash_byte - mLen - 1;

    // 사용할 sha에 맞게 arc4random함수를 통해 seed 생성
    *seed = arc4random_uniform(SHASIZE);

    // DB 생성 과정 
    execute_hash(sha2_ndx, label, strlen(label), Hlabel); // 사용할 해시 함수에 맞게 label을 해시한다
    memcpy(DB, Hlabel, mHash_byte); // 해시한 label을 DB에 추가한다 
    memset(DB + mHash_byte, 0x00, PS_byte); // Padding String 길이에 맞게 DB에 0x00을 추가한다 
    DB[mHash_byte + PS_byte] = 0x01; // Padding String과 Message를 구분하기 위해 0x01을 추가한다 
    memcpy(DB + mHash_byte + PS_byte + 1, m, mLen); // 마지막으로 DB에 Message를 추가한다

    // maskedDB 생성
    mgf(sha2_ndx, seed, mHash_byte, maskedDB, DB_byte);
    for (int i = 0; i < DB_byte; i++)
    {
        maskedDB[i] ^= DB[i];
    }

    // maskedSeed 생성
    mgf(sha2_ndx, maskedDB, DB_byte, maskedSeed, mHash_byte);
    for (int j = 0; j < mHash_byte; j++)
    {
        maskedSeed[j] ^= seed[j];
    }

    // EM 생성
    EM[0] = 0x00; // 초기에 0x00을 추가 
    memcpy(EM + 1, maskedSeed, mHash_byte); // EM에 Masked Seed를 추가 
    memcpy(EM + 1 + mHash_byte, maskedDB, DB_byte); // EM에 Masked DB를 추가 

    // Encoded Message 검사
    if (rsa_cipher(EM, e, n))
        return PKCS_MSG_OUT_OF_RANGE;

    // c에 EM 복사
    memcpy(c, EM, RSA_byte);

    return 0;
}

/*
 * rsaes_oaep_decrypt() - RSA decrytion with the EME-OAEP encoding method
 * 암호문 c를 개인키 (d,n)을 사용하여 원본 메시지 m과 길이 len을 회복한다.
 * label과 sha2_ndx는 암호화할 때 사용한 것과 일치해야 한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int rsaes_oaep_decrypt(void* m, size_t* mLen, const void* label, const void* d, const void* n, const void* c, int sha2_ndx)
{
    // 기본적으로 cipher encrypt를 반대로 하는 것이다.
    // 앞 과정이 seed 생성, Data Block 생성 (Padding 삽입 및 Hash 계산), MGF hash 통과 및 XOR, 00 붙여서 EM 생성, 공개키 암호화
    // 거꾸로 뒤집으면? 개인키 복호화, EM 앞 00 제거 및 masked seed와 masked message로 나누기,
    // MGF hash 통과 및 XOR로 원본 seed 계산 후 data block 복원,
    // 마지막으로 Data block에서 padding 제거 및 메시지 복원

    // 저렇게만 보면 lable hash mismatch 부분에 대한 내용이 없는데, data block 복원 후에 확인하면 된다.
    // 같은 label로 hash가 일치하는지 확인해야 한다는 것.

    int mHash_byte = select_hash(sha2_ndx); // label hash 분리 용도
    int SHASIZE = mHash_byte * 8;
    size_t RSA_byte = RSAKEYSIZE / 8;
    size_t DB_byte = DB_byte = RSAKEYSIZE / 8 - mHash_byte - 1;

    if (strlen(label) > SHASIZE / 8)
        return PKCS_LABEL_TOO_LONG;
    // msg len은 여기서 직접 구해야하는 것이라 확인할 수 있는 부분이 아니다.

    // 사용하는 SHA가 정해지니 seed의 길이와 label hash의 길이, Data Block의 길이까지는 문제없이 구해진다.
    // 다만 message의 길이 mLen을 직접 알아내야 하는 것이니, 마찬가지로 padding string의 길이 PS_byte는 구할 수 없다.
    unsigned char Hlabel[SHASIZE];
    unsigned char DB[DB_byte];
    unsigned char PS_message[DB_byte - mHash_byte];
    // PS_message는 padding string + 0x01 + message
    // 이것을 추가한 이유는? 메모리 겹침(overlap) 방지
    // memmove대신 사용
    unsigned char EM[RSAKEYSIZE], maskedDB[DB_byte];
    unsigned char seed[mHash_byte];
    unsigned char maskedSeed[mHash_byte];
    size_t PS_byte;

    // EM에 c 복사. 이때 EM은 아직 encrypt된 상태.
    memcpy(EM, c, RSA_byte);

    // decoded Message 검사
    if (rsa_cipher(EM, d, n))
        return PKCS_MSG_OUT_OF_RANGE;
    // 이 결과로 EM 안엔 온전히 decrypt 된 EM이 남게 된다.

    //일단 Encrypted message 첫 byte가 0x00인지 확인
    if (EM[0] != 0x00)
        return PKCS_INITIAL_NONZERO;

    // EM을 maskedSeed와 masekdDB로 분리
    memcpy(maskedSeed, EM + 1, mHash_byte);         // 앞 8bit짜리 0x00을 제외하고 남은 mHash_byte만큼
    memcpy(maskedDB, EM + 1 + mHash_byte, DB_byte); //마찬가지로 0x00과 masked seed를 제외한 나머지를 memcpy

    // 원본 seed 및 Data Block 복원
    // 원본 seed는 masked data block을 MGF 돌린 다음 XOR 하면 나온다.
    mgf(sha2_ndx, maskedDB, DB_byte, seed, mHash_byte); // maskedDB를 mgf 돌려서 seed에 넣어두고
    for (int i = 0; i < mHash_byte; i++)
    {
        seed[i] ^= maskedSeed[i]; // maskedSeed와 xor 하면 원본 seed가 나오게 된다.
    }

    // 원본 Data Block은 seed를 MGF 돌린 다음 XOR 하면 나온다.
    mgf(sha2_ndx, seed, mHash_byte, DB, DB_byte); // 얻어낸 원본 seed를 mgf 돌려 DB에 넣어두고
    for (int j = 0; j < DB_byte; j++)
    {
        DB[j] ^= maskedDB[j]; // maskedDB와 xor 해 원본 DB를 얻어낸다.
    }

    // DB를 Hlabel, Padding, 식별자 0x01, m으로 분리한다.
    execute_hash(sha2_ndx, label, strlen(label), Hlabel); // decode 자체에 필요한건 아니지만 hash mismatch 확인을 위해 만든 것.
    // hash를 돌려 그 값이 동일한지 본다.

    for (int i = 0; i < mHash_byte; i++)
    {
        if (DB[i] != Hlabel[i]) // char를 비교해 다르다면
        {
            return PKCS_HASH_MISMATCH;
        }
    }

    // hash에 문제가 없는 것을 확인했다면
    memcpy(PS_message, DB + mHash_byte, DB_byte - mHash_byte); // PS_message에 DB앞부분 Hlabel을 제외한 값을 넣고

    // 원래대로라면 PS_byte 만큼 0x00을 넣는 과정이나 여기선 확실하게 확인해야한다.
    // padding string을 제거할 때 식별자가 온전하지 않으면 에러.
    // padding string 0x00을 반복해 제거(정확히는 그 위치를 저장)하다가 0x00이 아닌 첫 위치를 확인하는 것.
    for (PS_byte = 0; PS_byte < DB_byte - mHash_byte; ++PS_byte)
    {
        if (PS_message[PS_byte] != 0x00)
        {
            if (PS_message[PS_byte] == 0x01) // 0x01을 만나면 나간다
            {
                break;
            }
            else // padding string이 끝났는데 0x01이 아닌 경우.
            {
                return PKCS_INVALID_PS;
            }
        }
    }

    // Data Block에서 Hash, padding과 식별자를 제외한 나머지가 message
    // mLen이 구해지니, 나머지를 m에 복사
    *mLen = DB_byte - mHash_byte - PS_byte - 1;
    memcpy(m, PS_message + PS_byte + 1, (size_t)*mLen);

    return 0;
}

/*
 * rsassa_pss_sign - RSA Signature Scheme with Appendix
 * 길이가 len 바이트인 메시지 m을 개인키 (d,n)으로 서명한 결과를 s에 저장한다.
 * s의 크기는 RSAKEYSIZE와 같아야 한다. 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int rsassa_pss_sign(const void* m, size_t mLen, const void* d, const void* n, void* s, int sha2_ndx)
{
    // 메시지 길이가 hash 함수의 input 제한(2^61-1)보다 크면 오류 반환
    if (mLen >= 0x2000000000000000)
        return PKCS_MSG_TOO_LONG;

    unsigned int h_len = select_hash(sha2_ndx), M_len = 2 * h_len + 8, DB_len = RSAKEYSIZE / 8 - h_len - 1, ps_len = DB_len - 1 - h_len;
    unsigned char mHash[h_len], M[M_len], DB[DB_len], mask[DB_len], EM[RSAKEYSIZE / 8], salt[h_len];

    // hash(m)을 연산해 mHash 생성
    execute_hash(sha2_ndx, m, mLen, mHash);

    // M'에 0s, mHash, salt 추가하기
    memset(M, 0, 8);             // 0s 추가
    memcpy(M + 8, mHash, h_len); // mHash 추가
    arc4random_buf(salt, h_len);
    memcpy(M + 8 + h_len, salt, h_len); // salt 추가

    // H(hash(M')) 생성
    execute_hash(sha2_ndx, M, M_len, mHash);

    // DB 생성
    memset(DB, 0, ps_len);                // PS(000..) 채우기
    memset(DB + ps_len, 0x01, 1);         // 0x01 채우기
    memcpy(DB + ps_len + 1, salt, h_len); // salt 채우기

    // maskedDB 생성
    mgf(sha2_ndx, mHash, h_len, mask, DB_len); // mgf(H) 연산
    for (int i = 0; i < DB_len; i++)
        DB[i] = DB[i] ^ mask[i]; // DB XOR MGF(H)

    // EM 생성
    memcpy(EM, DB, DB_len);                   // maskedDB 복사
    memcpy(EM + DB_len, mHash, h_len);        // hash(hash(M')) 복사
    memset(EM + RSAKEYSIZE / 8 - 1, 0xbc, 1); // TF 채우기
    EM[0] &= 0x7f;                            // EM의 MSB를 0으로 설정

    // RSA 암호화
    rsa_cipher(EM, d, n);
    memcpy(s, EM, RSAKEYSIZE / 8);

    return 0;
}

/*
 * rsassa_pss_verify - RSA Signature Scheme with Appendix
 * 길이가 len 바이트인 메시지 m에 대한 서명이 s가 맞는지 공개키 (e,n)으로 검증한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int rsassa_pss_verify(const void* m, size_t mLen, const void* e, const void* n, const void* s, int sha2_ndx)
{
    // 메시지 길이가 hash 함수의 input 제한(2^61-1)보다 크면 오류 반환
    if (mLen >= 0x2000000000000000)
        return PKCS_MSG_TOO_LONG;

    int h_len = select_hash(sha2_ndx), M_len = 2 * h_len + 8, DB_len = RSAKEYSIZE / 8 - h_len - 1, ps_len = DB_len - 1 - h_len;
    char mHash[h_len], H[h_len], M[M_len], DB[DB_len], mask[DB_len], EM[RSAKEYSIZE / 8], salt[h_len];

    // 서명으로부터 EM 추출
    memcpy(EM, s, RSAKEYSIZE / 8);
    rsa_cipher(EM, e, n);

    // TF가 0xbc와 일치하지 않으면 오류 반환
    if (EM[RSAKEYSIZE / 8 - 1] != (char)0xbc)
        return PKCS_INVALID_LAST;

    // EM의 첫 비트가 0이 아니면 오류 반환
    if ((EM[0] >> 7) & 1)
        return PKCS_INVALID_INIT;

    // maskedDB 추출
    memcpy(DB, EM, DB_len);

    // H(=hash(M')) 추출
    memcpy(H, EM + DB_len, h_len);

    // DB 추출
    mgf(sha2_ndx, H, h_len, mask, DB_len); // dbMask 생성
    for (int i = 0; i < DB_len; i++)
        DB[i] = DB[i] ^ mask[i]; // maskedDB XOR dbMask
    DB[0] &= 0x7f;               // DB의 MSB 0으로 설정

    // DB의 앞이 0x00...||0x01과 일치하지 않으면 오류 반환
    for (int i = 0; i < ps_len; i++)
        if (DB[i] != (char)0x00)
            return PKCS_INVALID_PD2;
    if (DB[ps_len] != (char)0x01)
        return PKCS_INVALID_PD2;

    // salt 추출
    memcpy(salt, DB + ps_len + 1, h_len);

    // part M'
    execute_hash(sha2_ndx, m, mLen, mHash);

    // M' 생성
    memset(M, 0, 8);                    // 0s 추가
    memcpy(M + 8, mHash, h_len);        // mHash 추가
    memcpy(M + 8 + h_len, salt, h_len); // salt 추가

    // H'(=hash(M')) 생성
    execute_hash(sha2_ndx, M, M_len, mHash);

    // H = H'이면 0을(검증 성공) 반환하고, 다르다면 오류 반환
    if (memcmp(mHash, H, h_len) == 0)
        return 0;
    else
        return PKCS_HASH_MISMATCH;
}