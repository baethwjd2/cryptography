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
#include "mRSA.h"

/*
 * mod_add() - computes a + b mod m
 */
static uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    // a와 b이 m보다 작아지도록 나머지 연산
    a %= m;
    b %= m;

    if (a >= m - b) return a - (m - b);  // a+b가 m보다 크거나 같다면 오버플로가 발생하지 않도록 a-(m-b)
    else return a + b;    // 그렇지 않은 경우 a+b
}

/*
 * mod_mul() - computes a * b mod m
 */
static uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t res = 0;

    // double addition 알고리즘을 사용해 오버플로 없이 multiply 진행
    while (b > 0){ // b가 0이 될 때까지 반복
        // b의 가장 낮은 비트가 1일 때
        if (b & 1) res = mod_add(res, a, m);

        b >>= 1; // b를 오른쪽으로 한 칸 쉬프트
        a = mod_add(a, a, m);   // a는 2배가 됨
    }

    return res; // 결과값 반환
}

/*
 * mod_pow() - computes a^b mod m
 */
static uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t res = 1;

    // square multiplication 알고리즘을 사용해 오버플로 없이 지수 연산 수행
    while (b > 0){ // b가 0이 될 때까지 반복
        // b의 가장 낮은 비트가 1일 때
        if (b & 1) res = mod_mul(res, a, m);

        b >>= 1;    // b를 오른쪽으로 한 칸 쉬프트
        a = mod_mul(a, a, m);   // a는 제곱이 됨
    }

    return res; // 결과값 반환
}

/*
 * gcd() - Euclidean algorithm
 */
static uint64_t gcd(uint64_t a, uint64_t b)
{
    while(a!=0 && b!=0){    // a 혹은 b가 0이 되면 반복을 멈춤
        uint64_t tmp;

        // 유클리드 알고리즘에 따라 a는 b, b는 a mod b로 값을 바꿈
        tmp = a;
        a = b;
        b = tmp%b;
    }

    // 만약 b가 0이라면 a를, a가 0이라면 b를 리턴함
    if(b==0) return a;  
    else return b;
}

/*
 * mul_inv() - computes multiplicative inverse a^-1 mod m
 * It returns 0 if no inverse exist.
 */
static uint64_t mul_inv(uint64_t a, uint64_t m)
{
    uint64_t d0=a, d1=m, x0=1, x1=0;
    int x0_pos=1, x1_pos=0; // x0과 x1이 부호를 표시해줄 변수 선언 (1이면 양수, 0이면 음수)
   
    while (d1>1){
        uint64_t q, tmp;
        int tmp2;
        
        // 확장유클리드 알고리즘에 따라 q, d를 연산
        q = d0/d1;
        d0 = d0-q*d1;
        
        // x0과 x1의 부호에 따라 x를 연산
        if(x0_pos==x1_pos){ // x0과 x1의 부호가 서로 같을 때
        	x0 = (x0>q*x1) ? x0-q*x1 : q*x1-x0;
        	tmp2 = x0_pos; x0_pos = x1_pos; 
        	x1_pos = (x0>q*x1) ? tmp2 : !tmp2;
        }else{  // x0과 x1의 부호가 서로 다를 때
        	x0 = x0+q*x1; 
        	tmp2 = x0_pos; x0_pos = x1_pos;x1_pos = tmp2;
        } 

        // d, x 값을 업데이트함
        tmp=d0; d0=d1; d1=tmp;
        tmp=x0; x0=x1; x1=tmp;
    }

    // 만약 d가 1이라면, 역이 항상 존재하므로 역을 리턴해주고, 그렇지 않다면 0을 리턴함
    if(d1==1) return (x1_pos ? x1 : m-x1);  // d가 1일 때 x의 값이 음수라면 m에서 x를 빼서 리턴
    else return 0;  
}

/*
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3317044064679887385961981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
static const uint64_t a[BASELEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns 1 if n is prime, 0 otherwise.
 */
static int miller_rabin(uint64_t n)
{
    if (n == 2) return PRIME;  // n이 2일 때 PRIME 반환
    else if (n % 2 == 0 || n == 1) return COMPOSITE;   // n이 1이거나 짝수일 때 COMPOSITE 반환

    uint64_t k = 0, q = n - 1, tmp = 1;

    // k와 q를 선택함
    while (q % 2 == 0){
        k += 1;
        q /= 2;
    }

    // 베이스 집합에 대해 검증
    for (int i = 0; i < BASELEN; i++){
        int flag = 0;

        if (a[i] >= n - 1) break;    // n-1보다 a가 크면 검증할 필요 없으므로 검증 중단

        // a^q mod n 연산
        tmp = mod_pow(a[i], q, n);

        if (tmp == 1) continue;    // 연산 결과가 1이면 다음 a에 대해 검증 

        // a^(q+2^j) mod n 연산
        for (int j = 0; j < k; j++){
            if (tmp == n - 1){   // 연산 결과가 n-1이면 연산 중단
                flag = 1;
                break;
            }
            tmp = mod_mul(tmp, tmp, n); // 현재 값을 mod n에 대해 제곱
        }

        if (flag) continue;  // 방금 반복문을 탈출했었다면 다음 a에 대해 검증
        else return COMPOSITE;  // 0~k-1에 대해 만족하는 j가 없었다면 합성수 판정
    }

    return PRIME;   // 모든 검증을 통과했다면 소수 판정
}

/*
 * mRSA_generate_key() - generates mini RSA keys e, d and n
 *
 * Carmichael's totient function Lambda(n) is used.
 */

void mRSA_generate_key(uint64_t *e, uint64_t *d, uint64_t *n)
{
    uint64_t p, q, c;
    *n = 0;
    
    while(*n<MINIMUM_N){    // n이 64비트가 될 때까지 반복
        p = 1; q = 1;
        // 32비트이면서 소수일 때까지 반복하여 p와 q를 랜덤하게 생성함
        while(miller_rabin(p)!=PRIME || p<2147483648) arc4random_buf(&p, sizeof(uint32_t));
        while(miller_rabin(q)!=PRIME || q<2147483648) arc4random_buf(&q, sizeof(uint32_t));
       
        // p와 q를 곱해 n을 구함
        *n = p*q;  
    }
    
    // 카마이클 함수 값을 구함
    c = (p-1)*(q-1) / gcd(p-1, q-1);   

    // 카마이클 함수 값과의 최대공약수 값이 1이 될 때까지 e를 랜덤하게 생성함
    do{
        arc4random_buf(e, sizeof(uint64_t));
    }while (gcd(*e, c) != 1);

    // e의 역수를 계산해 d를 구함
    *d = mul_inv(*e, c);
}

/*
 * mRSA_cipher() - compute m^k mod n
 *
 * If data >= n then returns 1 (error), otherwise 0 (success).
 */
int mRSA_cipher(uint64_t *m, uint64_t k, uint64_t n)
{
    if(*m>=n) return 1; // m>=n일 때 오류로 처리하기 위해 1을 반환

    // 암호화 혹은 복호화를 위해 메세지에 모듈로 n에 대해 k제곱을 연산함
    *m = mod_pow(*m, k, n);

    return 0;   // 오류가 없었으므로 0을 반환
}
