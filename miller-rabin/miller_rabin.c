/*
 * Copyright 2020-2022. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었다.
 */
#include "miller_rabin.h"

/*
 * mod_add() - computes a+b mod m
 * a와 b가 m보다 작다는 가정하에서 a+b >= m이면 결과에서 m을 빼줘야 하므로
 * 오버플로가 발생하지 않도록 a-(m-b)를 계산하고, 그렇지 않으면 그냥 a+b를 계산하면 된다.
 * a+b >= m을 검사하는 과정에서 오버플로가 발생할 수 있으므로 a >= m-b를 검사한다.
 */
uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    // a와 b이 m보다 작아지도록 나머지 연산
    a %= m;
    b %= m;

    if(a>=m-b) return a-(m-b);  // a+b가 m보다 크거나 같다면 오버플로가 발생하지 않도록 a-(m-b)
    else return a+b;    // 그렇지 않은 경우 a+b
    
}

/*
 * mod_sub() - computes a-b mod m
 * 만일 a < b이면 결과가 음수가 되므로 m을 더해서 양수로 만든다.
 */
uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t m)
{
    // a와 b이 m보다 작아지도록 나머지 연산
    a %= m;
    b %= m;

    if(a<b) return a+(m-b); // a가 b보다 작다면 m을 더해 양수 결과값 반환
    else return a-b;    // 그렇지 않다면 a-b
}

/*
 * mod_mul() - computes a*b mod m
 * a*b에서 오버플로가 발생할 수 있기 때문에 덧셈을 사용하여 빠르게 계산할 수 있는
 * "double addition" 알고리즘을 사용한다. 그 알고리즘은 다음과 같다.
 *     r = 0;
 *     while (b > 0) {
 *         if (b & 1)
 *             r = mod_add(r, a, m);
 *         b = b >> 1;
 *         a = mod_add(a, a, m);
 *     }
 */
uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t res = 0;

    // double addition 알고리즘을 사용해 오버플로 없이 multiply 진행
    while(b>0){ // b가 0이 될 때까지 반복
        // b의 가장 낮은 비트가 1일 때
        if(b&1) res = mod_add(res, a, m);   

        b >>= 1; // b를 오른쪽으로 한 칸 쉬프트
        a = mod_add(a, a, m);   // a는 2배가 됨
    }

    return res; // 결과값 반환
}

/*
 * mod_pow() - computes a^b mod m
 * a^b에서 오버플로가 발생할 수 있기 때문에 곱셈을 사용하여 빠르게 계산할 수 있는
 * "square multiplication" 알고리즘을 사용한다. 그 알고리즘은 다음과 같다.
 *     r = 1;
 *     while (b > 0) {
 *         if (b & 1)
 *             r = mod_mul(r, a, m);
 *         b = b >> 1;
 *         a = mod_mul(a, a, m);
 *     }
 */
uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t res = 1;

    // square multiplication 알고리즘을 사용해 오버플로 없이 지수 연산 수행
    while(b>0){ // b가 0이 될 때까지 반복
        // b의 가장 낮은 비트가 1일 때
        if(b&1) res = mod_mul(res, a, m);

        b >>= 1;    // b를 오른쪽으로 한 칸 쉬프트
        a = mod_mul(a, a, m);   // a는 제곱이 됨
    }

    return res; // 결과값 반환
 }

/*
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3,317,044,064,679,887,385,961,981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
const uint64_t a[BASELEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns PRIME if n is prime, COMPOSITE otherwise.
 */

int miller_rabin(uint64_t n)
{
    if(n==2) return PRIME;  // n이 2일 때 PRIME 반환
    else if(n%2==0 || n==1) return COMPOSITE;   // n이 1이거나 짝수일 때 COMPOSITE 반환

    uint64_t k=0, q=n-1, tmp=1;

    // k와 q를 선택함
    while(q%2==0){
        k += 1;
        q /= 2;
    }
    
    // 베이스 집합에 대해 검증
    for(int i=0;i<BASELEN;i++){
        int flag=0;

        if(a[i]>=n-1) break;    // n-1보다 a가 크면 검증할 필요 없으므로 검증 중단

        // a^q mod n 연산
        tmp = mod_pow(a[i], q, n); 

        if(tmp==1) continue;    // 연산 결과가 1이면 다음 a에 대해 검증 

        // a^(q+2^j) mod n 연산
        for(int j=0;j<k;j++){
            if(tmp==n-1){   // 연산 결과가 n-1이면 연산 중단
                flag = 1;
                break;
            } 
            tmp = mod_mul(tmp, tmp, n); // 현재 값을 mod n에 대해 제곱
        } 
        
        if(flag) continue;  // 방금 반복문을 탈출했었다면 다음 a에 대해 검증
        else return COMPOSITE;  // 0~k-1에 대해 만족하는 j가 없었다면 합성수 판정
    }

    return PRIME;   // 모든 검증을 통과했다면 소수 판정
}
