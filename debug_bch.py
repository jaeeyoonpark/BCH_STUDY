#!/usr/bin/env python3
"""BCH 코드 디버그"""

import sys
sys.path.insert(0, '/home/user/BCH_STUDY')

from bch_learning import GaloisField, BCHCode, BerlekampMassey, ChienSearch

# GF(2^4) 생성
gf = GaloisField(4)

# BCH(15, 11, 1) 생성
bch = BCHCode(m=4, t=1, field=gf)

print("=" * 60)
print("BCH 코드 디버그")
print("=" * 60)

# 간단한 테스트: 모든 0인 코드워드
print("\n테스트 1: 모든 0인 코드워드")
codeword = [0] * 15
print(f"코드워드: {''.join(map(str, codeword))}")

syndromes = bch.compute_syndromes(codeword)
print(f"신드롬: {[str(s) for s in syndromes]}")
print(f"모두 0? {all(s.is_zero() for s in syndromes)}")

# 위치 5에 오류 추가
print("\n테스트 2: 위치 5에 오류")
received = codeword[:]
received[5] = 1
print(f"수신: {''.join(map(str, received))}")
print(f"위치 5의 값: {received[5]}")

syndromes = bch.compute_syndromes(received)
print(f"\n신드롬:")
for i, s in enumerate(syndromes, 1):
    print(f"  S_{i} = {s} (power={s.power if not s.is_zero() else '-∞'})")

# BM 알고리즘
print("\nBerlekamp-Massey:")
bm = BerlekampMassey(syndromes, verbose=False)
elp, _ = bm.run()
print(f"오류 위치 다항식: {[str(c) for c in elp]}")

# Chien Search 상세
print("\nChien Search 상세:")
for i in range(15):
    alpha_inv = gf.alpha(-i)
    # Λ(α^(-i)) 계산
    result = gf.zero()
    for j, coeff in enumerate(elp):
        result = result + coeff * (alpha_inv ** j)

    print(f"위치 {i}: Λ(α^{-i}) = Λ({alpha_inv}) = {result}", end="")
    if result.is_zero():
        print(" ← 오류!")
    else:
        print()

print("\n" + "=" * 60)
