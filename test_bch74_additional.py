#!/usr/bin/env python3
"""
BCH(7,4) 추가 예제 검증 - 위치 2와 위치 0
"""

import sys
sys.path.insert(0, '/home/user/BCH_STUDY')

from bch_learning import GaloisField, BerlekampMassey, ChienSearch
from bch_learning.syndrome_calculator import SyndromeCalculator

gf = GaloisField(3)

print("=" * 80)
print("BCH(7,4) 추가 예제")
print("=" * 80)

# 예제 2: 위치 2에 에러
print("\n" + "=" * 80)
print("예제 2: 위치 2에 에러")
print("=" * 80)

received2 = [0, 0, 1, 0, 0, 0, 0]  # r₂ = 1
print(f"수신: {''.join(map(str, received2))}")

calc = SyndromeCalculator(gf, t=1, verbose=False)
result2 = calc.compute_by_evaluation(received2)

print(f"\n신드롬:")
print(f"  S₁ = {result2.syndromes[0]} (예상: α²)")
print(f"  일치: {result2.syndromes[0] == gf.alpha(2)}")

bm2 = BerlekampMassey(result2.syndromes, verbose=False)
elp2, _ = bm2.run()

print(f"\nELP: Λ(x) = 1 + {elp2[1]}·x")
print(f"예상: Λ(x) = 1 + α²·x")

cs2 = ChienSearch(elp2, 7, verbose=False)

# 근 찾기
print(f"\nChien Search:")
for j in range(7):
    val = elp2[0] + elp2[1] * gf.alpha(j)
    if val.is_zero():
        print(f"  근: α^{j}")
        print(f"  에러 위치 = 7 - {j} = {7-j} ✓")

error_pos2 = cs2.search()
print(f"\n찾은 오류 위치: {error_pos2}")
print(f"예상: [2]")
print(f"일치: {error_pos2 == [2]}")

# 예제 3: 위치 0에 에러
print("\n" + "=" * 80)
print("예제 3: 위치 0에 에러")
print("=" * 80)

received0 = [1, 0, 0, 0, 0, 0, 0]  # r₀ = 1
print(f"수신: {''.join(map(str, received0))}")

result0 = calc.compute_by_evaluation(received0)

print(f"\n신드롬:")
print(f"  S₁ = {result0.syndromes[0]} (예상: 1 = α⁰)")
print(f"  일치: {result0.syndromes[0] == gf.one()}")

bm0 = BerlekampMassey(result0.syndromes, verbose=False)
elp0, _ = bm0.run()

print(f"\nELP: Λ(x) = 1 + {elp0[1]}·x")
print(f"예상: Λ(x) = 1 + x")

cs0 = ChienSearch(elp0, 7, verbose=False)

# 근 찾기
print(f"\nChien Search:")
for j in range(7):
    val = elp0[0] + elp0[1] * gf.alpha(j)
    if val.is_zero():
        print(f"  근: α^{j}")
        print(f"  에러 위치 = 7 - {j} mod 7 = {(7-j)%7} ✓")

error_pos0 = cs0.search()
print(f"\n찾은 오류 위치: {error_pos0}")
print(f"예상: [0]")
print(f"일치: {error_pos0 == [0]}")

print("\n" + "=" * 80)
print("모든 예제 검증 완료!")
print("=" * 80)
