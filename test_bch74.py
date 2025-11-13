#!/usr/bin/env python3
"""
BCH(7,4) 예제 검증

제공된 예제와 현재 구현이 일치하는지 확인합니다.
"""

import sys
sys.path.insert(0, '/home/user/BCH_STUDY')

from bch_learning import GaloisField, BCHCode, BerlekampMassey, ChienSearch
from bch_learning.syndrome_calculator import SyndromeCalculator

print("=" * 80)
print("BCH(7,4) 예제 검증")
print("=" * 80)

# GF(2^3) 생성 - primitive polynomial: x^3 + x + 1 = 0b1011
gf = GaloisField(3)
print(f"\n생성된 유한체: {gf}")
print(f"원시 다항식: {bin(gf.primitive_poly)} (x³ + x + 1)")

# Power table 출력
gf.print_table()

# BCH(7, 4, 1) 코드
print("\n" + "=" * 80)
print("BCH(7, 4, 1) 코드 생성")
print("=" * 80)

bch = BCHCode(m=3, t=1, field=gf)
print(f"실제 k = {bch.k} (예상: 4)")

# 예제: 위치 5에 오류
print("\n" + "=" * 80)
print("예제: 위치 5에 오류")
print("=" * 80)

# All-zero 코드워드 (또는 임의의 유효한 코드워드)
received = [0, 0, 0, 0, 0, 1, 0]  # r = (r₀, r₁, r₂, r₃, r₄, r₅, r₆)
print(f"수신: {''.join(map(str, received))}")
print(f"인덱싱: r = (r₀, r₁, r₂, r₃, r₄, r₅, r₆)")
print(f"오류 위치: 5 (r₅ = 1)")

# 신드롬 계산
print("\n" + "-" * 80)
print("신드롬 계산")
print("-" * 80)

calc = SyndromeCalculator(gf, t=1, verbose=False)
result = calc.compute_by_evaluation(received)

print(f"S₁ = {result.syndromes[0]}")
print(f"S₂ = {result.syndromes[1]}")

# 예상값 확인
expected_s1 = gf.alpha(5)  # α⁵
print(f"\n예상 S₁ = α⁵ = {expected_s1}")
print(f"일치: {result.syndromes[0] == expected_s1}")

# S₂ = (S₁)² 확인
s2_calc = result.syndromes[0] ** 2
print(f"\nS₂ = (S₁)² = (α⁵)² = α¹⁰ = {s2_calc}")
print(f"실제 S₂ = {result.syndromes[1]}")
print(f"일치: {result.syndromes[1] == s2_calc}")

# BM 알고리즘
print("\n" + "-" * 80)
print("Berlekamp-Massey 알고리즘")
print("-" * 80)

bm = BerlekampMassey(result.syndromes, verbose=False)
elp, iterations = bm.run()

print(f"오류 위치 다항식: Λ(x) = {'+'.join([str(c) + ('·x' if i==1 else '') for i, c in enumerate(elp) if not c.is_zero()])}")
print(f"\n예상: Λ(x) = 1 + α⁵·x")

# 계수 확인
print(f"\nΛ₀ = {elp[0]} (예상: 1)")
print(f"Λ₁ = {elp[1]} (예상: α⁵)")

# Chien Search
print("\n" + "-" * 80)
print("Chien Search")
print("-" * 80)

cs = ChienSearch(elp, bch.n, verbose=False)

# 각 α^j를 대입
print(f"\nELP의 근 찾기:")
print(f"{'α^j':<8} {'Λ(α^j)':<15} {'결과':<10}")
print("-" * 35)

for j in range(7):
    alpha_j = gf.alpha(j)
    # Λ(α^j) 계산
    val = elp[0] + elp[1] * alpha_j
    is_root = val.is_zero()

    print(f"α^{j:<6} {str(val):<15} {'근 ✓' if is_root else ''}")

# 실제 Chien Search 실행
error_positions = cs.search()

print(f"\n찾은 오류 위치: {error_positions}")
print(f"예상 오류 위치: [5]")
print(f"일치: {error_positions == [5]}")

# 공식 확인
print("\n" + "-" * 80)
print("공식 확인: 근 α^j → 오류 위치 = n - j")
print("-" * 80)

print("\n예제에서:")
print("  근: α² (j=2)")
print("  오류 위치 = 7 - 2 = 5 ✓")

print("\n우리 구현:")
print("  위치 i에 대해 α^(-i) 대입")
print("  위치 5: α^(-5) = α^(7-5) = α²")
print("  Λ(α²) = 0 → 위치 5가 오류 ✓")

print("\n" + "=" * 80)
print("검증 완료!")
print("=" * 80)
