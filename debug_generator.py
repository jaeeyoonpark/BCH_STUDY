#!/usr/bin/env python3
"""BCH 생성 다항식 디버그"""

import sys
sys.path.insert(0, '/home/user/BCH_STUDY')

from bch_learning import GaloisField, BCHCode

# GF(2^4) 생성
gf = GaloisField(4)

print("=" * 60)
print("BCH 생성 다항식 검증")
print("=" * 60)

# BCH(15, 11, 1) - α를 근으로 가져야 함
bch = BCHCode(m=4, t=1, field=gf)

bch.print_generator_poly()

# 생성 다항식 g(x)가 α, α^2를 근으로 가지는지 확인
print(f"\n생성 다항식의 근 확인:")

g_poly = bch.generator_poly
print(f"g(x) 계수: {g_poly}")

# GF 원소로 변환
g_coeffs = []
for coeff in g_poly:
    if coeff == 0:
        g_coeffs.append(gf.zero())
    else:
        # 계수는 GF(2)이므로 1이면 gf.one()
        g_coeffs.append(gf.one())

print(f"GF 계수: {[str(c) for c in g_coeffs]}")

# α와 α^2를 대입
for i in [1, 2]:
    alpha_i = gf.alpha(i)
    result = gf.zero()
    for j, coeff in enumerate(g_coeffs):
        result = result + coeff * (alpha_i ** j)

    print(f"g(α^{i}) = g({alpha_i}) = {result}", end="")
    if result.is_zero():
        print(" ✓")
    else:
        print(" ✗ (0이어야 함)")

print("=" * 60)
