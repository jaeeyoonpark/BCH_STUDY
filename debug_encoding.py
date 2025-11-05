#!/usr/bin/env python3
"""BCH 인코딩 및 신드롬 디버그"""

import sys
sys.path.insert(0, '/home/user/BCH_STUDY')

from bch_learning import GaloisField, BCHCode

# GF(2^4) 생성
gf = GaloisField(4)

# BCH(15, 11, 1) 생성
bch = BCHCode(m=4, t=1, field=gf)

print("=" * 60)
print("BCH 인코딩 및 신드롬 디버그")
print("=" * 60)

# 간단한 메시지: 1000...
message = [1] + [0] * 10
print(f"\n메시지: {''.join(map(str, message))}")

codeword = bch.encode(message)
print(f"코드워드: {''.join(map(str, codeword))}")
print(f"길이: {len(codeword)}")

# 코드워드가 올바른지 확인 (신드롬이 0이어야 함)
syndromes = bch.compute_syndromes(codeword)
print(f"\n코드워드 검증 (신드롬이 모두 0이어야 함):")
for i, s in enumerate(syndromes, 1):
    print(f"  S_{i} = {s}")

# 각 위치에 오류를 추가하고 신드롬 확인
print(f"\n각 위치의 오류와 신드롬:")
print(f"{'위치':<6} {'신드롬 S_1':<15} {'예상 (α^i)':<15}")
print("-" * 40)

for pos in range(min(10, len(codeword))):
    received = codeword[:]
    received[pos] ^= 1  # 비트 반전

    syndromes = bch.compute_syndromes(received)
    s1 = syndromes[0]
    expected = gf.alpha(pos)

    match = "✓" if s1 == expected else "✗"
    print(f"{pos:<6} {str(s1):<15} {str(expected):<15} {match}")

print("=" * 60)
