#!/usr/bin/env python3
"""
BCH 코드 및 Berlekamp-Massey 알고리즘 기본 테스트
"""

import sys
sys.path.insert(0, '/home/user/BCH_STUDY')

from bch_learning import GaloisField, BCHCode, BerlekampMassey, ChienSearch
from bch_learning.chien_search import decode_bch


def test_galois_field():
    """GF(2^4) 테스트"""
    print("=" * 60)
    print("테스트 1: Galois Field GF(2^4)")
    print("=" * 60)

    gf = GaloisField(4)
    print(f"✓ {gf} 생성 성공")
    print(f"  원소 개수: {gf.size}")
    print(f"  원시 다항식: {bin(gf.primitive_poly)}")

    # 기본 연산 테스트
    a = gf.alpha(3)
    b = gf.alpha(5)
    c = a * b

    print(f"\n✓ 연산 테스트:")
    print(f"  {a} × {b} = {c}")
    print(f"  예상: α^8")
    assert c.power == 8, "곱셈 오류"
    print(f"  ✓ 통과")

    return True


def test_bch_code_single_error():
    """BCH 코드 단일 오류 테스트"""
    print("\n" + "=" * 60)
    print("테스트 2: BCH(15,11,1) - 단일 오류")
    print("=" * 60)

    gf = GaloisField(4)
    bch = BCHCode(m=4, t=1, field=gf)

    print(f"✓ BCH({bch.n}, {bch.k}, {bch.t}) 생성 성공")

    # 메시지 인코딩
    message = [1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1]
    codeword = bch.encode(message)
    print(f"\n✓ 인코딩 성공")
    print(f"  메시지:    {''.join(map(str, message))}")
    print(f"  코드워드:  {''.join(map(str, codeword))}")

    # 오류 추가
    error_pos = [5]
    received = bch.add_errors(codeword, error_pos)
    print(f"\n✓ 오류 추가 (위치 {error_pos})")
    print(f"  수신:      {''.join(map(str, received))}")

    # 신드롬 계산
    syndromes = bch.compute_syndromes(received)
    print(f"\n✓ 신드롬 계산 성공")
    print(f"  S_1 = {syndromes[0]}, S_2 = {syndromes[1]}")

    # BM 알고리즘
    print(f"\n✓ Berlekamp-Massey 실행...")
    bm = BerlekampMassey(syndromes, verbose=False)
    elp, _ = bm.run()
    print(f"  오류 위치 다항식: {[str(c) for c in elp]}")

    # Chien Search
    print(f"\n✓ Chien Search 실행...")
    cs = ChienSearch(elp, bch.n, verbose=False)
    found_errors = cs.search()
    print(f"  찾은 오류 위치: {found_errors}")

    # 검증
    assert found_errors == error_pos, "오류 위치 불일치"
    print(f"  ✓ 정확히 찾음!")

    # 오류 정정
    corrected = cs.correct_errors(received)
    assert corrected == codeword, "정정 실패"
    print(f"\n✓ 오류 정정 성공!")
    print(f"  정정 후: {''.join(map(str, corrected))}")
    print(f"  원본:    {''.join(map(str, codeword))}")

    return True


def test_bch_code_double_error():
    """BCH 코드 이중 오류 테스트"""
    print("\n" + "=" * 60)
    print("테스트 3: BCH(15,7,2) - 이중 오류")
    print("=" * 60)

    gf = GaloisField(4)
    bch = BCHCode(m=4, t=2, field=gf)

    print(f"✓ BCH({bch.n}, {bch.k}, {bch.t}) 생성 성공")

    # 메시지 인코딩
    message = [1, 0, 1, 1, 0, 1, 0]
    codeword = bch.encode(message)
    print(f"\n✓ 인코딩 성공")

    # 오류 추가
    error_pos = [3, 10]
    received = bch.add_errors(codeword, error_pos)
    print(f"✓ 오류 추가 (위치 {error_pos})")

    # 완전한 디코딩
    print(f"\n✓ 완전한 디코딩 실행...")
    corrected, success = decode_bch(received, bch, verbose=False)

    assert success, "디코딩 실패"
    assert corrected == codeword, "정정 실패"

    print(f"  ✓ 디코딩 성공!")
    print(f"  수신: {''.join(map(str, received))}")
    print(f"  정정: {''.join(map(str, corrected))}")
    print(f"  원본: {''.join(map(str, codeword))}")

    return True


def test_bm_iterations():
    """BM 알고리즘 반복 과정 테스트"""
    print("\n" + "=" * 60)
    print("테스트 4: BM 알고리즘 반복 과정")
    print("=" * 60)

    gf = GaloisField(4)

    # 간단한 신드롬
    syndromes = [gf.alpha(3), gf.alpha(6)]

    print(f"✓ 신드롬: {[str(s) for s in syndromes]}")

    bm = BerlekampMassey(syndromes, verbose=False)
    elp, iterations = bm.run()

    print(f"✓ {len(iterations)}번 반복 완료")
    print(f"✓ 최종 다항식: {[str(c) for c in elp]}")

    # 반복 정보 확인
    assert len(iterations) == 2, "반복 횟수 오류"
    print(f"✓ 반복 정보 수집 성공")

    # 요약 테이블
    print("\n" + bm.get_summary_table())

    return True


def main():
    """모든 테스트 실행"""
    print("\n" + "=" * 60)
    print("BCH 코드 및 Berlekamp-Massey 알고리즘 테스트")
    print("=" * 60 + "\n")

    tests = [
        ("Galois Field", test_galois_field),
        ("BCH 단일 오류", test_bch_code_single_error),
        ("BCH 이중 오류", test_bch_code_double_error),
        ("BM 반복 과정", test_bm_iterations),
    ]

    results = []

    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
            print(f"\n{'✓' * 30}")
            print(f"✓ {name} 테스트 통과!")
            print(f"{'✓' * 30}")
        except Exception as e:
            results.append((name, False))
            print(f"\n{'✗' * 30}")
            print(f"✗ {name} 테스트 실패!")
            print(f"✗ 오류: {e}")
            print(f"{'✗' * 30}")
            import traceback
            traceback.print_exc()

    # 최종 결과
    print("\n" + "=" * 60)
    print("테스트 결과 요약")
    print("=" * 60)

    for name, success in results:
        status = "✓ 통과" if success else "✗ 실패"
        print(f"{status}: {name}")

    passed = sum(1 for _, s in results if s)
    total = len(results)

    print("\n" + "=" * 60)
    print(f"전체: {passed}/{total} 통과")
    print("=" * 60 + "\n")

    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
