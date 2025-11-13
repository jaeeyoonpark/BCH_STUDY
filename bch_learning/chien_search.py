"""
Chien Search 알고리즘 구현 (galois 라이브러리 사용)

Chien Search는 오류 위치 다항식(Error Locator Polynomial)의 근을 찾는 알고리즘입니다.
Λ(x)의 근이 α^(-i)이면, 위치 i에 오류가 있습니다.
"""

import galois
from typing import List, Tuple
from .gf_utils import gf_poly_eval, gf_poly_str, format_gf_element


class ChienSearch:
    """Chien Search 알고리즘 클래스"""

    def __init__(self, error_locator_poly: List[galois.FieldArray],
                 code_length: int, verbose: bool = True):
        self.elp = error_locator_poly
        self.n = code_length
        self.verbose = verbose
        # galois 라이브러리의 필드는 클래스 속성으로 접근
        if error_locator_poly and len(error_locator_poly) > 0:
            self.gf = type(error_locator_poly[0])  # 필드 클래스
            self.alpha = self.gf.primitive_element
        else:
            self.gf = None
            self.alpha = None
        self.error_positions = []

    def search(self) -> List[int]:
        """오류 위치 탐색"""
        if self.verbose:
            print("\n" + "=" * 80)
            print("Chien Search 알고리즘 시작")
            print("=" * 80)
            print(f"오류 위치 다항식: Λ(x) = {gf_poly_str(self.elp)}")
            print(f"코드 길이: n = {self.n}")
            print("=" * 80)

        error_positions = []

        # 모든 위치에 대해 검사
        for i in range(self.n):
            # α^(-i) 계산
            alpha_inv_i = self.alpha ** (-i)

            # Λ(α^(-i)) 계산
            result = gf_poly_eval(self.elp, alpha_inv_i)

            if self.verbose:
                alpha_str = format_gf_element(alpha_inv_i, 'power')
                result_str = format_gf_element(result, 'power')
                print(f"\n위치 {i}: Λ(α^{-i}) = Λ({alpha_str}) = {result_str}", end="")

            # 근이면 오류 위치
            if int(result) == 0:
                error_positions.append(i)
                if self.verbose:
                    print(" ← 오류 발견!")
            else:
                if self.verbose:
                    print()

        self.error_positions = error_positions

        if self.verbose:
            print("\n" + "=" * 80)
            print(f"탐색 완료: {len(error_positions)}개의 오류 발견")
            print(f"오류 위치: {error_positions}")
            print("=" * 80)

        return error_positions

    def correct_errors(self, received: List[int]) -> List[int]:
        """오류 정정"""
        if len(received) != self.n:
            raise ValueError(f"Received codeword length must be {self.n}")

        corrected = received[:]

        for pos in self.error_positions:
            corrected[pos] ^= 1  # GF(2)에서는 비트 반전

        if self.verbose:
            print(f"\n오류 정정 완료:")
            print(f"  정정된 위치: {self.error_positions}")
            print(f"  정정 전: {''.join(map(str, received))}")
            print(f"  정정 후: {''.join(map(str, corrected))}")

        return corrected


def find_error_positions(error_locator_poly: List[galois.FieldArray],
                          code_length: int,
                          verbose: bool = True) -> List[int]:
    """오류 위치 찾기 (편의 함수)"""
    cs = ChienSearch(error_locator_poly, code_length, verbose=verbose)
    return cs.search()


def decode_bch(received: List[int], bch_code, verbose: bool = True) -> Tuple[List[int], bool]:
    """완전한 BCH 디코딩 (신드롬 → BM → Chien Search)"""
    from .berlekamp_massey import BerlekampMassey

    if verbose:
        print("\n" + "=" * 80)
        print("완전한 BCH 디코딩")
        print("=" * 80)

    # 1단계: 신드롬 계산
    if verbose:
        print("\n[1단계] 신드롬 계산")
        print("-" * 80)

    syndromes = bch_code.compute_syndromes(received)

    if verbose:
        for i, s in enumerate(syndromes, 1):
            print(f"  S_{i} = {format_gf_element(s, 'power')}")

    # 신드롬이 모두 0이면 오류 없음
    if all(int(s) == 0 for s in syndromes):
        if verbose:
            print("\n모든 신드롬이 0입니다. 오류가 없습니다.")
        return received, True

    # 2단계: Berlekamp-Massey 알고리즘
    if verbose:
        print("\n[2단계] Berlekamp-Massey 알고리즘")
        print("-" * 80)

    bm = BerlekampMassey(syndromes, verbose=verbose)
    elp, _ = bm.run()

    # 3단계: Chien Search
    if verbose:
        print("\n[3단계] Chien Search")
        print("-" * 80)

    cs = ChienSearch(elp, bch_code.n, verbose=verbose)
    error_positions = cs.search()

    # 4단계: 오류 정정
    corrected = cs.correct_errors(received)

    success = len(error_positions) <= bch_code.t

    if verbose:
        print("\n" + "=" * 80)
        print(f"디코딩 {'성공' if success else '실패'}")
        if not success:
            print(f"  오류 개수 ({len(error_positions)})가 정정 능력 (t={bch_code.t})을 초과했습니다.")
        print("=" * 80)

    return corrected, success
