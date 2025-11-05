"""
Chien Search 알고리즘 구현

Chien Search는 오류 위치 다항식(Error Locator Polynomial)의 근을 찾는 알고리즘입니다.
Λ(x)의 근이 α^(-i)이면, 위치 i에 오류가 있습니다.

알고리즘 원리:
- Λ(x) = 1 + Λ_1·x + Λ_2·x^2 + ... + Λ_L·x^L
- GF(2^m)의 모든 원소 α^(-i) (i = 0, 1, ..., n-1)에 대해 Λ(α^(-i))를 계산
- Λ(α^(-i)) = 0이면 위치 i에 오류 존재

효율적 구현:
- 각 반복에서 α를 곱하면서 다항식을 평가 (Horner's method 변형)
"""

from typing import List, Tuple
from .galois_field import GFElement, gf_poly_eval, gf_poly_str


class ChienSearch:
    """
    Chien Search 알고리즘 클래스

    오류 위치 다항식의 근을 찾아 오류 위치를 결정합니다.
    """

    def __init__(self, error_locator_poly: List[GFElement],
                 code_length: int, verbose: bool = True):
        """
        초기화

        Args:
            error_locator_poly: 오류 위치 다항식 Λ(x)
            code_length: 코드 길이 n
            verbose: 상세 출력 여부
        """
        self.elp = error_locator_poly
        self.n = code_length
        self.verbose = verbose
        self.field = error_locator_poly[0].field if error_locator_poly else None

        self.error_positions = []

    def search(self) -> List[int]:
        """
        오류 위치 탐색

        Returns:
            오류 위치 리스트 (0부터 시작)
        """
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
            alpha_inv_i = self.field.alpha(-i)

            # Λ(α^(-i)) 계산
            result = gf_poly_eval(self.elp, alpha_inv_i)

            if self.verbose:
                print(f"\n위치 {i}: Λ(α^{-i}) = Λ({alpha_inv_i}) = {result}", end="")

            # 근이면 오류 위치
            if result.is_zero():
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

    def efficient_search(self) -> List[int]:
        """
        효율적인 Chien Search 구현

        각 반복에서 이전 값에 α를 곱하여 계산량을 줄입니다.

        Returns:
            오류 위치 리스트
        """
        if self.verbose:
            print("\n" + "=" * 80)
            print("Chien Search (효율적 버전)")
            print("=" * 80)

        error_positions = []

        # 초기값: Λ_i (i = 0, 1, ..., L)
        components = self.elp[:]

        # 각 위치에 대해
        for i in range(self.n):
            # Λ(α^(-i)) = Σ Λ_j · (α^(-i))^j
            result = self.field.zero()
            for comp in components:
                result = result + comp

            if self.verbose and i < 10:  # 처음 10개만 출력
                print(f"위치 {i}: ", end="")
                print(f"Λ(α^{-i}) = {result}", end="")

            if result.is_zero():
                error_positions.append(i)
                if self.verbose and i < 10:
                    print(" ← 오류!")
            else:
                if self.verbose and i < 10:
                    print()

            # 다음 반복을 위해 각 성분에 α^j를 곱함
            for j in range(1, len(components)):
                components[j] = components[j] * self.field.alpha(j)

        if self.verbose and self.n > 10:
            print(f"... ({self.n - 10}개 더 검사)")

        self.error_positions = error_positions

        if self.verbose:
            print("\n" + "=" * 80)
            print(f"탐색 완료: {len(error_positions)}개의 오류 발견")
            print(f"오류 위치: {error_positions}")
            print("=" * 80)

        return error_positions

    def correct_errors(self, received: List[int]) -> List[int]:
        """
        오류 정정

        Args:
            received: 수신된 코드워드

        Returns:
            정정된 코드워드
        """
        if len(received) != self.n:
            raise ValueError(f"Received codeword length must be {self.n}")

        corrected = received[:]

        for pos in self.error_positions:
            # GF(2)에서는 비트 반전
            corrected[pos] ^= 1

        if self.verbose:
            print(f"\n오류 정정 완료:")
            print(f"  정정된 위치: {self.error_positions}")
            print(f"  정정 전: {''.join(map(str, received))}")
            print(f"  정정 후: {''.join(map(str, corrected))}")

        return corrected


def find_error_positions(error_locator_poly: List[GFElement],
                          code_length: int,
                          verbose: bool = True) -> List[int]:
    """
    오류 위치 찾기 (편의 함수)

    Args:
        error_locator_poly: 오류 위치 다항식
        code_length: 코드 길이
        verbose: 상세 출력 여부

    Returns:
        오류 위치 리스트
    """
    cs = ChienSearch(error_locator_poly, code_length, verbose=verbose)
    return cs.search()


def decode_bch(received: List[int], bch_code, verbose: bool = True) -> Tuple[List[int], bool]:
    """
    완전한 BCH 디코딩 (신드롬 → BM → Chien Search)

    Args:
        received: 수신된 코드워드
        bch_code: BCHCode 객체
        verbose: 상세 출력 여부

    Returns:
        (정정된 코드워드, 정정 성공 여부)
    """
    from .berlekamp_massey import BerlekampMassey

    if verbose:
        print("\n" + "=" * 80)
        print("BCH 디코딩 시작")
        print("=" * 80)
        print(f"수신: {''.join(map(str, received))}")

    # 1단계: 신드롬 계산
    if verbose:
        print("\n[1단계] 신드롬 계산")
        print("-" * 80)

    syndromes = bch_code.compute_syndromes(received)

    if verbose:
        print(f"신드롬: {[str(s) for s in syndromes]}")

    # 모든 신드롬이 0이면 오류 없음
    if all(s.is_zero() for s in syndromes):
        if verbose:
            print("모든 신드롬이 0 → 오류 없음")
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

    # 오류 개수 확인
    num_errors = len(error_positions)
    expected_errors = len(elp) - 1

    if num_errors != expected_errors:
        if verbose:
            print(f"\n경고: 찾은 오류 개수({num_errors})가 예상({expected_errors})과 다릅니다.")
        return received, False

    # 오류 정정
    corrected = cs.correct_errors(received)

    # 검증: 정정된 코드워드의 신드롬 확인
    syndromes_after = bch_code.compute_syndromes(corrected)
    success = all(s.is_zero() for s in syndromes_after)

    if verbose:
        print("\n" + "=" * 80)
        if success:
            print("디코딩 성공!")
        else:
            print("디코딩 실패 - 신드롬이 0이 아님")
        print("=" * 80)

    return corrected, success


def simple_chien_search_example():
    """
    간단한 Chien Search 예시
    """
    from .galois_field import GaloisField

    print("\n간단한 Chien Search 예시")
    print("=" * 60)

    # GF(2^4) 생성
    gf = GaloisField(4)

    # 오류 위치 다항식: Λ(x) = 1 + α^3·x
    # 근은 α^(-3) = α^12 (α^15 = 1이므로)
    # 따라서 위치 12에 오류
    elp = [gf.one(), gf.alpha(3)]

    print(f"오류 위치 다항식: Λ(x) = {gf_poly_str(elp)}")

    # Chien Search
    cs = ChienSearch(elp, code_length=15, verbose=True)
    error_positions = cs.search()

    print(f"\n결과: 오류 위치 = {error_positions}")


if __name__ == "__main__":
    simple_chien_search_example()
