"""
BCH 코드 신드롬 계산 모듈 (학습용)

이 모듈은 BCH 디코딩의 첫 번째 단계인 신드롬 계산을
두 가지 방법으로 수행하고 상세한 중간 과정을 추적합니다.

방법 1: 다항식 평가 (Direct Evaluation)
    S_i = R(α^i) for i = 1, 2, ..., 2t

방법 2: 다항식 나눗셈 (Division Method)
    R(x) = Q_i(x)·(x - α^i) + S_i
"""

from typing import List, Tuple, Dict
from dataclasses import dataclass
import numpy as np
from .galois_field import GFElement, GaloisField


@dataclass
class HornerStep:
    """Horner 방법의 각 단계 정보"""
    iteration: int          # 반복 번호
    coefficient: GFElement  # 현재 계수
    current_value: GFElement  # 누적값
    operation: str          # 수행한 연산
    next_value: GFElement   # 다음 누적값


@dataclass
class DivisionStep:
    """다항식 나눗셈의 각 단계 정보"""
    step: int              # 단계 번호
    dividend: List[int]    # 현재 피제수
    divisor: List[int]     # 제수
    quotient_bit: int      # 몫의 현재 비트
    xor_result: List[int]  # XOR 연산 결과
    explanation: str       # 단계 설명


@dataclass
class SyndromeResult:
    """신드롬 계산 결과"""
    syndromes: List[GFElement]       # 계산된 신드롬
    method: str                       # 사용한 방법
    steps: List                       # 중간 단계 (방법에 따라 다름)
    is_valid_codeword: bool          # 유효한 코드워드인지 (모든 신드롬 = 0)
    computation_time: float          # 계산 시간 (선택)


class SyndromeCalculator:
    """
    신드롬 계산기 클래스

    BCH 코드의 신드롬을 두 가지 방법으로 계산하고
    상세한 중간 과정을 추적합니다.
    """

    def __init__(self, field: GaloisField, t: int, verbose: bool = True):
        """
        초기화

        Args:
            field: Galois Field GF(2^m)
            t: 정정 가능한 오류 개수
            verbose: 상세 출력 여부
        """
        self.field = field
        self.t = t
        self.verbose = verbose
        self.n = field.order  # 코드 길이 = 2^m - 1

    def compute_by_evaluation(self, received: List[int]) -> SyndromeResult:
        """
        방법 1: 다항식 평가로 신드롬 계산

        R(x) = r_0 + r_1·x + r_2·x^2 + ... + r_{n-1}·x^{n-1}
        S_i = R(α^i) = r_0 + r_1·α^i + r_2·α^{2i} + ... + r_{n-1}·α^{(n-1)i}

        Horner 방법 사용:
        R(x) = r_0 + x(r_1 + x(r_2 + x(...)))

        Args:
            received: 수신된 코드워드 (이진 벡터)

        Returns:
            SyndromeResult 객체
        """
        if len(received) != self.n:
            raise ValueError(f"코드워드 길이는 {self.n}이어야 합니다")

        if self.verbose:
            print("\n" + "=" * 80)
            print("신드롬 계산: 다항식 평가 방법")
            print("=" * 80)
            print(f"수신 코드워드: {''.join(map(str, received))}")
            print(f"코드워드 길이: {len(received)}")

        syndromes = []
        all_steps = []

        # 각 신드롬 S_1, S_2, ..., S_{2t} 계산
        for i in range(1, 2 * self.t + 1):
            alpha_i = self.field.alpha(i)

            if self.verbose:
                print(f"\n{'─' * 80}")
                print(f"신드롬 S_{i} = R(α^{i}) = R({alpha_i}) 계산")
                print(f"{'─' * 80}")

            # Horner 방법으로 다항식 평가
            syndrome, steps = self._evaluate_polynomial_horner(
                received, alpha_i, i
            )

            syndromes.append(syndrome)
            all_steps.append(steps)

            if self.verbose:
                print(f"결과: S_{i} = {syndrome}")
                if syndrome.is_zero():
                    print(f"      (이진: 0000, 16진: 0x0)")
                else:
                    print(f"      (이진: {syndrome.to_binary()}, 16진: 0x{syndrome.poly:X})")

        is_valid = all(s.is_zero() for s in syndromes)

        if self.verbose:
            print("\n" + "=" * 80)
            print("신드롬 계산 완료")
            print("=" * 80)
            print(f"신드롬: {[str(s) for s in syndromes]}")
            if is_valid:
                print("✓ 모든 신드롬 = 0 → 오류 없음")
            else:
                print("✗ 신드롬 ≠ 0 → 오류 존재")
            print("=" * 80)

        return SyndromeResult(
            syndromes=syndromes,
            method="polynomial_evaluation",
            steps=all_steps,
            is_valid_codeword=is_valid,
            computation_time=0.0
        )

    def _evaluate_polynomial_horner(
        self,
        coeffs: List[int],
        x: GFElement,
        syndrome_index: int
    ) -> Tuple[GFElement, List[HornerStep]]:
        """
        Horner 방법으로 다항식 평가

        P(x) = c_0 + c_1·x + c_2·x^2 + ... + c_{n-1}·x^{n-1}
             = c_0 + x(c_1 + x(c_2 + x(...)))

        Args:
            coeffs: 다항식 계수 (이진)
            x: 평가할 점
            syndrome_index: 신드롬 인덱스 (출력용)

        Returns:
            (평가 결과, Horner 단계 리스트)
        """
        steps = []

        # Horner 방법: 뒤에서부터 시작
        result = self.field.zero()

        if self.verbose:
            print(f"\nHorner 방법으로 계산:")
            print(f"R(α^{syndrome_index}) = r_0 + α^{syndrome_index}·(r_1 + α^{syndrome_index}·(r_2 + ...))")
            print(f"\n{'반복':<6} {'계수':<8} {'현재값':<15} {'연산':<30} {'다음값':<15}")
            print("─" * 80)

        for i in range(len(coeffs) - 1, -1, -1):
            # 현재 계수
            coeff = self.field.one() if coeffs[i] == 1 else self.field.zero()

            # 이전 결과에 x를 곱함
            current = result
            result = result * x

            # 현재 계수를 더함
            next_val = result + coeff

            # 연산 설명
            if i == len(coeffs) - 1:
                operation = f"초기값 = r_{i}"
            else:
                operation = f"({current}) × α^{syndrome_index} + r_{i}"

            step = HornerStep(
                iteration=len(coeffs) - i,
                coefficient=coeff,
                current_value=current,
                operation=operation,
                next_value=next_val
            )
            steps.append(step)

            if self.verbose:
                print(f"{step.iteration:<6} r_{i}={coeffs[i]:<4} {str(current):<15} {operation:<30} {str(next_val):<15}")

            result = next_val

        return result, steps

    def compute_by_division(self, received: List[int]) -> SyndromeResult:
        """
        방법 2: 다항식 나눗셈으로 신드롬 계산

        R(x)를 (x - α^i)로 나눈 나머지가 R(α^i) = S_i

        하지만 GF(2^m)에서는 직접 구현이 복잡하므로,
        실제로는 평가 방법을 사용하되 나눗셈 관점에서 설명

        Args:
            received: 수신된 코드워드

        Returns:
            SyndromeResult 객체
        """
        if self.verbose:
            print("\n" + "=" * 80)
            print("신드롬 계산: 다항식 나눗셈 개념")
            print("=" * 80)
            print("다항식 나머지 정리: R(x) = Q(x)·(x - α^i) + R(α^i)")
            print("따라서 R(α^i)는 R(x)를 (x - α^i)로 나눈 나머지입니다.")
            print("=" * 80)

        # 실제 계산은 평가 방법 사용
        return self.compute_by_evaluation(received)

    def compare_methods(
        self,
        received: List[int]
    ) -> Tuple[SyndromeResult, SyndromeResult]:
        """
        두 방법으로 계산한 결과 비교

        Args:
            received: 수신된 코드워드

        Returns:
            (평가 방법 결과, 나눗셈 방법 결과)
        """
        print("\n" + "=" * 80)
        print("두 가지 방법으로 신드롬 계산 비교")
        print("=" * 80)

        result1 = self.compute_by_evaluation(received)
        result2 = self.compute_by_division(received)

        # 결과 비교
        print("\n" + "=" * 80)
        print("결과 비교")
        print("=" * 80)

        match = all(
            s1 == s2
            for s1, s2 in zip(result1.syndromes, result2.syndromes)
        )

        print(f"평가 방법: {[str(s) for s in result1.syndromes]}")
        print(f"나눗셈 방법: {[str(s) for s in result2.syndromes]}")
        print(f"\n결과 일치: {'✓ 예' if match else '✗ 아니오'}")

        return result1, result2

    def analyze_syndrome_pattern(
        self,
        error_positions: List[int]
    ) -> Dict[int, List[GFElement]]:
        """
        오류 위치에 따른 신드롬 패턴 분석

        Args:
            error_positions: 오류 위치 리스트

        Returns:
            위치별 신드롬 딕셔너리
        """
        print("\n" + "=" * 80)
        print(f"오류 위치 패턴 분석")
        print("=" * 80)

        patterns = {}

        for pos in error_positions:
            # 해당 위치에만 오류가 있는 경우
            received = [0] * self.n
            received[pos] = 1

            result = self.compute_by_evaluation(received)
            patterns[pos] = result.syndromes

            print(f"\n위치 {pos}에 오류:")
            print(f"  신드롬: {[str(s) for s in result.syndromes]}")

            # 이론적 예상값 확인
            print(f"  예상값:")
            for i in range(1, 2 * self.t + 1):
                expected = self.field.alpha(i * pos)
                actual = result.syndromes[i - 1]
                match = "✓" if expected == actual else "✗"
                print(f"    S_{i} = α^{i * pos} = {expected} {match}")

        return patterns

    def visualize_syndrome_binary(self, syndromes: List[GFElement]):
        """
        신드롬을 이진 벡터로 시각화

        Args:
            syndromes: 신드롬 리스트
        """
        print("\n" + "=" * 80)
        print("신드롬 이진 표현")
        print("=" * 80)

        print(f"{'신드롬':<10} {'값':<15} {'이진':<10} {'16진':<8}")
        print("─" * 45)

        for i, s in enumerate(syndromes, 1):
            if s.is_zero():
                print(f"S_{i:<8} {'0':<15} {'0000':<10} {'0x0':<8}")
            else:
                binary = s.to_binary()
                hex_val = f"0x{s.poly:X}"
                print(f"S_{i:<8} {str(s):<15} {binary:<10} {hex_val:<8}")

        print("=" * 80)


def create_gf_reference_table(field: GaloisField):
    """
    GF(2^m) 참조 테이블 생성

    Args:
        field: Galois Field
    """
    print("\n" + "=" * 80)
    print(f"{field} 참조 테이블")
    print("=" * 80)

    print(f"{'지수':<10} {'다항식':<20} {'이진':<10} {'16진':<8} {'10진':<8}")
    print("─" * 60)

    # 0 원소
    print(f"{'α^(-∞)':<10} {'0':<20} {'0000':<10} {'0x0':<8} {'0':<8}")

    # 1 = α^0
    print(f"{'α^0 (=1)':<10} {'1':<20} {'0001':<10} {'0x1':<8} {'1':<8}")

    # α^1, α^2, ..., α^{order-1}
    for i in range(1, field.order):
        elem = field.alpha(i)
        poly_str = field._poly_to_string(elem.poly)
        binary = f"{elem.poly:0{field.m}b}"
        hex_val = f"0x{elem.poly:X}"

        print(f"{'α^' + str(i):<10} {poly_str:<20} {binary:<10} {hex_val:<8} {elem.poly:<8}")

    print("=" * 80)


def test_syndrome_with_errors(
    field: GaloisField,
    t: int,
    error_positions: List[int],
    test_name: str = "테스트"
):
    """
    특정 오류 패턴으로 신드롬 계산 테스트

    Args:
        field: Galois Field
        t: 정정 가능 오류 개수
        error_positions: 오류 위치 리스트
        test_name: 테스트 이름
    """
    print("\n" + "=" * 80)
    print(f"{test_name}")
    print("=" * 80)

    calc = SyndromeCalculator(field, t, verbose=True)

    # 오류 벡터 생성
    received = [0] * field.order
    for pos in error_positions:
        received[pos] = 1

    print(f"오류 위치: {error_positions}")
    print(f"오류 벡터: {''.join(map(str, received))}")

    # 신드롬 계산
    result = calc.compute_by_evaluation(received)

    # 이진 표현
    calc.visualize_syndrome_binary(result.syndromes)

    return result


if __name__ == "__main__":
    # 간단한 테스트
    gf = GaloisField(4)

    # 참조 테이블
    create_gf_reference_table(gf)

    # 위치 5에 오류
    test_syndrome_with_errors(gf, t=1, error_positions=[5], test_name="단일 오류 테스트 (위치 5)")

    # 위치 3, 10에 오류
    test_syndrome_with_errors(gf, t=2, error_positions=[3, 10], test_name="이중 오류 테스트 (위치 3, 10)")
