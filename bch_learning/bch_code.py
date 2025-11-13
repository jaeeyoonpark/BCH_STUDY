"""
BCH 코드 구현 (galois 라이브러리 사용)

BCH (Bose-Chaudhuri-Hocquenghem) 코드는 강력한 오류 정정 코드입니다.
이 모듈은 BCH 인코딩과 신드롬 계산을 제공합니다.
"""

import galois
import numpy as np
from typing import List, Tuple
from .gf_utils import (
    create_gf,
    get_minimal_polynomial,
    get_conjugates,
    alpha_power_to_element,
    binary_to_gf_array,
    gf_array_to_binary
)


class BCHCode:
    """
    BCH 코드 클래스

    BCH(n, k, t): n은 코드 길이, k는 정보 비트 수, t는 정정 가능한 오류 개수
    n = 2^m - 1 (원시 BCH 코드)

    생성 다항식 g(x)의 근은 α, α^2, ..., α^(2t)
    """

    def __init__(self, m: int, t: int, gf: galois.FieldArray = None):
        """
        BCH 코드 초기화

        Args:
            m: GF(2^m)의 차수
            t: 정정 가능한 오류 개수
            gf: Galois Field (None이면 자동 생성)
        """
        self.m = m
        self.t = t
        self.n = 2 ** m - 1  # 코드 길이

        # Galois Field 생성
        if gf is None:
            self.gf = create_gf(m)
        else:
            self.gf = gf

        self.alpha = self.gf.primitive_element

        # 생성 다항식 계산
        self.generator_poly = self._compute_generator_polynomial()
        self.k = self.n - len(self.generator_poly) + 1  # 정보 비트 수

        print(f"BCH({self.n}, {self.k}, {self.t}) 코드 생성")
        print(f"생성 다항식 차수: {len(self.generator_poly) - 1}")

    def _compute_generator_polynomial(self) -> List[int]:
        """
        생성 다항식 g(x) 계산

        g(x)는 α, α^2, ..., α^(2t)를 근으로 가지는 최소 다항식들의 곱

        Returns:
            생성 다항식의 계수 (이진 계수, 낮은 차수부터)
        """
        # 초기화: g(x) = 1
        g_poly = galois.Poly([1], field=galois.GF(2))

        # 이미 처리한 conjugate 클래스 추적
        processed = set()

        # α^i를 근으로 가지는 최소 다항식을 순차적으로 곱함
        for i in range(1, 2 * self.t + 1):
            # 이미 처리된 conjugate 클래스는 건너뛰기
            if i in processed:
                continue

            # Conjugate 집합 찾기
            conjugates = get_conjugates(i, self.m)
            processed.update(conjugates)

            # 최소 다항식 계산
            element = self.alpha ** i
            min_poly = element.minimal_poly()

            # g(x) = g(x) * min_poly(x)
            g_poly = g_poly * min_poly

        # 계수를 리스트로 변환 (낮은 차수부터)
        coeffs = [int(c) for c in g_poly.coeffs[::-1]]

        return coeffs

    def encode(self, message: List[int]) -> List[int]:
        """
        메시지를 BCH 코드로 인코딩

        체계적 인코딩: [parity bits | message]

        Args:
            message: 정보 비트 (길이 k)

        Returns:
            코드워드 (길이 n)
        """
        if len(message) != self.k:
            raise ValueError(f"Message length must be {self.k}, got {len(message)}")

        # galois 다항식으로 변환
        m_poly = galois.Poly(message[::-1], field=galois.GF(2))  # 높은 차수부터
        g_poly = galois.Poly(self.generator_poly[::-1], field=galois.GF(2))

        # x^(n-k) * m(x)
        shifted = m_poly * galois.Poly.Degrees([self.n - self.k], field=galois.GF(2))

        # parity = shifted mod g(x)
        _, remainder = divmod(shifted, g_poly)

        # 나머지 계수 추출 (낮은 차수부터)
        parity_coeffs = [int(c) for c in remainder.coeffs[::-1]]

        # 길이 맞추기
        parity = parity_coeffs + [0] * (self.n - self.k - len(parity_coeffs))
        parity = parity[:self.n - self.k]

        # 체계적 코드워드: [parity | message]
        codeword = parity + message

        return codeword

    def compute_syndromes(self, received: List[int]) -> List[galois.FieldArray]:
        """
        수신된 코드워드의 신드롬 계산

        신드롬 S_i = R(α^i) for i = 1, 2, ..., 2t
        여기서 R(x)는 수신된 다항식

        Args:
            received: 수신된 비트 시퀀스 (길이 n)

        Returns:
            신드롬 리스트 [S_1, S_2, ..., S_2t]
        """
        if len(received) != self.n:
            raise ValueError(f"Received vector length must be {self.n}")

        # R(x)를 galois 다항식으로 변환 (높은 차수부터)
        r_array = binary_to_gf_array(received, self.gf)
        r_poly = galois.Poly(r_array[::-1], field=self.gf)

        syndromes = []

        # S_i = R(α^i) for i = 1, ..., 2t
        for i in range(1, 2 * self.t + 1):
            alpha_i = self.alpha ** i
            syndrome = r_poly(alpha_i)
            syndromes.append(syndrome)

        return syndromes

    def add_errors(self, codeword: List[int], error_positions: List[int]) -> List[int]:
        """
        코드워드에 오류 추가 (테스트용)

        Args:
            codeword: 원본 코드워드
            error_positions: 오류 위치 리스트 (0부터 시작)

        Returns:
            오류가 추가된 코드워드
        """
        received = codeword[:]
        for pos in error_positions:
            if 0 <= pos < len(received):
                received[pos] ^= 1  # 비트 반전

        return received

    def print_generator_poly(self):
        """생성 다항식 출력"""
        print(f"\n생성 다항식 g(x):")
        terms = []
        for i, coeff in enumerate(self.generator_poly):
            if coeff == 1:
                if i == 0:
                    terms.append("1")
                elif i == 1:
                    terms.append("x")
                else:
                    terms.append(f"x^{i}")

        if terms:
            print(" + ".join(reversed(terms)))
        else:
            print("0")
        print(f"이진 표현: {self.generator_poly[::-1]}")
        print(f"차수: {len(self.generator_poly) - 1}")


def syndrome_to_string(syndromes: List[galois.FieldArray]) -> str:
    """신드롬을 문자열로 변환"""
    from .gf_utils import format_gf_element
    return "[" + ", ".join(format_gf_element(s, 'power') for s in syndromes) + "]"


def codeword_to_string(codeword: List[int]) -> str:
    """코드워드를 문자열로 변환"""
    return "".join(str(bit) for bit in codeword)


def hamming_weight(vector: List[int]) -> int:
    """해밍 가중치 계산 (1의 개수)"""
    return sum(vector)


def hamming_distance(v1: List[int], v2: List[int]) -> int:
    """해밍 거리 계산"""
    if len(v1) != len(v2):
        raise ValueError("Vectors must have the same length")
    return sum(b1 ^ b2 for b1, b2 in zip(v1, v2))
