"""
BCH 코드 구현

BCH (Bose-Chaudhuri-Hocquenghem) 코드는 강력한 오류 정정 코드입니다.
이 모듈은 BCH 인코딩과 신드롬 계산을 제공합니다.
"""

import numpy as np
from typing import List, Tuple
from .galois_field import GaloisField, GFElement, gf_poly_eval


class BCHCode:
    """
    BCH 코드 클래스

    BCH(n, k, t): n은 코드 길이, k는 정보 비트 수, t는 정정 가능한 오류 개수
    n = 2^m - 1 (원시 BCH 코드)

    생성 다항식 g(x)의 근은 α, α^2, ..., α^(2t)
    """

    def __init__(self, m: int, t: int, field: GaloisField = None):
        """
        BCH 코드 초기화

        Args:
            m: GF(2^m)의 차수
            t: 정정 가능한 오류 개수
            field: Galois Field (None이면 자동 생성)
        """
        self.m = m
        self.t = t
        self.n = 2 ** m - 1  # 코드 길이

        # Galois Field 생성
        if field is None:
            self.field = GaloisField(m)
        else:
            self.field = field

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
            생성 다항식의 계수 (이진 계수)
        """
        # 초기화: g(x) = 1
        g = [1]

        # 이미 처리한 conjugate 클래스 추적
        processed = set()

        # α^i를 근으로 가지는 최소 다항식을 순차적으로 곱함
        for i in range(1, 2 * self.t + 1):
            # 이미 처리된 conjugate 클래스는 건너뛰기
            if i in processed:
                continue

            # 최소 다항식 계산: (x - α^i)와 그 conjugate들
            min_poly, conjugates = self._minimal_polynomial(i)

            # 이 conjugate 클래스의 모든 원소를 처리된 것으로 표시
            processed.update(conjugates)

            # g(x) = g(x) * min_poly(x)
            g = self._poly_multiply_binary(g, min_poly)

        return g

    def _minimal_polynomial(self, power: int) -> Tuple[List[int], set]:
        """
        α^power의 최소 다항식 계산

        GF(2^m)에서 α^power의 최소 다항식은
        (x - α^power)(x - α^(power*2))(x - α^(power*4))...를 모두 곱한 것

        Returns:
            (최소 다항식의 이진 계수, conjugate 집합)
        """
        # Conjugate 집합 찾기: {power, power*2, power*4, ..., power*2^k} mod (2^m - 1)
        conjugates = set()
        current = power
        order = self.field.order

        while current not in conjugates:
            conjugates.add(current % order)
            current = (current * 2) % order

        # 최소 다항식 = ∏(x - α^i) for i in conjugates
        result = [1]  # x^0 계수
        for conj in conjugates:
            # (x - α^conj) = x + α^conj (GF(2)에서 - = +)
            # 이진 표현으로 곱셈
            alpha_poly = self.field.alpha_to_poly[conj]
            term = [alpha_poly, 1]  # α^conj + x
            result = self._poly_multiply_binary(result, term)

        return result, conjugates

    def _poly_multiply_binary(self, p1: List[int], p2: List[int]) -> List[int]:
        """
        이진 계수 다항식 곱셈 (GF(2) 계수)

        Args:
            p1, p2: 다항식 계수 리스트

        Returns:
            곱셈 결과
        """
        result = [0] * (len(p1) + len(p2) - 1)

        for i, c1 in enumerate(p1):
            for j, c2 in enumerate(p2):
                # GF(2)의 곱셈과 덧셈
                result[i + j] ^= self._gf_multiply_binary(c1, c2)

        return result

    def _gf_multiply_binary(self, a: int, b: int) -> int:
        """
        GF(2^m)에서 두 원소의 곱셈 (이진 표현)

        Args:
            a, b: GF 원소 (이진 표현)

        Returns:
            a * b
        """
        if a == 0 or b == 0:
            return 0

        # α의 지수로 변환
        if a == 1:
            return b
        if b == 1:
            return a

        power_a = self.field.poly_to_alpha.get(a, None)
        power_b = self.field.poly_to_alpha.get(b, None)

        if power_a is None or power_b is None:
            return 0

        # 지수 덧셈
        power_result = (power_a + power_b) % self.field.order

        return self.field.alpha_to_poly[power_result]

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
            raise ValueError(f"Message length must be {self.k}")

        # x^(n-k) * m(x) - 메시지를 왼쪽으로 시프트
        # 다항식 표현: m(x) * x^(n-k) = m_0·x^(n-k) + m_1·x^{n-k+1} + ...
        # 리스트 표현: [0, 0, ..., 0, m_0, m_1, ...] (앞에 n-k개의 0)
        shifted = [0] * (self.n - self.k) + message

        # parity = x^(n-k) * m(x) mod g(x)
        parity = self._poly_mod_binary(shifted, self.generator_poly)

        # 체계적 코드워드: c(x) = x^(n-k)·m(x) - r(x) (GF(2)에서 - = +)
        # 리스트 표현: [r_0, r_1, ..., r_{n-k-1}, m_0, m_1, ..., m_{k-1}]
        codeword = parity + message

        return codeword

    def _poly_mod_binary(self, dividend: List[int], divisor: List[int]) -> List[int]:
        """
        이진 다항식 나머지 연산

        Args:
            dividend: 피제수
            divisor: 제수

        Returns:
            나머지
        """
        # 복사본 생성
        result = dividend[:]
        divisor_len = len(divisor)

        for i in range(len(result) - divisor_len, -1, -1):
            if result[i + divisor_len - 1] == 1:
                for j in range(divisor_len):
                    result[i + j] ^= divisor[j]

        # 앞의 0 제거하고 나머지 반환
        remainder = result[:divisor_len - 1]
        return remainder

    def compute_syndromes(self, received: List[int]) -> List[GFElement]:
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

        syndromes = []

        # R(x)를 GF 원소로 변환
        r_poly = []
        for bit in received:
            if bit == 1:
                r_poly.append(self.field.one())
            else:
                r_poly.append(self.field.zero())

        # S_i = R(α^i) for i = 1, ..., 2t
        for i in range(1, 2 * self.t + 1):
            alpha_i = self.field.alpha(i)
            syndrome = gf_poly_eval(r_poly, alpha_i)
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

        print(" + ".join(terms[::-1]))
        print(f"이진 표현: {self.generator_poly[::-1]}")
        print(f"차수: {len(self.generator_poly) - 1}")


def syndrome_to_string(syndromes: List[GFElement]) -> str:
    """신드롬을 문자열로 변환"""
    return "[" + ", ".join(str(s) for s in syndromes) + "]"


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
