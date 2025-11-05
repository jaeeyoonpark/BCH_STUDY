"""
Galois Field GF(2^m) 구현

이 모듈은 BCH 코드에 필요한 유한체 연산을 제공합니다.
GF(2^m)는 2^m개의 원소를 가지는 유한체로, 각 원소는 m비트 이진수로 표현됩니다.
"""

import numpy as np
from typing import List, Union


class GFElement:
    """
    Galois Field의 원소를 나타내는 클래스

    GF(2^m)의 원소는 다항식으로 표현됩니다:
    a_{m-1} * α^{m-1} + ... + a_1 * α + a_0
    여기서 a_i ∈ {0, 1}

    내부적으로는 정수(power representation) 또는 다항식 계수로 저장됩니다.
    """

    def __init__(self, value: Union[int, List[int]], field: 'GaloisField',
                 is_power: bool = False):
        """
        GF 원소 초기화

        Args:
            value: 정수 값 또는 다항식 계수 리스트
            field: 속한 Galois Field
            is_power: True면 α의 지수로 해석, False면 다항식 계수로 해석
        """
        self.field = field

        if is_power:
            # α^value 형태로 주어진 경우
            if value == -np.inf or value is None:
                self.power = -np.inf  # 0 원소
                self.poly = 0
            else:
                # 음수 지수도 modulo 연산으로 처리
                value = value % field.order
                self.power = value
                self.poly = field.alpha_to_poly[value]
        else:
            # 다항식 계수로 주어진 경우
            if isinstance(value, list):
                self.poly = sum(bit << i for i, bit in enumerate(value))
            else:
                self.poly = value

            # 0인 경우 특별 처리
            if self.poly == 0:
                self.power = -np.inf
            else:
                self.power = field.poly_to_alpha.get(self.poly, None)
                if self.power is None:
                    raise ValueError(f"Invalid polynomial {self.poly} for this field")

    def __add__(self, other: 'GFElement') -> 'GFElement':
        """GF 덧셈: 다항식 계수의 XOR 연산"""
        if not isinstance(other, GFElement):
            raise TypeError("Can only add GFElement objects")
        if self.field != other.field:
            raise ValueError("Elements must be from the same field")

        result_poly = self.poly ^ other.poly  # XOR 연산 (GF(2)에서 덧셈)
        return GFElement(result_poly, self.field, is_power=False)

    def __mul__(self, other: 'GFElement') -> 'GFElement':
        """GF 곱셈: 지수의 덧셈 (modulo order)"""
        if not isinstance(other, GFElement):
            raise TypeError("Can only multiply GFElement objects")
        if self.field != other.field:
            raise ValueError("Elements must be from the same field")

        # 0과의 곱셈
        if self.power == -np.inf or other.power == -np.inf:
            return GFElement(0, self.field, is_power=False)

        # α^a * α^b = α^(a+b mod order)
        result_power = (self.power + other.power) % self.field.order
        return GFElement(result_power, self.field, is_power=True)

    def __truediv__(self, other: 'GFElement') -> 'GFElement':
        """GF 나눗셈: 지수의 뺄셈"""
        if not isinstance(other, GFElement):
            raise TypeError("Can only divide GFElement objects")
        if self.field != other.field:
            raise ValueError("Elements must be from the same field")
        if other.power == -np.inf:
            raise ZeroDivisionError("Division by zero in GF")

        if self.power == -np.inf:
            return GFElement(0, self.field, is_power=False)

        result_power = (self.power - other.power) % self.field.order
        return GFElement(result_power, self.field, is_power=True)

    def __pow__(self, exp: int) -> 'GFElement':
        """GF 거듭제곱: 지수의 곱셈"""
        if self.power == -np.inf:
            return GFElement(0, self.field, is_power=False)

        result_power = (self.power * exp) % self.field.order
        return GFElement(result_power, self.field, is_power=True)

    def inverse(self) -> 'GFElement':
        """GF 역원: α^(-a) = α^(order - a)"""
        if self.power == -np.inf:
            raise ZeroDivisionError("Zero has no inverse")

        result_power = (-self.power) % self.field.order
        return GFElement(result_power, self.field, is_power=True)

    def __eq__(self, other) -> bool:
        """같음 비교"""
        if not isinstance(other, GFElement):
            return False
        return self.poly == other.poly and self.field == other.field

    def __repr__(self) -> str:
        """문자열 표현"""
        if self.power == -np.inf:
            return "0"
        elif self.power == 0:
            return "1"
        else:
            return f"α^{int(self.power)}"

    def __str__(self) -> str:
        return self.__repr__()

    def to_binary(self) -> str:
        """이진수 표현"""
        return f"{self.poly:0{self.field.m}b}"

    def is_zero(self) -> bool:
        """0인지 확인"""
        return self.power == -np.inf


class GaloisField:
    """
    Galois Field GF(2^m) 클래스

    원시 다항식(primitive polynomial)을 사용하여 GF(2^m)를 생성합니다.
    모든 0이 아닌 원소는 원시원소 α의 거듭제곱으로 표현됩니다.
    """

    # 일반적인 원시 다항식들
    PRIMITIVE_POLYNOMIALS = {
        2: 0b111,      # x^2 + x + 1
        3: 0b1011,     # x^3 + x + 1
        4: 0b10011,    # x^4 + x + 1
        5: 0b100101,   # x^5 + x^2 + 1
        6: 0b1000011,  # x^6 + x + 1
        7: 0b10001001, # x^7 + x^3 + 1
        8: 0b100011101 # x^8 + x^4 + x^3 + x^2 + 1
    }

    def __init__(self, m: int, primitive_poly: int = None):
        """
        GF(2^m) 초기화

        Args:
            m: 체의 차수 (원소 개수는 2^m개)
            primitive_poly: 원시 다항식 (주어지지 않으면 기본값 사용)
        """
        self.m = m
        self.size = 2 ** m  # 체의 크기
        self.order = self.size - 1  # α의 위수 (0 제외)

        # 원시 다항식 설정
        if primitive_poly is None:
            if m in self.PRIMITIVE_POLYNOMIALS:
                self.primitive_poly = self.PRIMITIVE_POLYNOMIALS[m]
            else:
                raise ValueError(f"No default primitive polynomial for m={m}")
        else:
            self.primitive_poly = primitive_poly

        # α^i -> 다항식 표현 변환 테이블 생성
        self.alpha_to_poly = {}
        self.poly_to_alpha = {}

        self._build_tables()

    def _build_tables(self):
        """
        α^i와 다항식 표현 간의 변환 테이블 구축

        α^0 = 1
        α^1 = α (이진수로 0010)
        α^2 = α^2 (이진수로 0100)
        ...
        원시 다항식을 사용하여 α^m 이상은 축약
        """
        # α^0 = 1
        poly = 1
        self.alpha_to_poly[0] = poly
        self.poly_to_alpha[poly] = 0

        # α^1, α^2, ..., α^(2^m - 2) 계산
        for i in range(1, self.order):
            # 다음 거듭제곱: α를 곱함 (왼쪽 시프트)
            poly <<= 1

            # 최상위 비트가 넘어가면 원시 다항식으로 모듈로 연산
            if poly & (1 << self.m):
                poly ^= self.primitive_poly

            self.alpha_to_poly[i] = poly
            self.poly_to_alpha[poly] = i

        # 0 원소 (α^(-∞))
        self.alpha_to_poly[-np.inf] = 0

    def element(self, value: Union[int, List[int]], is_power: bool = False) -> GFElement:
        """
        GF 원소 생성

        Args:
            value: 정수 값 또는 다항식 계수 리스트
            is_power: True면 α의 지수, False면 다항식 계수

        Returns:
            GFElement 객체
        """
        return GFElement(value, self, is_power)

    def zero(self) -> GFElement:
        """0 원소 반환"""
        return GFElement(0, self, is_power=False)

    def one(self) -> GFElement:
        """1 원소 반환"""
        return GFElement(0, self, is_power=True)

    def alpha(self, power: int = 1) -> GFElement:
        """α^power 반환"""
        return GFElement(power, self, is_power=True)

    def random_element(self) -> GFElement:
        """랜덤 GF 원소 생성"""
        poly = np.random.randint(0, self.size)
        return GFElement(poly, self, is_power=False)

    def all_elements(self) -> List[GFElement]:
        """모든 GF 원소 반환"""
        elements = [self.zero()]
        for i in range(self.order):
            elements.append(self.alpha(i))
        return elements

    def __eq__(self, other) -> bool:
        """필드 비교"""
        if not isinstance(other, GaloisField):
            return False
        return self.m == other.m and self.primitive_poly == other.primitive_poly

    def __repr__(self) -> str:
        return f"GF(2^{self.m})"

    def print_table(self):
        """GF 테이블 출력 (학습용)"""
        print(f"\n{self} 원소 테이블:")
        print(f"원시 다항식: {bin(self.primitive_poly)}")
        print("-" * 60)
        print(f"{'Power':<10} {'Polynomial':<15} {'Binary':<10} {'Decimal':<10}")
        print("-" * 60)

        print(f"{'α^(-∞)':<10} {'0':<15} {0:0{self.m}b:<10} {0:<10}")

        for i in range(self.order):
            poly = self.alpha_to_poly[i]
            # 다항식 문자열 표현
            poly_str = self._poly_to_string(poly)
            print(f"{'α^' + str(i):<10} {poly_str:<15} {poly:0{self.m}b:<10} {poly:<10}")
        print("-" * 60)

    def _poly_to_string(self, poly: int) -> str:
        """다항식을 문자열로 변환"""
        if poly == 0:
            return "0"
        if poly == 1:
            return "1"

        terms = []
        for i in range(self.m - 1, -1, -1):
            if poly & (1 << i):
                if i == 0:
                    terms.append("1")
                elif i == 1:
                    terms.append("α")
                else:
                    terms.append(f"α^{i}")

        return " + ".join(terms) if terms else "0"


def gf_poly_eval(coeffs: List[GFElement], x: GFElement) -> GFElement:
    """
    다항식 평가: P(x) = c_0 + c_1*x + c_2*x^2 + ...

    Args:
        coeffs: 다항식 계수 리스트 [c_0, c_1, c_2, ...]
        x: 평가할 점

    Returns:
        P(x) 값
    """
    if not coeffs:
        return x.field.zero()

    # Horner's method
    result = coeffs[-1]
    for i in range(len(coeffs) - 2, -1, -1):
        result = result * x + coeffs[i]

    return result


def gf_poly_multiply(p1: List[GFElement], p2: List[GFElement]) -> List[GFElement]:
    """
    두 다항식의 곱셈

    Args:
        p1: 첫 번째 다항식 계수
        p2: 두 번째 다항식 계수

    Returns:
        곱셈 결과 다항식 계수
    """
    if not p1 or not p2:
        return [p1[0].field.zero()]

    field = p1[0].field
    result = [field.zero()] * (len(p1) + len(p2) - 1)

    for i, c1 in enumerate(p1):
        for j, c2 in enumerate(p2):
            result[i + j] = result[i + j] + c1 * c2

    return result


def gf_poly_str(coeffs: List[GFElement]) -> str:
    """다항식을 문자열로 변환"""
    if not coeffs or all(c.is_zero() for c in coeffs):
        return "0"

    terms = []
    for i, coeff in enumerate(coeffs):
        if not coeff.is_zero():
            if i == 0:
                terms.append(str(coeff))
            elif i == 1:
                terms.append(f"{coeff}·x")
            else:
                terms.append(f"{coeff}·x^{i}")

    return " + ".join(terms) if terms else "0"
