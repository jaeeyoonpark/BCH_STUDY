"""
Galois Field 유틸리티 모듈
galois 라이브러리를 사용하여 GF(2^m) 연산 및 유틸리티 제공
"""

import galois
import numpy as np
from typing import Union, List, Tuple


def create_gf(m: int) -> galois.FieldArray:
    """
    GF(2^m) 유한체 생성

    Args:
        m: 유한체 차수 (GF(2^m))

    Returns:
        galois.FieldArray: GF(2^m) 유한체 클래스
    """
    return galois.GF(2**m)


def print_power_table(gf: galois.FieldArray, show_all: bool = True):
    """
    GF(2^m)의 power representation 테이블 출력

    Args:
        gf: galois 유한체
        show_all: True면 모든 원소, False면 0이 아닌 원소만
    """
    m = int(np.log2(gf.order))
    n = gf.order - 1  # α의 차수 (2^m - 1)

    print(f"\nGF(2^{m}) Power Table (Primitive polynomial: {gf.irreducible_poly})")
    print("=" * 80)
    print(f"{'Power':<12} {'Polynomial':<20} {'Binary':<15} {'Decimal':<10}")
    print("-" * 80)

    # 0 원소
    if show_all:
        zero = gf(0)
        print(f"{'0':<12} {'0':<20} {'000...':<15} {int(zero):<10}")

    # α의 거듭제곱
    alpha = gf.primitive_element

    for i in range(n):
        element = alpha ** i

        # Power 표현
        if i == 0:
            power_str = "1 (α^0)"
        else:
            power_str = f"α^{i}"

        # Polynomial 표현
        poly_str = element_to_poly_str(element, var='α')

        # Binary 표현
        binary_str = format(int(element), f'0{m}b')

        # Decimal 표현
        decimal = int(element)

        print(f"{power_str:<12} {poly_str:<20} {binary_str:<15} {decimal:<10}")

    print("=" * 80)


def element_to_poly_str(element: galois.FieldArray, var: str = 'x', gf=None) -> str:
    """
    GF 원소를 다항식 문자열로 변환

    Args:
        element: galois 유한체 원소
        var: 변수명 (기본값: 'x')
        gf: galois field (optional)

    Returns:
        str: 다항식 문자열 (예: "α^2 + α + 1")

    Examples:
        >>> gf = galois.GF(2**3)
        >>> element = gf(5)  # binary: 101
        >>> element_to_poly_str(element)
        'x^2 + 1'
    """
    if int(element) == 0:
        return "0"

    # 필드 정보 가져오기
    if gf is None:
        # element로부터 필드 정보 추출
        order = element.__class__.order
    else:
        order = gf.order

    m = int(np.log2(order))
    binary = format(int(element), f'0{m}b')

    terms = []
    for i, bit in enumerate(reversed(binary)):
        if bit == '1':
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append(var)
            else:
                terms.append(f"{var}^{i}")

    return " + ".join(reversed(terms)) if terms else "0"


def poly_str_to_element(poly_str: str, gf: galois.FieldArray) -> galois.FieldArray:
    """
    다항식 문자열을 GF 원소로 변환

    Args:
        poly_str: 다항식 문자열 (예: "x^2 + x + 1")
        gf: galois 유한체

    Returns:
        galois.FieldArray: GF 원소
    """
    # 간단한 파싱 (향후 개선 가능)
    if poly_str == "0":
        return gf(0)

    # TODO: 더 robust한 파싱 구현
    raise NotImplementedError("다항식 문자열 파싱은 향후 구현 예정")


def gf_poly_eval(coeffs: List[galois.FieldArray], x: galois.FieldArray) -> galois.FieldArray:
    """
    GF 상에서 다항식 평가

    Args:
        coeffs: 계수 리스트 [c0, c1, c2, ...] (c0 + c1*x + c2*x^2 + ...)
        x: 평가할 점

    Returns:
        galois.FieldArray: 다항식 값 p(x)

    Examples:
        >>> gf = galois.GF(2**3)
        >>> alpha = gf.primitive_element
        >>> coeffs = [gf(1), alpha**2]  # 1 + α^2 * x
        >>> gf_poly_eval(coeffs, alpha)
    """
    if not coeffs:
        return type(x)(0)

    result = type(x)(0)
    x_power = type(x)(1)

    for coeff in coeffs:
        result += coeff * x_power
        x_power *= x

    return result


def gf_poly_str(coeffs: List[galois.FieldArray], var: str = 'x') -> str:
    """
    GF 다항식을 문자열로 변환

    Args:
        coeffs: 계수 리스트
        var: 변수명

    Returns:
        str: 다항식 문자열

    Examples:
        >>> gf = galois.GF(2**3)
        >>> coeffs = [gf(1), gf.primitive_element]
        >>> gf_poly_str(coeffs)
        '1 + α*x'
    """
    if not coeffs or all(int(c) == 0 for c in coeffs):
        return "0"

    terms = []
    for i, coeff in enumerate(coeffs):
        if int(coeff) == 0:
            continue

        # 계수 표현
        if int(coeff) == 1 and i > 0:
            coeff_str = ""
        else:
            coeff_str = element_to_poly_str(coeff, 'α')

        # 변수 표현
        if i == 0:
            term = coeff_str if coeff_str else "1"
        elif i == 1:
            term = f"{coeff_str}·{var}" if coeff_str else var
        else:
            term = f"{coeff_str}·{var}^{i}" if coeff_str else f"{var}^{i}"

        terms.append(term)

    return " + ".join(terms) if terms else "0"


def binary_to_gf_array(bits: List[int], gf: galois.FieldArray) -> galois.FieldArray:
    """
    이진 리스트를 GF 배열로 변환

    Args:
        bits: 이진 비트 리스트 [0, 1, 0, 1, ...]
        gf: galois 유한체

    Returns:
        galois.FieldArray: GF 배열
    """
    return gf(bits)


def gf_array_to_binary(arr: galois.FieldArray) -> List[int]:
    """
    GF 배열을 이진 리스트로 변환

    Args:
        arr: galois 배열

    Returns:
        List[int]: 이진 비트 리스트
    """
    return [int(x) for x in arr]


def get_minimal_polynomial(alpha_power: int, m: int) -> List[int]:
    """
    α^i의 최소 다항식 계산

    Args:
        alpha_power: α의 지수
        m: 유한체 차수

    Returns:
        List[int]: 최소 다항식 계수 (낮은 차수부터)
    """
    gf = galois.GF(2**m)
    alpha = gf.primitive_element
    element = alpha ** alpha_power

    # galois 라이브러리의 minimal_poly 사용
    min_poly = galois.minimal_poly(element)

    # 계수를 리스트로 변환 (낮은 차수부터)
    coeffs = [int(c) for c in min_poly.coeffs[::-1]]

    return coeffs


def get_conjugates(alpha_power: int, m: int) -> List[int]:
    """
    α^i의 conjugate class 계산

    Args:
        alpha_power: α의 지수
        m: 유한체 차수

    Returns:
        List[int]: conjugate 지수들
    """
    order = 2**m - 1
    conjugates = []
    current = alpha_power % order

    while current not in conjugates:
        conjugates.append(current)
        current = (current * 2) % order

    return sorted(conjugates)


def print_gf_element_details(element: galois.FieldArray, gf: galois.FieldArray = None):
    """
    GF 원소의 상세 정보 출력

    Args:
        element: galois 유한체 원소
        gf: galois field
    """
    if gf is None:
        order = element.__class__.order
        alpha = element.__class__.primitive_element
    else:
        order = gf.order
        alpha = gf.primitive_element

    m = int(np.log2(order))
    n = order - 1

    print(f"\nGF(2^{m}) Element Details:")
    print("-" * 60)
    print(f"  Decimal:    {int(element)}")
    print(f"  Binary:     {format(int(element), f'0{m}b')}")
    print(f"  Polynomial: {element_to_poly_str(element, 'α', gf)}")

    # α의 거듭제곱으로 표현
    if int(element) == 0:
        print(f"  Power:      0")
    else:
        for i in range(n):
            if element == alpha ** i:
                print(f"  Power:      α^{i}")
                break
    print("-" * 60)


def create_gf_reference_table(gf: galois.FieldArray):
    """
    완전한 GF(2^m) 참조 테이블 생성 (시각화 포함)

    Args:
        gf: galois 유한체
    """
    print_power_table(gf, show_all=True)

    # 추가 정보
    m = int(np.log2(gf.order))
    print(f"\nAdditional Information:")
    print(f"  Field:              GF(2^{m})")
    print(f"  Order:              {gf.order}")
    print(f"  Primitive element:  α = {gf.primitive_element}")
    print(f"  Primitive poly:     {gf.irreducible_poly}")
    print()


def alpha_power_to_element(power: int, gf: galois.FieldArray) -> galois.FieldArray:
    """
    α^i를 GF 원소로 변환

    Args:
        power: α의 지수
        gf: galois 유한체

    Returns:
        galois.FieldArray: α^power
    """
    if power == float('-inf'):
        return gf(0)

    alpha = gf.primitive_element
    n = gf.order - 1

    # 음수 지수 처리
    power = power % n if power != 0 else 0

    return alpha ** power


def element_to_alpha_power(element: galois.FieldArray, gf: galois.FieldArray = None) -> int:
    """
    GF 원소를 α의 거듭제곱으로 변환

    Args:
        element: galois 유한체 원소
        gf: galois field

    Returns:
        int: α의 지수 (0이면 -inf 반환)
    """
    if int(element) == 0:
        return float('-inf')

    if gf is None:
        alpha = element.__class__.primitive_element
        n = element.__class__.order - 1
    else:
        alpha = gf.primitive_element
        n = gf.order - 1

    for i in range(n):
        if element == alpha ** i:
            return i

    raise ValueError(f"원소 {element}는 α의 거듭제곱이 아닙니다")


# 편의 함수: 문자열 표현
def format_gf_element(element: galois.FieldArray, format_type: str = 'power', gf: galois.FieldArray = None) -> str:
    """
    GF 원소를 다양한 형식으로 포맷

    Args:
        element: galois 유한체 원소
        format_type: 'power', 'poly', 'binary', 'decimal' 중 하나
        gf: galois field

    Returns:
        str: 포맷된 문자열
    """
    if gf is None:
        order = element.__class__.order
    else:
        order = gf.order

    if format_type == 'decimal':
        return str(int(element))
    elif format_type == 'binary':
        m = int(np.log2(order))
        return format(int(element), f'0{m}b')
    elif format_type == 'poly':
        return element_to_poly_str(element, 'α', gf)
    elif format_type == 'power':
        if int(element) == 0:
            return '0'
        power = element_to_alpha_power(element, gf)
        if power == 0:
            return '1'
        return f'α^{power}'
    else:
        raise ValueError(f"Unknown format type: {format_type}")
