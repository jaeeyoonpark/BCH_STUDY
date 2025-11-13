"""
BCH 코드 학습을 위한 모듈 (galois 라이브러리 사용)
Berlekamp-Massey 알고리즘을 중심으로 BCH 디코딩을 학습합니다.
"""

__version__ = "2.0.0"
__author__ = "BCH Study"

# galois 라이브러리 기반 유틸리티
from .gf_utils import (
    create_gf,
    print_power_table,
    element_to_poly_str,
    gf_poly_eval,
    gf_poly_str,
    format_gf_element,
    create_gf_reference_table,
    alpha_power_to_element,
    element_to_alpha_power,
    get_minimal_polynomial,
    get_conjugates
)

# BCH 코드 및 알고리즘
from .bch_code import BCHCode
from .berlekamp_massey import BerlekampMassey
from .chien_search import ChienSearch, decode_bch

__all__ = [
    # GF 유틸리티
    'create_gf',
    'print_power_table',
    'element_to_poly_str',
    'gf_poly_eval',
    'gf_poly_str',
    'format_gf_element',
    'create_gf_reference_table',
    'alpha_power_to_element',
    'element_to_alpha_power',
    'get_minimal_polynomial',
    'get_conjugates',
    # BCH 관련
    'BCHCode',
    'BerlekampMassey',
    'ChienSearch',
    'decode_bch'
]
