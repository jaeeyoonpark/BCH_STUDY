"""
BCH 코드 학습을 위한 모듈
Berlekamp-Massey 알고리즘을 중심으로 BCH 디코딩을 학습합니다.
"""

__version__ = "1.1.0"
__author__ = "BCH Study"

from .galois_field import GaloisField, GFElement
from .bch_code import BCHCode
from .berlekamp_massey import BerlekampMassey
from .chien_search import ChienSearch
from .syndrome_calculator import SyndromeCalculator

__all__ = [
    'GaloisField',
    'GFElement',
    'BCHCode',
    'BerlekampMassey',
    'ChienSearch',
    'SyndromeCalculator'
]
