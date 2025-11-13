# BCH Learning 마이그레이션 가이드 (v1.x → v2.0)

## 주요 변경사항

### galois 라이브러리로 전환

v2.0부터는 자체 구현된 `GaloisField` 대신 표준 `galois` 라이브러리를 사용합니다.

### 설치

```bash
pip install galois
```

## API 변경사항

### 1. Galois Field 생성

**Before (v1.x):**
```python
from bch_learning import GaloisField

gf = GaloisField(4)  # GF(2^4)
alpha = gf.alpha(3)  # α^3
```

**After (v2.0):**
```python
from bch_learning import create_gf
import galois

gf = create_gf(4)  # GF(2^4) 생성
alpha = gf.primitive_element ** 3  # α^3
```

### 2. BCH 코드 생성

**Before (v1.x):**
```python
from bch_learning import BCHCode, GaloisField

field = GaloisField(4)
bch = BCHCode(m=4, t=1, field=field)
```

**After (v2.0):**
```python
from bch_learning import BCHCode, create_gf

gf = create_gf(4)
bch = BCHCode(m=4, t=1, gf=gf)  # 파라미터명 변경: field → gf
```

### 3. 원소 접근 및 표현

**Before (v1.x):**
```python
element = gf.alpha(5)
is_zero = element.is_zero()
str_repr = str(element)
```

**After (v2.0):**
```python
from bch_learning import format_gf_element, element_to_alpha_power

element = gf.primitive_element ** 5  # 또는 gf(5)
is_zero = int(element) == 0
power_str = format_gf_element(element, 'power')  # 'α^5'
poly_str = format_gf_element(element, 'poly')    # 'α^2 + 1'
```

## 새로운 유틸리티 함수

### Power Table 출력

```python
from bch_learning import print_power_table, create_gf

gf = create_gf(3)
print_power_table(gf, show_all=True)
```

### 완전한 참조 테이블

```python
from bch_learning import create_gf_reference_table

create_gf_reference_table(gf)
```

### 원소 포맷팅

```python
from bch_learning import format_gf_element

element = gf(5)

# 다양한 형식으로 출력
print(format_gf_element(element, 'power'))    # 'α^2'
print(format_gf_element(element, 'poly'))     # 'α^2 + 1'
print(format_gf_element(element, 'binary'))   # '101'
print(format_gf_element(element, 'decimal'))  # '5'
```

### 최소 다항식 및 Conjugate

```python
from bch_learning import get_minimal_polynomial, get_conjugates

min_poly = get_minimal_polynomial(alpha_power=3, m=4)
conjugates = get_conjugates(alpha_power=3, m=4)
```

## 호환성 노트

### 코어 알고리즘

- **Berlekamp-Massey**: API 동일, 내부만 galois 라이브러리 사용
- **Chien Search**: API 동일
- **BCH 인코딩/디코딩**: API 동일

### 타입 변경

| v1.x | v2.0 |
|------|------|
| `GFElement` | `galois.FieldArray` |
| `GaloisField` | `galois.FieldArray` (클래스) |

## 완전한 예제

```python
from bch_learning import (
    create_gf,
    BCHCode,
    decode_bch,
    print_power_table,
    format_gf_element
)

# GF(2^3) 생성
gf = create_gf(3)
print_power_table(gf)

# BCH(7,4,1) 코드
bch = BCHCode(m=3, t=1, gf=gf)

# 인코딩
message = [1, 0, 1, 1]
codeword = bch.encode(message)

# 오류 추가
received = bch.add_errors(codeword, [2])

# 디코딩
corrected, success = decode_bch(received, bch, verbose=True)

print(f"디코딩 성공: {success}")
print(f"정정된 코드워드: {corrected}")
```

## 장점

1. **표준 라이브러리**: 검증된 galois 라이브러리 사용
2. **성능**: JIT 컴파일(numba) 지원으로 빠른 연산
3. **풍부한 기능**: 최소 다항식, conjugate 계산 등 내장
4. **유틸리티 확장**: power table, 다양한 포맷팅 옵션
5. **유지보수**: 커뮤니티 지원 및 업데이트

## 문제 해결

### ImportError: No module named 'galois'

```bash
pip install galois
```

### 필드 접근 오류

galois 라이브러리의 원소는 `.field` 속성이 없습니다. 대신:

```python
# 필드 클래스 가져오기
field_class = type(element)
order = field_class.order
alpha = field_class.primitive_element
```

## 추가 리소스

- galois 라이브러리 문서: https://galois.readthedocs.io/
- BCH 코드 학습: `notebooks/` 디렉토리 참조
