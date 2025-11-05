# BCH 코드 Berlekamp-Massey 알고리즘 학습 프로젝트

BCH 코드 디코딩의 핵심인 **Berlekamp-Massey (BM) 알고리즘**을 수식적, 개념적으로 이해하기 위한 학습 프로젝트입니다.

## 프로젝트 개요

이 프로젝트는 다음을 제공합니다:

1. **GF(2^m) 유한체 구현**: BCH 코드의 수학적 기반
2. **BCH 인코더/디코더**: 완전한 BCH 코드 시스템
3. **신드롬 계산**: 두 가지 방법 (다항식 평가, 나눗셈) 상세 구현
4. **Berlekamp-Massey 알고리즘**: 상세한 단계별 구현
5. **Chien Search**: 오류 위치 탐색 알고리즘
6. **학습용 Jupyter 노트북**: 대화형 학습 환경 (신드롬 계산 + BM 알고리즘)

## 디렉토리 구조

```
BCH_STUDY/
├── bch_learning/           # 핵심 모듈
│   ├── __init__.py
│   ├── galois_field.py     # GF(2^m) 유한체 구현
│   ├── bch_code.py         # BCH 코드 및 인코딩
│   ├── syndrome_calculator.py  # 신드롬 계산 (두 가지 방법)
│   ├── berlekamp_massey.py # BM 알고리즘 (상세 버전)
│   └── chien_search.py     # Chien Search 및 완전 디코딩
├── notebooks/              # Jupyter 노트북
│   ├── syndrome_calculation_learning.ipynb  # 신드롬 계산 학습
│   └── berlekamp_massey_learning.ipynb      # BM 알고리즘 학습
├── README.md
└── requirements.txt
```

## 설치 방법

### 1. 저장소 클론

```bash
git clone <repository-url>
cd BCH_STUDY
```

### 2. 의존성 설치

```bash
pip install -r requirements.txt
```

### 3. Jupyter 노트북 실행

```bash
jupyter notebook notebooks/berlekamp_massey_learning.ipynb
```

## 사용 방법

### 1. 기본 예제: Python 스크립트

```python
from bch_learning import GaloisField, BCHCode, BerlekampMassey, ChienSearch

# GF(2^4) 생성
gf = GaloisField(4)

# BCH(15, 11, 1) 코드 생성 (1비트 오류 정정)
bch = BCHCode(m=4, t=1, field=gf)

# 메시지 인코딩
message = [1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1]
codeword = bch.encode(message)

# 오류 추가 (위치 5)
received = bch.add_errors(codeword, [5])

# 신드롬 계산
syndromes = bch.compute_syndromes(received)

# Berlekamp-Massey 알고리즘
bm = BerlekampMassey(syndromes, verbose=True)
elp, iterations = bm.run()

# Chien Search로 오류 위치 찾기
cs = ChienSearch(elp, bch.n, verbose=True)
error_positions = cs.search()

# 오류 정정
corrected = cs.correct_errors(received)
```

### 2. 완전한 디코딩 (한 번에)

```python
from bch_learning import BCHCode
from bch_learning.chien_search import decode_bch

# BCH(15, 7, 2) 코드 (2비트 오류 정정)
bch = BCHCode(m=4, t=2)

# 메시지 인코딩 및 오류 추가
message = [1, 0, 1, 1, 0, 1, 0]
codeword = bch.encode(message)
received = bch.add_errors(codeword, [3, 10])

# 완전한 디코딩 (신드롬 → BM → Chien Search)
corrected, success = decode_bch(received, bch, verbose=True)
```

### 3. Jupyter 노트북에서 학습

Jupyter 노트북(`berlekamp_massey_learning.ipynb`)은 다음을 포함합니다:

1. **배경 이론**
   - BCH 코드 개요
   - GF(2^m) 유한체 기초
   - BCH 디코딩 3단계

2. **Berlekamp-Massey 알고리즘 상세 설명**
   - 수학적 원리
   - Discrepancy 개념
   - 다항식 갱신 메커니즘
   - 의사코드

3. **단계별 예시**
   - **예시 1**: 오류 1개 (t=1)
   - **예시 2**: 오류 2개 (t=2) - 각 반복 상세 분석
   - **예시 3**: 완전한 디코딩 과정

4. **실험 및 시각화**
   - 다양한 오류 패턴 테스트
   - BM 알고리즘 수렴 시각화
   - GF 크기 비교

## 핵심 개념

### BCH 디코딩 3단계

```
수신 코드워드 R(x)
    ↓
[1단계] 신드롬 계산
    S_i = R(α^i) for i = 1, 2, ..., 2t
    ↓
[2단계] Berlekamp-Massey 알고리즘
    신드롬 → 오류 위치 다항식 Λ(x)
    ↓
[3단계] Chien Search
    Λ(x)의 근 찾기 → 오류 위치
    ↓
오류 정정
```

### Berlekamp-Massey 알고리즘

**목적**: 신드롬 시퀀스를 생성하는 최소 길이의 선형 피드백 시프트 레지스터(LFSR) 찾기

**핵심 요소**:
1. **Discrepancy (Δ)**: 현재 다항식의 예측 정확도
2. **다항식 갱신**: Δ ≠ 0일 때 다항식 수정
3. **길이 갱신**: 2L ≤ k 조건으로 다항식 길이 증가

**수식**:
```
Δ_k = S_k + Σ(Λ_i · S_{k-i})
Λ^(k)(x) = Λ^(k-1)(x) + (Δ_k/b) · x^m · B^(k-1)(x)
```

## 주요 특징

### 1. 상세한 디버깅 출력

BM 알고리즘 실행 시 각 단계를 상세히 출력:
- 현재 다항식 상태
- Discrepancy 계산 과정
- 다항식 갱신 이유
- 길이 변화

### 2. 검증 기능

- 신드롬 방정식 검증
- 오류 정정 결과 검증
- 각 단계별 중간 결과 확인

### 3. 교육적 설계

- 명확한 변수명과 주석
- 수학적 표기와 코드의 일대일 대응
- 단계별 진행 상황 추적

## 예시 출력

### Berlekamp-Massey 알고리즘 실행 예시

```
================================================================================
Berlekamp-Massey 알고리즘 시작
================================================================================
입력 신드롬: ['α^7', 'α^14']
초기 상태: Λ(x) = 1, B(x) = 1, L = 0, m = 1
================================================================================

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
반복 1: 신드롬 S_1 = α^7 처리
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
현재 다항식: Λ(x) = 1
현재 길이 L = 0, m = 1

Discrepancy 계산:
  Δ_1 = S_1
  Δ_1 = α^7
  Δ_1 = α^7

Δ ≠ 0 이므로 다항식 갱신 필요

다항식 갱신:
  Λ_new(x) = Λ(x) + (Δ/b)·x^1·B(x)
  여기서 Δ/b = α^7/1 = α^7
  x^1·B(x) = α^0·x
  Λ_new(x) = 1 + α^7·x

길이 갱신 조건 만족 (2L ≤ k):
  새 길이 L = k + 1 - L = 1 - 0 = 1
  B(x) ← Λ_old(x) = 1
  b ← Δ = α^7
  m ← 1

결과: Δ ≠ 0이고 2L ≤ k이므로, 길이 갱신: L = 1
...
```

## 학습 목표

이 프로젝트를 통해 다음을 이해할 수 있습니다:

1. **GF(2^m) 유한체**의 구조와 연산
2. **BCH 코드**의 생성과 인코딩 과정
3. **신드롬 계산**의 의미와 방법
4. **Berlekamp-Massey 알고리즘**의 작동 원리:
   - Discrepancy의 수학적 의미
   - 다항식 갱신 메커니즘
   - 길이 갱신 조건의 의미
5. **Chien Search**의 효율적 구현
6. 전체 **BCH 디코딩 흐름**

## 테스트

간단한 테스트 실행:

```bash
python -c "from bch_learning.berlekamp_massey import simple_bm_example; simple_bm_example()"
python -c "from bch_learning.chien_search import simple_chien_search_example; simple_chien_search_example()"
```

## 참고 자료

1. **교재**
   - "Error Control Coding" by Shu Lin and Daniel J. Costello
   - "Algebraic Codes for Data Transmission" by Richard E. Blahut

2. **논문**
   - Berlekamp, E. R. (1968). "Algebraic Coding Theory"
   - Massey, J. L. (1969). "Shift-register synthesis and BCH decoding"

3. **온라인 리소스**
   - Wikipedia: BCH code, Berlekamp-Massey algorithm
   - MIT OpenCourseWare: Coding Theory

## 라이센스

이 프로젝트는 교육 목적으로 제작되었습니다.

## 기여

버그 리포트나 개선 제안은 언제나 환영합니다!

## 문의

질문이나 피드백이 있으시면 이슈를 등록해주세요.
