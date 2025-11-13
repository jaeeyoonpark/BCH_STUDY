"""
Berlekamp-Massey 알고리즘 구현 (galois 라이브러리 사용)

이 알고리즘은 선형 피드백 시프트 레지스터(LFSR)의 최소 길이를 찾는 알고리즘입니다.
BCH 디코딩에서는 신드롬으로부터 오류 위치 다항식(Error Locator Polynomial, ELP)을 찾는데 사용됩니다.

알고리즘의 핵심 개념:
1. Discrepancy (Δ): 현재 다항식이 신드롬을 얼마나 잘 예측하는지 측정
2. 다항식 갱신: Discrepancy가 0이 아니면 다항식을 수정
3. 길이 갱신: 필요시 다항식의 길이(L)를 증가

수학적 배경:
- 신드롬: S = [S_1, S_2, ..., S_2t]
- 오류 위치 다항식: Λ(x) = 1 + Λ_1·x + Λ_2·x^2 + ... + Λ_L·x^L
- 관계식: S_i + Λ_1·S_{i-1} + ... + Λ_L·S_{i-L} = 0 for i > L
"""

import galois
from typing import List, Tuple
from dataclasses import dataclass
from .gf_utils import gf_poly_str, format_gf_element


@dataclass
class BMIteration:
    """BM 알고리즘의 각 반복 정보"""
    iteration: int  # 반복 번호 (k)
    syndrome_used: galois.FieldArray  # 사용된 신드롬 S_k
    current_poly: List[galois.FieldArray]  # 현재 다항식 Λ^(k)
    discrepancy: galois.FieldArray  # 불일치 Δ_k
    length: int  # 현재 다항식 길이 L_k
    update_occurred: bool  # 업데이트 발생 여부
    temp_poly: List[galois.FieldArray]  # 임시 저장 다항식 B^(k)
    m_value: int  # m 값 (마지막 업데이트 이후 반복 수)
    explanation: str  # 단계 설명


class BerlekampMassey:
    """
    Berlekamp-Massey 알고리즘 클래스

    신드롬 시퀀스로부터 오류 위치 다항식(Error Locator Polynomial)을 찾습니다.
    """

    def __init__(self, syndromes: List[galois.FieldArray], verbose: bool = True):
        """
        초기화

        Args:
            syndromes: 신드롬 리스트 [S_1, S_2, ..., S_2t]
            verbose: 상세 출력 여부
        """
        self.syndromes = syndromes
        self.verbose = verbose
        # galois 라이브러리의 필드는 클래스 속성으로 접근
        if syndromes and len(syndromes) > 0:
            self.gf = type(syndromes[0])  # 필드 클래스
        else:
            self.gf = None

        # 알고리즘 상태
        self.iterations: List[BMIteration] = []
        self.final_poly: List[galois.FieldArray] = None

    def run(self) -> Tuple[List[galois.FieldArray], List[BMIteration]]:
        """
        BM 알고리즘 실행

        Returns:
            (오류 위치 다항식, 반복 정보 리스트)
        """
        if not self.syndromes:
            return [self.gf(1)], []

        # 초기화
        Lambda = [self.gf(1)]  # Λ^(0)(x) = 1
        B = [self.gf(1)]  # B^(0)(x) = 1
        L = 0  # 다항식 길이
        m = 1  # 마지막 업데이트 이후 반복 수
        b = self.gf(1)  # 마지막 업데이트 시의 discrepancy

        if self.verbose:
            print("\n" + "=" * 80)
            print("Berlekamp-Massey 알고리즘 시작")
            print("=" * 80)
            syn_strs = [format_gf_element(s, 'power') for s in self.syndromes]
            print(f"입력 신드롬: {syn_strs}")
            print(f"초기 상태: Λ(x) = 1, B(x) = 1, L = 0, m = 1")
            print("=" * 80)

        # 각 신드롬에 대해 반복
        for k in range(len(self.syndromes)):
            S_k = self.syndromes[k]

            # 1단계: Discrepancy 계산
            # Δ_k = S_{k+1} + Σ(Λ_i · S_{k+1-i}) for i=1 to L
            discrepancy = S_k
            for i in range(1, len(Lambda)):
                if k - i >= 0:
                    discrepancy = discrepancy + Lambda[i] * self.syndromes[k - i]

            if self.verbose:
                print(f"\n{'━' * 80}")
                print(f"반복 {k + 1}: 신드롬 S_{k + 1} = {format_gf_element(S_k, 'power')} 처리")
                print(f"{'━' * 80}")
                print(f"현재 다항식: Λ(x) = {gf_poly_str(Lambda)}")
                print(f"현재 길이 L = {L}, m = {m}")

            # Discrepancy 계산 과정 출력
            if self.verbose:
                print(f"\nDiscrepancy 계산:")
                print(f"  Δ_{k + 1} = S_{k + 1}", end="")
                for i in range(1, len(Lambda)):
                    if k - i >= 0 and int(Lambda[i]) != 0:
                        print(f" + Λ_{i}·S_{k + 1 - i}", end="")
                print()
                print(f"  Δ_{k + 1} = {format_gf_element(S_k, 'power')}", end="")
                for i in range(1, len(Lambda)):
                    if k - i >= 0 and int(Lambda[i]) != 0:
                        lamb_str = format_gf_element(Lambda[i], 'power')
                        syn_str = format_gf_element(self.syndromes[k - i], 'power')
                        print(f" + {lamb_str}·{syn_str}", end="")
                print()
                print(f"  Δ_{k + 1} = {format_gf_element(discrepancy, 'power')}")

            # 반복 정보 저장 (업데이트 전)
            update_occurred = False
            explanation = ""

            # 2단계: 다항식 갱신 여부 결정
            if int(discrepancy) != 0:
                # Discrepancy가 0이 아니면 업데이트 필요
                if self.verbose:
                    print(f"\nΔ ≠ 0 이므로 다항식 갱신 필요")

                # 임시 저장
                T = Lambda[:]

                # Λ(x) = Λ(x) - (Δ_k / b) · x^m · B(x)
                factor = discrepancy / b

                # x^m · B(x) 계산
                xm_B = [self.gf(0)] * m + B

                # Λ(x) 업데이트
                max_len = max(len(Lambda), len(xm_B))
                new_Lambda = Lambda + [self.gf(0)] * (max_len - len(Lambda))
                xm_B_extended = xm_B + [self.gf(0)] * (max_len - len(xm_B))

                for i in range(max_len):
                    new_Lambda[i] = new_Lambda[i] + factor * xm_B_extended[i]

                if self.verbose:
                    print(f"\n다항식 갱신:")
                    print(f"  Λ_new(x) = Λ(x) + (Δ/b)·x^{m}·B(x)")
                    disc_str = format_gf_element(discrepancy, 'power')
                    b_str = format_gf_element(b, 'power')
                    factor_str = format_gf_element(factor, 'power')
                    print(f"  여기서 Δ/b = {disc_str}/{b_str} = {factor_str}")
                    print(f"  x^{m}·B(x) = {gf_poly_str(xm_B)}")
                    print(f"  Λ_new(x) = {gf_poly_str(new_Lambda)}")

                # 길이 갱신 조건: 2L ≤ k
                if 2 * L <= k:
                    L = k + 1 - L
                    B = T
                    b = discrepancy
                    m = 1

                    explanation = f"Δ ≠ 0이고 2L ≤ k이므로, 길이 갱신: L = {L}"
                    update_occurred = True

                    if self.verbose:
                        print(f"\n길이 갱신 조건 만족 (2L ≤ k):")
                        print(f"  새 길이 L = k + 1 - L = {k + 1} - {k + 1 - L} = {L}")
                        print(f"  B(x) ← Λ_old(x) = {gf_poly_str(B)}")
                        print(f"  b ← Δ = {format_gf_element(b, 'power')}")
                        print(f"  m ← 1")
                else:
                    m += 1
                    explanation = f"Δ ≠ 0이지만 2L > k이므로, 길이 유지. m 증가: m = {m}"

                    if self.verbose:
                        print(f"\n길이 갱신 조건 불만족 (2L > k):")
                        print(f"  길이 유지: L = {L}")
                        print(f"  m 증가: m = {m}")

                Lambda = new_Lambda
            else:
                # Discrepancy가 0이면 업데이트 불필요
                m += 1
                explanation = f"Δ = 0이므로 다항식 유지. m 증가: m = {m}"

                if self.verbose:
                    print(f"\nΔ = 0 이므로 다항식 유지")
                    print(f"m 증가: m = {m}")

            # 반복 정보 저장
            iteration_info = BMIteration(
                iteration=k + 1,
                syndrome_used=S_k,
                current_poly=Lambda[:],
                discrepancy=discrepancy,
                length=L,
                update_occurred=update_occurred,
                temp_poly=B[:],
                m_value=m,
                explanation=explanation
            )
            self.iterations.append(iteration_info)

            if self.verbose:
                print(f"\n결과: {explanation}")

        # 최종 다항식
        self.final_poly = Lambda

        if self.verbose:
            print(f"\n" + "=" * 80)
            print(f"알고리즘 완료")
            print(f"=" * 80)
            print(f"최종 오류 위치 다항식: Λ(x) = {gf_poly_str(Lambda)}")
            print(f"다항식 차수: {L}")
            print(f"예상 오류 개수: {L}")
            print(f"=" * 80)

        return Lambda, self.iterations

    def get_summary_table(self) -> str:
        """
        반복 과정을 표로 요약

        Returns:
            표 형식의 문자열
        """
        if not self.iterations:
            return "알고리즘이 아직 실행되지 않았습니다."

        table = "\n" + "=" * 120 + "\n"
        table += "Berlekamp-Massey 알고리즘 반복 과정 요약\n"
        table += "=" * 120 + "\n"

        # 헤더
        header = f"{'반복':<6} {'신드롬':<12} {'Δ':<12} {'L':<4} {'m':<4} {'갱신':<6} {'Λ(x)':<50}\n"
        table += header
        table += "-" * 120 + "\n"

        # 각 반복
        for it in self.iterations:
            update_mark = "✓" if it.update_occurred else ""
            poly_str = gf_poly_str(it.current_poly)
            if len(poly_str) > 45:
                poly_str = poly_str[:42] + "..."

            syn_str = format_gf_element(it.syndrome_used, 'power')
            disc_str = format_gf_element(it.discrepancy, 'power')

            row = f"{it.iteration:<6} {syn_str:<12} {disc_str:<12} "
            row += f"{it.length:<4} {it.m_value:<4} {update_mark:<6} {poly_str:<50}\n"
            table += row

        table += "=" * 120 + "\n"

        return table

    def verify_polynomial(self) -> bool:
        """
        찾은 다항식이 신드롬 방정식을 만족하는지 검증

        Λ(x) = 1 + Λ_1·x + ... + Λ_L·x^L 일 때,
        S_i + Λ_1·S_{i-1} + ... + Λ_L·S_{i-L} = 0 for all valid i

        Returns:
            검증 성공 여부
        """
        if self.final_poly is None:
            return False

        Lambda = self.final_poly
        L = len(Lambda) - 1

        print(f"\n다항식 검증:")
        print(f"Λ(x) = {gf_poly_str(Lambda)}")
        print(f"차수 L = {L}")
        print(f"\n신드롬 방정식 확인:")

        all_valid = True

        for i in range(L, len(self.syndromes)):
            # S_{i+1} + Λ_1·S_i + Λ_2·S_{i-1} + ... + Λ_L·S_{i-L+1} = 0
            result = self.syndromes[i]
            for j in range(1, len(Lambda)):
                if i - j >= 0:
                    result = result + Lambda[j] * self.syndromes[i - j]

            is_valid = int(result) == 0
            all_valid = all_valid and is_valid

            status = "✓" if is_valid else "✗"
            result_str = format_gf_element(result, 'power')
            print(f"  i={i + 1}: S_{i + 1} + ... = {result_str} {status}")

        return all_valid


def compute_error_locator_polynomial(syndromes: List[galois.FieldArray],
                                      verbose: bool = True) -> List[galois.FieldArray]:
    """
    신드롬으로부터 오류 위치 다항식 계산 (편의 함수)

    Args:
        syndromes: 신드롬 리스트
        verbose: 상세 출력 여부

    Returns:
        오류 위치 다항식
    """
    bm = BerlekampMassey(syndromes, verbose=verbose)
    poly, _ = bm.run()
    return poly


def simple_bm_example():
    """
    간단한 BM 알고리즘 예시
    """
    from .gf_utils import create_gf, alpha_power_to_element

    print("\n간단한 Berlekamp-Massey 예시")
    print("=" * 60)

    # GF(2^4) 생성
    gf = create_gf(4)
    print(f"사용 필드: GF(2^4)")

    # 예시 신드롬: 오류 1개 케이스
    # S_1 = α^3
    alpha = gf.primitive_element
    syndromes = [alpha ** 3]

    print(f"신드롬: S_1 = α^3")

    # BM 알고리즘 실행
    bm = BerlekampMassey(syndromes, verbose=True)
    elp, iterations = bm.run()

    print(f"\n결과:")
    print(f"오류 위치 다항식: {gf_poly_str(elp)}")

    # 검증
    bm.verify_polynomial()


if __name__ == "__main__":
    simple_bm_example()
