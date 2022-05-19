# Module for PS2 (,PM6)
# De novo variants

import pytest
import __main__


@pytest.fixture
def proband_gt_dic():

    gt_dic = {
        "1-909073-C-T": "1/1",
        "1-982941-T-C": "1/1",
        "6-32151443-C-T": "1/0",
        "6-32170464-AG-A": "1/0",
        "8-76468309-GTT-G": "0/1",
        "9-132897231-G-A": "0/1",
        "15-48807637-C-T": "./1",
        "15-48807637-C-A": "1/.",
    }

    return gt_dic


@pytest.fixture
def father_gt_dic():

    gt_dic = {
        "1-909073-C-T": "0/1",
        "1-982941-T-C": "0/1",
        "6-32151443-C-T": "0/1",
    }

    return gt_dic


@pytest.fixture
def mother_gt_dic():

    gt_dic = {
        "9-132897231-G-A": "0/1",
        "15-48807637-C-T": "1/1",
    }

    return gt_dic


def execute(
    proband_var_df: object,
    proband_gt_dic: dict,
    father_gt_dic: dict,
    mother_gt_dic: dict,
):
    """_summary_
    Note: ACMG rule 중에서, ps2 에 해당하는 룰을 구현한 모듈이다. 각각의 rule에 대한 설명은 다음과 같다.

        ps2: 부모에게서 유전되지 않으 새로운 변이.

        본 모듈에서는 GT 값을 얻기 위해 따로 VCF 파일을 파징하여 얻은 dictionary를 사용하였다.
        - dict: {var_id : GT}
        {
            - 1-866319-G-A: '1/1'
            - 1-897325-G-C: '0/1'
            - 1-866511-C-CCCCT: './1'
        }

        본 모듈의 작동방식은 다음과 같다.
        먼저, 환자의 variant_id를 GT 정보를 확인한다. 다음으로, 아버지또는 어머니의 VCF 정보에서 해당 변이가 존재하는지
        확인한다. 만약 양쪽 부모에게 모두 없다면, ps2 rule을 부여한다.

        GT값을 면밀하게 분석하여, 좀 더 엄밀한 기준을 적용할 방법도 있겠지만, 본 모듈에서는 가장 간단한 방법으로 ps2를 부여하였다.

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        proband_gt_dic (dict): 환자의 {var_id : GT}
        father_gt_dic (dict): 환자 아버지의 {var_id : GT}
        mother_gt_dic (dict): 환자 어머니의 {var_id : GT}

    Returns:
        object: proband_variant_DF의 variant_dic["evidnece_score_dic"]이 update된 object.
    """
    variant_dic = proband_var_df.variant_dic

    # 부모 정보가 확인된 경우에 수행됨.
    for var_id in proband_gt_dic:

        if proband_gt_dic[var_id] == "0/1":
            # 아버지에게 해당 변이가 있는지 확인.
            if father_gt_dic.get(var_id):
                f_flag = 1
            else:
                f_flag = 0
            # 어머니에게 해당 변이가 있는지 확인.
            if mother_gt_dic.get(var_id):
                m_flag = 1
            else:
                m_flag = 0

            # 부모에게 없으면 de novo rule 부여
            if (f_flag + m_flag) == 0:
                ps2 = 1
            else:
                ps2 = 0

        elif proband_gt_dic[var_id] == "1/1":
            # 아버지에게 해당 변이가 있는지 확인.
            if father_gt_dic.get(var_id):
                f_flag = 1
            else:
                f_flag = 0
            # 어머니에게 해당 변이가 있는지 확인.
            if mother_gt_dic.get(var_id):
                m_flag = 1
            else:
                m_flag = 0

            # 부모에게 없으면 de novo rule 부여.
            if (f_flag + m_flag) == 0:
                ps2 = 1
            elif (f_flag + m_flag) == 1:  # 한쪽 부모에게만 발견, 해당 변이가 homozygous인 경우?
                ps2 = 0
                pass
            else:
                ps2 = 0

        else:  # GT = ./1 ,  1/. or X-chr
            # 아버지에게 해당 변이가 있는지 확인.
            if father_gt_dic.get(var_id):
                f_flag = 1
            else:
                f_flag = 0
            # 어머니에게 해당 변이가 있는지 확인.
            if mother_gt_dic.get(var_id):
                m_flag = 1
            else:
                m_flag = 0

            # 부모에게 없으면 de novo rule 부여
            if (f_flag + m_flag) == 0:
                ps2 = 1
            else:
                ps2 = 0

        # ps2 rule을 각 VEP annotated variant에 저장
        for var_feature in variant_dic[var_id]:
            variant_dic[var_id][var_feature]["evidence_score_dic"]["ps2"] = ps2

    proband_var_df.variant_dic = variant_dic

    return proband_var_df
