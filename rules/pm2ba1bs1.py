# Module for PM2, BA1, BS1
# Population data, vs normal people variants

from ..helper import *

GNOMAD_DB = "/data/projects/ACMG/database/gnomad.exomes.r2.1.1.sites.vcf.gz"


def add_gnomad_into_var_infos(variant_dic: dict, df_col2idx: dict) -> dict:
    """_summary_
    Note:
        gnomad.gz db의 용량이 매우 크므로(약 60~70gb), 다른 데이터베이스처럼 따로 저장하지 않고,
        파일의 각 줄을 읽어서 INFO 란의 AF (AC/AN) 값을 variant_infos 리스트의 gnomad_ratio 항목에 update한다.

        ##CHROM	    POS	    ID      REF	    ALT	    QUAL    FILTER	INFO
        # 1	    12198  rs62635282	G	    C	    9876.24	AC0     AC=3;AN=2000;AF=~

        먼저, gnomad db를 한 줄씩 읽어들이면서, 각 라인에 해당하는 variant_id가 환자의 variant 목록에
        있는지 확인한다. AN 값이 0인 경우, AF 값은 존재하지 않으므로, AF의 여부를 먼저 확인한 뒤, variant_id에
        해당하는 모든 variant_id-feature 에 gnomad_ac, an, af 값을 업데이트 한다. Default 값으로 None을 부여함.

    Args:
        variant_dic (dict): class VariantDF.variant_dic
        df_col2idx (dict): DF column index dictionary

    Returns:
        dict: 각각의 variant_info에 variant_id에 대한 gnomad 비율값이 update된 variant_dic

    Examples:
        >>> ['12-49500509-C-T', 'ENSG00000139636', 'ENST00000551169',
            'missense_variant,NMD_transcript_variant', '72/771', '73/120',
            '25/39', 'V/I', 'Gta/Ata', 'LMBR1L', None(revel), 0(AC), 0(AN), None(AF)]
            ->
            ['12-49500509-C-T', 'ENSG00000139636', 'ENST00000551169',
            'missense_variant,NMD_transcript_variant', '72/771', '73/120',
            '25/39', 'V/I', 'Gta/Ata', 'LMBR1L', 0.051, 240, 100000, 0.241109]
    """

    # gnomad 각 라인 파징
    for line in dbparser.read_big_gz_file(GNOMAD_DB):
        if line.startswith("#") and (not line.startswith("##")):
            f_col2idx = {
                val: idx
                for idx, val in enumerate(line.strip("#").strip().split("\t"))
            }
            continue
        if not line.startswith("#"):  # new gnomad variation
            row = line.strip().split("\t")
            chr = row[f_col2idx["CHROM"]]
            pos = row[f_col2idx["POS"]]
            ref = row[f_col2idx["REF"]]
            alt = row[f_col2idx["ALT"]]
            var_id = f"{chr}-{pos}-{ref}-{alt}"  # 1-985445-G-GT

            # 환자 변이 정보에 gnomad ratio를 할당
            if var_id in variant_dic:
                # ac, an, af 할당하는 과정. AN = 0 인 경우, AF가 없음.
                gnomad_ac = row[f_col2idx["INFO"]].split("AC=")[1].split(";")[0]
                gnomad_an = row[f_col2idx["INFO"]].split("AN=")[1].split(";")[0]
                if "AF=" in row[f_col2idx["INFO"]]:
                    gnomad_af = (
                        row[f_col2idx["INFO"]].split("AF=")[1].split(";")[0]
                    )
                else:
                    gnomad_af = None
                # 각 transcript 정보 리스트에 저장.
                for var_feature in variant_dic[var_id]:
                    if gnomad_ac:
                        variant_dic[var_id][var_feature]["var_infos"][
                            df_col2idx["gnomad_ac"]
                        ] = int(gnomad_ac)

                    if gnomad_an:
                        variant_dic[var_id][var_feature]["var_infos"][
                            df_col2idx["gnomad_an"]
                        ] = int(gnomad_an)

                    if gnomad_af:
                        variant_dic[var_id][var_feature]["var_infos"][
                            df_col2idx["gnomad_af"]
                        ] = float(gnomad_af)

    return variant_dic


def assign_ba1_rule(gnomad_an: int, gnomad_af: float) -> int:
    """_summary_
    Note:
        다음의 두 가지 조건을 만족하면 BA1 rule을 부여한다.
        1. > 1000 individuals
        2. 일반인들에게서 5% 이상 발견

    Args:
        gnomad_an (int): 특정 변이를 검사한 일반인 총 수
        gnomad_af (float): 특정 변이에 대해 일반인에서의 비율 (AC/AN)

    Returns:
        int: 1 or 0 (ba1 rule)
    """

    # 조건 1 > 1,000 individuals
    if gnomad_an < 1000:
        return 0
    elif gnomad_af >= 0.05:  # 5%
        return 1
    else:
        return 0


def compare_with_disease_inheritence(
    gnomad_an: int, gnomad_af: float, inheritence: str
) -> tuple:
    """_summary_
    Note:
        pm2: 일반인 대조군에서는 잘 발견되지 않는 변이. 'absent' or 'extremely low frequency'
        bs1: 특정한 질병에서 예상되는 빈도보다 높은 빈도의 변이

        bs1 rule의 경우, 질병에 따라 발병률이 다르므로 조건을 달리하는게 정확하지만, 각 질병의 발병률을
        구하기가 어려워서 모든 질병에 대해 일괄적인 기준을 적용하였다. ACMG guide line에서는 'absent'라고
        표현하지만, 실제로는 질병마다 조금씩 다른 기준을 적용한 경우가 많다. 따라서, '희귀질환'의 기준 및
        다른 연구자들의 기준을 참고하여 cutoff를 설정하였다.

        국내의 희귀질환 환자규모 상한: 2만명
        10만명당 환자 수 상한: 42.5명
        -> 약 50/100000 = 0.05% (0.0005)

        또한, 여러 질병에서 pm2에 대해 0.0001 ~ 0.00001 사이의 값을 cut off로 사용하고 있다.
        따라서, recessive에 대해서는 0.0001, dominant에 대해서는 그보다 낮은 0.00002를 적용하였다.
        bs1 rule에 대해서는 pm2 의 10배 값을 cut off로 적용하였다.

        <판단 기준>
        1. > 1000 individuals
        2.
            - Dominant:
                1. pm2 = 일반인에게서 0.002% 이하로 관찰 (0.00002)
                2. bs1 = 일반인 에게서 0.02% 이상으로 관찰 (pm2의 10배) (0.0002)
            - Recessive:
                1. pm2 = 일반인에게서 0.01% 이하로 관찰 (0.0001)
                2. bs1 = 일반인 에게서 0.1% 이상으로 관찰 (pm2의 10배) (0.001)

    Args:
        gnomad_an (int): 특정 변이를 검사한 일반인 총 수
        gnomad_af (float): 특정 변이에 대해 일반인에서의 비율 (AC/AN)
        inheritence (str): AD or AR or XD or XR

    Returns:
        tuple(int, int): (pm2, bs1)

    References:
        1. ClinGen Hearing Loss Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 1
        2. ClinGen Platelet Disorders Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 2.1
        3. https://varsome.com/about/resources/acmg-implementation/
        4. Overview of specifications to the ACMG/AMP variant interpretation guidelines
    """

    pm2, bs1 = 0, 0

    # > 1000 individuals
    if gnomad_an >= 1000:
        pass
    else:
        return (pm2, bs1)

    if inheritence in ["AD", "XD"]:  # dominant disease
        if gnomad_af < 0.00002:  # 0.002%
            pm2 = 1
        elif gnomad_af > 0.0002:  # 0.02%
            bs1 = 1
    elif inheritence in ["AR", "XR"]:  # recessive disease
        if gnomad_af < 0.0001:  # 0.01%
            pm2 = 1
        elif gnomad_af > 0.001:  # 0.1%
            bs1 = 1
    else:  # Y case?
        pass

    return (pm2, bs1)


def assign_pm2bs1_rule(
    gnomad_an: int,
    gnomad_af: float,
    gene_symbol: str,
    disease_db_dic: dict,
    disease_col2idx: dict,
) -> tuple:
    """_summary_
    Note:
        Gene_symbol을 키워드로 유전자-질병정보 dictionary의 정보를 확인하여, 해당 유전자의 inheritence를
        확인한다. 각 유전자에 대한 질병 정보는, 여러 보고들을 포함하기 때문에 inheritence 정보가 다른 경우가 있다.
        그 경우, 각각의 수를 세어서 dominant인지 recessive인지 비율로 기준을 정한다. 만약, X-linked 질환이
        있다면, X-linked 질환을 우선으로 판정한다. 이후, inhertence에 따라 기준을 정해서 gnomad의 AF 정보를
        바탕으로 pm2, bs1 rule을 부여한다(compare_with_disease_inheritence()).

        disease_db_dic: {
        "TBCE": [
            ["Ciliary dyskinesia, primary, 14" ,["Autosomal recessive"], ["Pediatric"]],
            ["Encephalopathy, progressive, with amyotrophy and optic atrophy" ,["Autosomal recessive"], ["Infancy", "Neonatal"]]
        ]}
        disease_col2idx: dict = {"title": 0, "inheritance": 1, "onsetAges": 2}

    Args:
        gnomad_an (int): 특정 변이를 검사한 일반인 총 수
        gnomad_af (float): 특정 변이에 대해 일반인에서의 비율 (AC/AN)
        gene_symbol (str): gene_symbol
        clinvar_db_dic (dict): gene_symbol 기준으로 정리된 clinvar_db
        clinvar_col2idx (dict): column index dictionary

    Returns:
        tuple(int, int): (pm2, bs1)

    """

    pm2, bs1 = 0, 0

    if disease_db_dic.get(gene_symbol):

        AD_count, AR_count = 0, 0
        XD_count, XR_count = 0, 0

        # 같은 유전자에 여러 질병과, 여러 inheritence가 있는 경우, 모든 AD, AR counting
        for each_disease_infos in disease_db_dic[gene_symbol]:
            AD_count += each_disease_infos[
                disease_col2idx["inheritance"]
            ].count("Autosomal dominant")

            AR_count += each_disease_infos[
                disease_col2idx["inheritance"]
            ].count("Autosomal recessive")

            XD_count += each_disease_infos[
                disease_col2idx["inheritance"]
            ].count("X-linked dominant")

            XR_count += each_disease_infos[
                disease_col2idx["inheritance"]
            ].count("X-linked recessive")

        # 보고된 비율이 더 많은 쪽을 inheritence로 설정.
        # 같은 수인 경우 dominant로 설정하며, X-linked 우선.
        if (XD_count + XR_count) > 0:
            if XD_count >= XR_count:
                inheritence = "XD"
            else:
                inheritence = "XR"
        elif AD_count == 0 and AR_count == 0:
            return (pm2, bs1)
        elif AD_count >= AR_count:
            inheritence = "AD"
        else:
            inheritence = "AR"

        # pm2, bs1 할당
        pm2, bs1 = compare_with_disease_inheritence(
            gnomad_an, gnomad_af, inheritence
        )

    return (pm2, bs1)


def execute(
    proband_var_df: object, disease_db_dic: dict, disease_col2idx: dict
) -> object:
    """_summary_
    Note: ACMG rule 중에서, pm2/ba1/bs1 에 해당하는 룰을 구현한 모듈이다. 각각의 rule에 대한 설명은 다음과 같다.

        ba1: 일반인 대조군에서 5% 이상의 비율로 발견되는 변이
        pm2: 일반인 대조군에서는 잘 발견되지 않는 변이 (dominant/recessive)
        bs1: 특정한 질병에서 예상되는 빈도보다 높은 빈도의 변이 (dominant/recessive)

        본 모듈에서는 1가지의 데이터베이스가 사용되었다.
        - disease_db: 동일한 유전자를 기준으로 하여 유전자에 해당하는 질병들을 정리하였다.
        {
            "TBCE": [
                ["Ciliary dyskinesia, primary, 14" ,["Autosomal recessive"], ["Pediatric"]],
                ["Encephalopathy, progressive, with amyotrophy and optic atrophy" ,["Autosomal recessive"], ["Infancy", "Neonatal"]]
            ]
        }
        disease_col2idx: dict = {"title": 0, "inheritance": 1, "onsetAges": 2}

        본 모듈의 작동방식은 다음과 같다.
        먼저, gnomAD database의 정보를 환자의 variant_info에 추가한다(add_gnomad_into_var_infos()).
        이 때, 변이에 따라 AF 값이 없는 경우(AN = 0), gnomad에 정보가 없는 경우, 환자의 변이 정보란에 추가가
        안되는 경우도 있다(약 2/3의 비율만 추가됨). 이 경우, 'NoneType' 으로 저장되어 있다.

        이후, 각 variant에 대하여 AF 값이 있는 경우, ba1의 할당여부를 판단한다(assign_ba1_rule()). ACMG guideline
        에서, >1000 individuals 인 경우에 적용하라고 권장하고 있으므로, 먼저 AN값이 1000 이상인지 확인하고, AF 값을
        확인한다. 만약, AF가 0.05 이상(5% 이상 발견)이라면, ba1 rule을 부여한다. BA1 rule이 부여된 경우, 중복 할당을 피하기
        위해 pm2, bs1 rule은 판단하지 않는다.

        각 변이가 일으킬 수 있는 질병의 dominant, recessive에 따라 일반인에서 발견되는 변이의 AF 판단기준을 달리하여
        pm2, bs1 rule을 부여한다(assign_pm2bs1_rule()). 이 때, gene symbol이 없으면 해당하는 disease를 확인할 수
        없으므로, gene_symbol이 있는 경우에만 본 함수를 수행한다. 본 함수에 사용한 cut off 값은 여러 레퍼런스를 참고하여
        결정하였다. 이 때, 질환마다 동일한 기준을 적용하였다. (각 함수 설명 참고)

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        disease_db_dic (dict): 각 gene_symbol에 해당하는 질병들을 정리한 dictionary. 여러개의 질병은 []로 구분.
        disease_col2idx (dict): disease_db_dic의 각 column을 index로 변환하는 dictionary

    Returns:
        object: proband_variant_DF의 variant_dic["evidnece_score_dic"]이 update된 object.

    References:
        1. ClinGen Hearing Loss Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 1
        2. ClinGen Platelet Disorders Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 2.1
        3. https://varsome.com/about/resources/acmg-implementation/
        4. Overview of specifications to the ACMG/AMP variant interpretation guidelines
    """

    df_col2idx = proband_var_df.df_col2idx
    variant_dic = proband_var_df.variant_dic
    variant_dic = add_gnomad_into_var_infos(variant_dic, df_col2idx)

    for var_id in variant_dic:
        for var_feature in variant_dic[var_id]:

            pm2, ba1, bs1 = 0, 0, 0

            gnomad_an = variant_dic[var_id][var_feature]["var_infos"][
                df_col2idx["gnomad_an"]
            ]
            gnomad_af = variant_dic[var_id][var_feature]["var_infos"][
                df_col2idx["gnomad_af"]
            ]
            gene_symbol = variant_dic[var_id][var_feature]["var_infos"][
                df_col2idx["symbol"]
            ]

            if gnomad_af:  # not NoneType

                ba1 = assign_ba1_rule(gnomad_an, gnomad_af)  # > 5%

                if gene_symbol != "-":
                    if not ba1:  # ba1이 이미 할당된 경우, pm2, bs1 not assigned
                        pm2, bs1 = assign_pm2bs1_rule(
                            gnomad_an,
                            gnomad_af,
                            gene_symbol,
                            disease_db_dic,
                            disease_col2idx,
                        )
            else:
                continue

            variant_dic[var_id][var_feature]["evidence_score_dic"]["pm2"] = pm2
            variant_dic[var_id][var_feature]["evidence_score_dic"]["ba1"] = ba1
            variant_dic[var_id][var_feature]["evidence_score_dic"]["bs1"] = bs1

    proband_var_df.variant_dic = variant_dic

    return proband_var_df
