# Module for PP2, BP1 rule
# Disease mechanism, missense mutation


def is_missense_var(var_infos_dic: dict, df_col2idx: dict) -> bool:
    """_summary_
    Note:
        해당 변이의 결과에 missense가 포함되어 있는지 확인하는 함수.

    Args:
        var_infos_dic (dict): variation-feature의 세부 정보 딕셔너리.
        df_col2idx (dict): column index dictionary

    Returns:
        bool: True or False
    """

    var_infos = var_infos_dic["var_infos"]  # uploaded_id

    if "missense_variant" in var_infos[df_col2idx["consequence"]]:
        return True
    else:
        return False


def check_gene_pathogenic_mechanism(
    gene_variant_dic: dict, clinvar_col2idx: dict
) -> str:
    """_summary_
    Note:
        아래와 같이 gene_symbol을 기준으로 정리된 ClinVar dictionary를 바탕으로, 해당 유전자에 생긴 변이가
        pathogenic한 지, 주된 원인이 무엇인지 판단하는 함수. 한 RCV 안에 여러 보고가 포함되어 있다. 먼저, 각
        변이의 pathogenicity와 변이의 종류의 개수를 센 뒤, 기준에 따라 원인을 규명한다.

        - gene_variant_dic:
        {
            "gene_symbol":{
                "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}]
                "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}]
            }
        }
        <병의 원인 판단 기준>
        - Missense: 아래 기준 모두 포함
            1. pathogenic missense variant 2개 이상
            2. pathogenic variant 중, pathogenic missense가 50% 이상.

        - Null variant: 아래 기준 모두 포함.
        {nonsense, frameshift variant, splice donor variant, splice acceptor variant, initiatior codon variant}
            1. pathogenic null variant 2개 이상
            2. pathogenic variant 중, pathogenic null variant가 50% 이상.

    Args:
        gene_variant_dic (dict): gene_symbol을 기준으로 정리한 ClinVar data dictionary
        clinvar_col2idx (dict): ClinVar_column index dictionary

    Returns:
        mechanism (str): "missesne" or "null variant" or ""
    """

    pathogenic_set = {"Pathogenic", "Likely pathogenic"}

    pathogenic_missense_count = 0
    pathogenic_null_count = 0
    pathogenic_count = 0

    # 각 variant_id의 pathogenicity와 variant의 종류 counting
    for var_infos in gene_variant_dic.values():
        var_patho = var_infos[clinvar_col2idx["pathogenicity"]]
        var_cons_dic = var_infos[clinvar_col2idx["consequence_dic"]]

        missense_count = var_cons_dic.get("missense variant", 0)
        null_count = (
            var_cons_dic.get("nonsense", 0)
            + var_cons_dic.get("frameshift variant", 0)
            + var_cons_dic.get("splice donor variant", 0)
            + var_cons_dic.get("splice acceptor variant", 0)
            + var_cons_dic.get("initiatior codon variant", 0)
        )
        if var_patho in pathogenic_set:  # pathogenic variant
            pathogenic_count += 1

            if (missense_count / sum(var_cons_dic.values())) > 0.5:
                pathogenic_missense_count += 1
            elif (null_count / sum(var_cons_dic.values())) > 0.5:
                pathogenic_null_count += 1

    # 기준에 따라 병의 원인 메커니즘 부여.
    if (pathogenic_missense_count >= 2) and (
        pathogenic_missense_count / pathogenic_count
    ) > 0.5:
        return "missense variant"
    elif (pathogenic_null_count >= 2) and (
        pathogenic_null_count / pathogenic_count
    ) > 0.5:
        return "null variant"

    return ""


def assign_pp2bp1_rule(
    gene_symbol: str, clinvar_db_dic: dict, clinvar_col2idx: dict
) -> tuple:
    """_summary_
    Note:
        각 variant가 missense variant인 경우, 관련된 유전자-질병의 원인을 확인하여
        "missense"면 pp2, 'null variant" 면 bp1을 부여한다.

    Args:
        gene_symbol (str): gene_symbol
        clinvar_db_dic (dict): gene_symbol 기준으로 정리된 clinvar_db
        clinvar_col2idx (dict): column index dictionary

    Returns:
        tuple(int, int): (pp2, bp1)

    """

    pp2, bp1 = 0, 0
    if clinvar_db_dic.get(gene_symbol):

        patho_mechanism: str = check_gene_pathogenic_mechanism(
            clinvar_db_dic[gene_symbol], clinvar_col2idx
        )
        if patho_mechanism == "missense variant":
            pp2 = 1
        elif patho_mechanism == "null variant":
            bp1 = 1
        else:
            pass

    return (pp2, bp1)


def execute(
    proband_var_df: object, clinvar_db_dic: dict, clinvar_col2idx: dict
) -> object:
    """_summary_
    Note: ACMG rule 중에서, pp2/bp1 에 해당하는 룰을 구현한 모듈이다. 각각의 rule에 대한 설명은 다음과 같다.

        pp2: Missense 변이가 해당 질병의 일반적인 원인이며, 주로 pathogenic 한 변이가 보고된 유전자에서의
        missense 변이
        bp1: Truncating variant (null variant)가 질병의 유일하게 알려진 원인이라고 보고된 유전자에서의
        missense 변이. pathogenic_null_count / pathogenic count == 1?

        본 모듈에서는 1가지의 데이터베이스가 사용되었다.
        - ClinVar: 동일한 유전자를 기준으로 하여 변이들을 정리하였다(RCV).
        {
            "BRCA2": {
                "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}]
                "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}]
            }
        }

        본 모듈의 작동방식은 다음과 같다.
        먼저, 각 variant에 대하여 missense variant인지 아닌지 조사한다(is_missense_var()).
        만약 해당 변이의 전사체가 missense variant라면, 해당 변이의 gene_symbol과 일치하는 ClinVar 데이터들을 조사하여,
        해당 유전자에서 발생하는 질병의 원인이 missense variant인지, null variant인지 계산한다
        (check_gene_pathogenic_mechanism()). 만약, 해당 유전자에서 주로 pathogenic 한 변이가 보고되고,
        missense 변이가 일반적인 원인이면 pp2를 부여하며, null variant가 질병의 주된 원인이라면 bp1을 부여한다
        (assign_pp2bp1_rule()).

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        clinvar_db_dic (dict): ClinVar 데이터베이스를 미리 가공한 dictionary
        clinvar_col2idx (dict): Clinvar_db_dic-index dictionary

    Returns:
        object: proband_variant_DF의 variant_dic["evidnece_score_dic"]이 update된 object.
    """

    df_col2idx = proband_var_df.df_col2idx
    variant_dic = proband_var_df.variant_dic

    for var_id in variant_dic:
        for var_feature in variant_dic[var_id]:

            pp2, bp1 = 0, 0
            if is_missense_var(variant_dic[var_id][var_feature], df_col2idx):
                gene_symbol = variant_dic[var_id][var_feature]["var_infos"][
                    df_col2idx["symbol"]
                ]
                if gene_symbol != "-":  # gene_symbol이 없는 경우, 비교 불가.
                    pp2, bp1 = assign_pp2bp1_rule(
                        gene_symbol, clinvar_db_dic, clinvar_col2idx
                    )

            variant_dic[var_id][var_feature]["evidence_score_dic"]["pp2"] = pp2
            variant_dic[var_id][var_feature]["evidence_score_dic"]["bp1"] = bp1

    proband_var_df.variant_dic = variant_dic

    return proband_var_df
