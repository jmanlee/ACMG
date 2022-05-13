# Module for PVS1
# 질병의 원인으로 알려진 유전자의 Null variants (nonsense, frameshift, splicing site, exon deletion etc..)


def is_null_var(var_infos_dic: dict, df_col2idx: dict) -> bool:
    """_summary_
    Note:
        해당 변이의 결과에 null variant가 포함되어 있는지 확인하는 함수.

    Args:
        var_infos_dic (dict): variation-feature의 세부 정보 딕셔너리.
        df_col2idx (dict): column index dictionary

    Returns:
        bool: True or False
    """

    var_infos = var_infos_dic["var_infos"]  # uploaded_id
    # splice 관련 변이는 제외함.
    if ("stop_gained" or "start_lost" or "frameshift_variant") in var_infos[
        df_col2idx["consequence"]
    ]:
        return True
    else:
        return False


def check_disease_cause_is_null(
    gene_variant_dic: dict, clinvar_col2idx: dict
) -> str:
    """_summary_
    Note:
        아래와 같이 gene_symbol을 기준으로 정리된 ClinVar dictionary를 바탕으로, 해당 유전자에 생긴 변이가
        pathogenic한 지, 주된 원인이 null variant인지 판단하는 함수. 한 RCV 안에 여러 보고가 포함되어 있다.
        먼저, 각 변이의 pathogenicity와 변이의 종류의 개수를 센 뒤, 기준에 따라 원인을 규명한다.

        - gene_variant_dic:
        {
            "gene_symbol":{
                "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}]
                "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}]
            }
        }

        <병의 원인 판단 기준>
        - Null variant: 아래 기준 모두 포함.
        {nonsense, frameshift variant, splice donor variant, splice acceptor variant, initiatior codon variant}
            1. 보고된 pathogenic null variant 2개 이상
            2. pathogenic variant 중, pathogenic null variant가 50% 이상.

    Args:
        gene_variant_dic (dict): gene_symbol을 기준으로 정리한 ClinVar data dictionary
        clinvar_col2idx (dict): ClinVar_column index dictionary

    Returns:
        mechanism (bool): True or False (null variant)
    """

    pathogenic_set = {"Pathogenic", "Likely pathogenic"}

    pathogenic_null_count = 0
    pathogenic_count = 0

    # 각 variant_id의 pathogenicity와 variant의 종류 counting
    for var_infos in gene_variant_dic.values():
        var_patho = var_infos[clinvar_col2idx["pathogenicity"]]
        var_cons_dic = var_infos[clinvar_col2idx["consequence_dic"]]

        null_count = (
            var_cons_dic.get("nonsense", 0)
            + var_cons_dic.get("frameshift variant", 0)
            + var_cons_dic.get("splice donor variant", 0)
            + var_cons_dic.get("splice acceptor variant", 0)
            + var_cons_dic.get("initiatior codon variant", 0)
        )
        if var_patho in pathogenic_set:  # pathogenic variant
            pathogenic_count += 1

            if (null_count / sum(var_cons_dic.values())) > 0.5:
                pathogenic_null_count += 1

    # 기준에 따라 병의 원인 메커니즘이 null variant 인지 아닌지 판단.
    if (pathogenic_null_count >= 2) and (
        pathogenic_null_count / pathogenic_count
    ) > 0.5:
        return True
    else:
        return False


def check_disease_inheritence(gene_disease_infos: list, disease_col2idx: dict):

    pass

    return

    # dominant 우선으로 // 수정 필요
    if (
        "Autosomal dominant"
        in gene_disease_infos[disease_col2idx["inheritance"]]
    ):
        return "Autosomal dominant"
    elif (
        "Autosomal recessive"
        in gene_disease_infos[disease_col2idx["inheritance"]]
    ):
        return "Autosomal recessive"


def execute(
    proband_var_df: object,
    clinvar_db_dic: dict,
    clinvar_col2idx: dict,
    disease_db_dic: dict,
    disease_col2idx: dict,
) -> object:
    """_summary_
    Note: ACMG rule 중에서, pvs1 에 해당하는 룰을 구현한 모듈이다. 각각의 rule에 대한 설명은 다음과 같다.

        pvs1: Null variant가 질병의 원인으로 알려진 유전자에서 발견된 null variant
        1. Null variant: nonsense, frameshift, canonical ±1 or 2 splice sites,
            initiation codon, single exon or multiexon deletion. splice관련은 제외함.
        2. NMD: 마지막 exon 또는 penultimate exon의 끝 50bp에서 stop codon이 발생했을 때, NMD X.
        3. Splice_site: require functional assay.
        4. whether variant remove 10% of protein.

        본 모듈에서는 2가지의 데이터베이스가 사용되었다. (disease data에 사용에 대한 수정필요.)
        - ClinVar: 동일한 유전자를 기준으로 하여 변이들을 정리하였다(RCV).
        {
            "BRCA2": {
                "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}]
                "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}]
            }
        }
        - Disease: Gene symbol을 key로 하여 질병의 정보를 정리하였다.
        {
            "TBCE": [
                ["Ciliary dyskinesia, primary, 14" ,["Autosomal recessive"], ["Pediatric"]],
                ["Encephalopathy, progressive, with amyotrophy and optic atrophy" ,["Autosomal recessive"], ["Infancy", "Neonatal"]]
            ]
        }

        본 모듈의 작동방식은 다음과 같다.
        먼저, 각 variant에 대하여 변이가 null variant인지 아닌지 조사한다(is_null_var()).
        만약 해당 변이의 전사체가 vep에 의해 null variant라고 판정되었다면, 해당 변이의 gene_symbol과
        일치하는 ClinVar 데이터들을 조사하여, 해당 유전자에서 발생하는 질병의 원인이 null variant인지 계산한다
        (check_disease_cause_is_null()). 만약, 해당 유전자에서 주로 pathogenic 한 null 변이가 보고었다면,
        유전형태를 파악하고(반영x), protein에서 변이가 생긴 위치에 따라 (기준 90%) pvs1 을 차등부여한다.

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        clinvar_db_dic (dict): ClinVar 데이터베이스를 미리 가공한 dictionary
        clinvar_col2idx (dict): Clinvar_db_dic-index dictionary
        disease_db_dic (dict): Disease 데이터베이스를 미리 가공한 dictionary
        disease_col2idx (dict): Disease_db_list-index dictionary

    Returns:
        object: proband_variant_DF의 variant_dic["evidnece_score_dic"]이 update된 object.

    References:
        1. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion
          (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6185798/)
    """

    df_col2idx = proband_var_df.df_col2idx
    variant_dic = proband_var_df.variant_dic

    for var_id in variant_dic:
        for var_feature in variant_dic[var_id]:

            pvs1 = 0
            if is_null_var(variant_dic[var_id][var_feature], df_col2idx):
                gene_symbol = variant_dic[var_id][var_feature]["var_infos"][
                    df_col2idx["symbol"]
                ]
                if gene_symbol != "-":  # gene_symbol이 없는 경우, 비교 불가.
                    if clinvar_db_dic.get(gene_symbol):
                        if check_disease_cause_is_null(
                            clinvar_db_dic[gene_symbol], clinvar_col2idx
                        ):
                            check_disease_inheritence(  # 추후 필요하면 사용
                                disease_db_dic[gene_symbol], disease_col2idx
                            )

                            # check whether variant remove 10% of protein.
                            protein_pos = variant_dic[var_id][var_feature][
                                "var_infos"
                            ][df_col2idx["protein_pos"]]
                            # e.g. 77/97, 987-988/4911
                            var_pos, protein_len = protein_pos.split("/")
                            var_pos = var_pos.split("-")[0]

                            # variant remove 기준: 90%
                            try:
                                if (int(var_pos) / int(protein_len)) < 0.9:
                                    pvs1 = 1
                                else:
                                    pvs1 = 1
                            except TypeError:  # '-'
                                pvs1 = 1

            variant_dic[var_id][var_feature]["evidence_score_dic"][
                "pvs1"
            ] = pvs1

    proband_var_df.variant_dic = variant_dic

    return proband_var_df
