# module for calculating the possibility of pathogenicity
# Based on bayesian rule


def calculate_acmg(proband_var_df: object, disease_db_dic: dict):

    """_summary_
    Note:
        각 환자의 variant의 판별된 ACMG rule을 bayesian score로 환산한 뒤, 
        별도의 outfile을 제작한다. 

        이 때, BA1 만을 가진 경우는 제외하였다. 

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        disease_db_dic (dict): Disease에 대한 정보를 유전자 별로 정리한 dict.
    """

    df_col2idx = proband_var_df.df_col2idx
    variant_dic = proband_var_df.variant_dic

    Opvs = 350
    prior_p = 0.1

    with open(
        "/data/projects/ACMG/output/ACMG_result_0520_2.txt", "w"
    ) as outfile:
        # contents
        print(
            "#Variant_ID",
            "Feature",
            "Consequence",
            "Symbol",
            "ACMG_rule",
            "ACMG_bayesian",
            "Pathogenicity",
            "Related disease_info",
            sep="\t",
            file=outfile,
        )

        # counting
        for var_id in variant_dic:
            for var_feature in variant_dic[var_id]:

                pp, pm, ps, pvs = 0, 0, 0, 0
                bp, bs = 0, 0
                ACMG_rules = []
                diseases = []

                # 각각의 항목에 대한 rule의 수를 counting
                for k, v in variant_dic[var_id][var_feature][
                    "evidence_score_dic"
                ].items():
                    if k.startswith("pvs") and v:
                        pvs += 1
                        ACMG_rules.append(v)  # flag_type
                    elif k.startswith("ps") and v:
                        ps += 1
                        ACMG_rules.append(k.upper())
                    elif k.startswith("pm") and v:
                        pm += 1
                        ACMG_rules.append(k.upper())
                    elif k.startswith("pp") and v:
                        pp += 1
                        ACMG_rules.append(k.upper())
                    elif k.startswith("bs") and v:
                        bs += 1
                        ACMG_rules.append(k.upper())
                    elif k.startswith("bp") and v:
                        bp += 1
                        ACMG_rules.append(k.upper())
                    elif k.startswith("ba") and v:
                        ACMG_rules.append(k.upper())

                ACMG_rule = "||".join(ACMG_rules)

                if ACMG_rule == "BA1":
                    continue

                var_infos = variant_dic[var_id][var_feature]["var_infos"]
                consequence = var_infos[df_col2idx["consequence"]]
                gene_symbol = var_infos[df_col2idx["symbol"]]

                # Bayesian framework
                OP = Opvs ** (
                    pp / 8 + pm / 4 + ps / 2 + pvs / 1 - bp / 8 - bs / 2
                )
                post_p = (OP * prior_p) / ((OP - 1) * prior_p + 1)

                if post_p >= 0.99:
                    pathogenicity = "Pathogenic"
                elif post_p >= 0.9:
                    pathogenicity = "Likely pathogenic"
                elif post_p <= 0.01:
                    pathogenicity = "Benign"
                elif post_p <= 0.1:
                    pathogenicity = "Likely benign"
                else:
                    pathogenicity = "VUS"

                # 질병 정보.
                if gene_symbol != "-":
                    diseases = disease_db_dic[gene_symbol]

                print(
                    var_id,
                    var_feature,
                    consequence,
                    gene_symbol,
                    ACMG_rule,
                    post_p,
                    pathogenicity,
                    diseases,
                    sep="\t",
                    file=outfile,
                )
