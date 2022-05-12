# Module for PS1, PM5, PP5, BP6
# reference data of ClinVar

DNA_codon_pair = """TTT F      CTT L      ATT I      GTT V
TTC F      CTC L      ATC I      GTC V
TTA L      CTA L      ATA I      GTA V
TTG L      CTG L      ATG M      GTG V
TCT S      CCT P      ACT T      GCT A
TCC S      CCC P      ACC T      GCC A
TCA S      CCA P      ACA T      GCA A
TCG S      CCG P      ACG T      GCG A
TAT Y      CAT H      AAT N      GAT D
TAC Y      CAC H      AAC N      GAC D
TAA *      CAA Q      AAA K      GAA E
TAG *      CAG Q      AAG K      GAG E
TGT C      CGT R      AGT S      GGT G
TGC C      CGC R      AGC S      GGC G
TGA *      CGA R      AGA R      GGA G
TGG W      CGG R      AGG R      GGG G """.split()

DNA_CODON_TABLE = dict(zip(DNA_codon_pair[0::2], DNA_codon_pair[1::2]))


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


def check_amino_acid_change_in_clinvar(
    proband_var_infos: list,
    df_col2idx: dict,
    clinvar_gene_var_dic: dict,
    clinvar_col2idx: dict,
):
    """_summary_
    Note:
        해당 missense 변이가 ClinVar에 보고된 아미노산 변이와 동일한 변화를 일으키는지 화긴하는 함수.

        - proband_var_infos: [var_id, gene, feature, consequence, cDNA_pos, CDS_pos,
                            Protein_pos, AA_change, codon_change, strand, symbol, etc...]

        - clinvar_gene_var_dic ("gene_symbol"):
        {
            "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}, G245R::G289R]
            "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}, K570N]
        }

        VEP annotation 파일에는 각 variant 의 아미노산 및 codon의 전후 정보가 담겨져 있다. 따라서,
        해당 변이가 코돈에서 몇 번째에 위치하고 있는 지 알수 있으며, transcript 마다의 protein 길이가
        다르더라도, 레퍼런스 지놈상의 위치를 이용하여 위치를 특정지을 수 있다. 이 때, strand 의 방향에 따라
        코돈을 읽는 방향이 달라지므로 이에 주의해야 한다.

        먼저, 해당 변이의 염색체, 위치, ref, alt 정보 및 codon, AA 변화도 저장한다. 이 때, strand가 -1
        이면, 해당 코돈을 뒤집은 뒤 상보적인 서열로 바꿔준다. 코돈 상에서 바뀐 서열은 대문자로 표시되어 있으므로,
        코돈에서 몇번째 염기가 바뀌어 있는지 알 수 있고, 이를 바탕으로 해당 코돈의 시작점을 찾고 각 위치에 해당하는
        염기들을 dictionary 로 만든다. aAAGgt -> {1000: 'A', 1001: 'A', .. ,1005: 'T'}

        다음으로, ClinVar에서 해당 변이와 같은 위치, 혹은 비슷한 위치에 있는 변이를 검색한다. 이 때,
        보고된 변이가 SNV가 아니라 여러개의 서열이 바뀌어 있을 수 있으므로, 해당 코돈 구간보다 최대 3만큼 넓은
        범위를 탐색한다. 탐색이 되었다면, 마찬가지로 해당변이도 인덱스-염기 dictionary를 만든다.

        마지막으로, 변이가 포함된 코돈 영역에 해당하는 clinvar 변이가 있다면, 염기를 치환한 뒤 아미노산으로 변경한다.
        만약, pathogenic으로 보고된 경우, 환자의 바뀐 아미노산와 같다면 ps1을, 다르다면 pm5를 부여한다.

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        df_col2idx (dict): DF column index dictionary
        clinvar_db_dic (dict): ClinVar 데이터베이스를 미리 가공한 dictionary
        clinvar_col2idx (dict): Clinvar_db_dic-index dictionary

    Returns:
        tuple: (ps1, pm5)

    Examples:
        >>> ( proband_var_id, pb_ref_codon, pb_alt_codon, clinvar_var_id, strand,
            cv_altered_codon, pb_ref_aa, pb_alt_aa, cv_alt_aa ) -> if pathogenic -> ps1 or pm5

            ( 1-31349647-C-T, taC, taT, 1-31349647-C-T, -1, ATA, V, I, I) -> ps1
            ( 3-133475812-G-A, Ggc, Agc, 3-133475813-G-A, 1, GAC, G, S, D) -> pm5
            ( 17-39742898-GC-AT, GCa ATa, 17-39742898-G-A, -1, TGT, C, Y, C) -> 부여 x


    """

    pathogenic_set = {"Pathogenic", "Likely pathogenic"}
    ps1, pm5 = 0, 0

    # 13-32950906-AAG-GAC
    proband_var_id = proband_var_infos[df_col2idx["var_id"]]
    strand = proband_var_infos[df_col2idx["strand"]]

    # variant_ID 관련 정보 할당. "13", "32950906", "AAG", "GAC"
    pb_chrom, pb_chr_pos, pb_ref_nuc, pb_alt_nuc = proband_var_id.split("-")
    pb_chr_pos = int(pb_chr_pos)

    # codons 변화 정보. e.g. aAAGgt/aGACgt -> "aAAGgt", "aGACgt"
    pb_codon_changes = proband_var_infos[df_col2idx["codon_change"]].split("/")
    pb_ref_codon, pb_alt_codon = (pb_codon_changes[0], pb_codon_changes[1])
    codon_len = len(pb_ref_codon)  # 6

    # 반대 Strand를 읽는 경우, 반대 방향으로 염기 수정.
    if strand == "-1":
        pb_ref_codon = pb_ref_codon[::-1].translate(
            str.maketrans("ATCGatcg", "TAGCtagc")
        )
        pb_alt_codon = pb_alt_codon[::-1].translate(
            str.maketrans("ATCGatcg", "TAGCtagc")
        )

    # Amino acids 변화 정보. e.g. K/E, DE/AK -> ["K", "E"], ["DE", "AK"]
    pb_aa_changes = proband_var_infos[df_col2idx["AA_change"]].split("/")
    pb_ref_aa, pb_alt_aa = pb_aa_changes[0], pb_aa_changes[1]

    # 바뀐 코돈 영역 (3의배수)의 시작지점, 끝지점 찾기. "aAAGgt".find("AAG") = 1
    alt_start_index = pb_alt_codon.find(pb_alt_nuc)
    first_codon_start = pb_chr_pos - alt_start_index
    last_codon_end = first_codon_start + (codon_len - 1)

    # print(pb_alt_codon, pb_alt_nuc, alt_start_index)

    # e.g. aAAGgt -> {1000: 'A', 1001: 'A', .. ,1005: 'T'}
    ref_nuc_idx_dic = dict()
    for i in range(codon_len):
        ref_nuc_idx_dic[first_codon_start + i] = pb_ref_codon[i].upper()

    # Clinvar_dic 에서 환자 근처에 있는 데이터 찾기. 해당 gene_symbol에 대한 정보만 있음.
    for clinvar_var_id, var_infos in clinvar_gene_var_dic.items():

        if clinvar_var_id == "-":
            continue
        # e.g. clinvar_id information 13-32950906-CA-AT
        cv_chrom, cv_chr_pos, cv_ref_nuc, cv_alt_nuc = clinvar_var_id.split("-")
        cv_chr_pos = int(cv_chr_pos)  # 변이 시작 위치
        cv_var_len = len(cv_ref_nuc)  # missense, indel 변이의 길이

        # missense가 아닌경우 제외
        if len(cv_ref_nuc) != len(cv_alt_nuc):
            continue

        # 변이의 위치정보 저장. e.g. {1000: 'A', 1001: 'A', .. ,1005: 'T'}
        cv_nuc_idx_dic = dict()
        for i in range(cv_var_len):
            cv_nuc_idx_dic[cv_chr_pos + i] = cv_alt_nuc[i].upper()

        # variant가 포함된 코돈 구간을 기준으로, +-3 영역까지 포함된 clinvar data 선정
        # 그후, 환자의 변이가 포함된 ref 코돈을 clinvar data로 치환한다.
        cv_altered_codons = list()
        if pb_chrom == cv_chrom:
            # 변이가 포함된 코돈 영역에 대한 Clinvar data check 후 ref codon 치환
            if first_codon_start - 3 <= cv_chr_pos <= last_codon_end + 3:
                for i in range(codon_len):
                    position = first_codon_start + i
                    if cv_nuc_idx_dic.get(position):
                        cv_altered_codons.append(cv_nuc_idx_dic[position])
                    else:
                        cv_altered_codons.append(ref_nuc_idx_dic[position])
            else:
                continue
        else:
            continue
        # e.g. "GAC", "GACAAT"
        cv_altered_codon: str = "".join(cv_altered_codons)

        # 반대 Strand를 읽는 경우, 반대 방향으로 codon 수정.
        if strand == "-1":
            cv_altered_codon = cv_altered_codon[::-1].translate(
                str.maketrans("ATCG", "TAGC")
            )

        # codon을 amino acid로 변환.
        cv_altered_AAs = list()
        for i in range(0, codon_len, 3):
            codon = cv_altered_codon[i : i + 3]
            cv_altered_AAs.append(DNA_CODON_TABLE[codon])

        # e.g. "R", "DA"
        cv_alt_aa: str = "".join(cv_altered_AAs)

        # clinvar 변이로 치환한 뒤, reference AA와 같은 경우
        if pb_ref_aa == cv_alt_aa:
            pass
        else:
            cv_var_patho = var_infos[clinvar_col2idx["pathogenicity"]]
            # amino acid is changed, and if pathogenic
            if cv_var_patho in pathogenic_set:
                if pb_alt_aa == cv_alt_aa:
                    ps1 = 1
                else:
                    pm5 = 1

        """
        print(
            proband_var_id,
            pb_ref_codon,
            pb_alt_codon,
            clinvar_var_id,
            cv_altered_codon,
            pb_ref_aa,
            pb_alt_aa,
            var_infos[clinvar_col2idx["aa_change"]],
            cv_alt_aa,
        )
        """

    return (ps1, pm5)


def check_same_variant_in_clinvar(
    proband_variant_id: str,
    gene_symbol: str,
    clinvar_db_dic: dict,
    clinvar_col2idx: dict,
) -> tuple:
    """_summary_
    Note:
        아래와 같이 gene_symbol을 기준으로 정리된 ClinVar dictionary를 바탕으로, 해당 유전자에 생긴 변이가
        기존에 보고된 바가 있는지 확인하고, pathogenic/benign 여부에 따라 pp5, bp6 를 부여한다.

        - gene_variant_dic:
        {
            "gene_symbol":{
                "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}, G245R::G289R]
                "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}, K570N]
            }
        }

    Args:
        gene_variant_dic (dict): gene_symbol을 기준으로 정리한 ClinVar data dictionary
        gene_symbol (str): HGVS gene symbol
        clinvar_db_dic (dict): ClinVar 데이터베이스를 미리 가공한 dictionary
        clinvar_col2idx (dict): Clinvar_db_dic-index dictionary

    Returns:
        tuple: (pp5, bp6)
    """

    pathogenic_set = {"Pathogenic", "Likely pathogenic"}
    benign_set = {"Benign", "Likely benign"}

    pp5, bp6 = 0, 0

    if clinvar_db_dic.get(gene_symbol):
        if clinvar_db_dic[gene_symbol].get(proband_variant_id):
            var_patho = clinvar_db_dic[gene_symbol][proband_variant_id][
                clinvar_col2idx["pathogenicity"]
            ]
            if var_patho in pathogenic_set:
                pp5 = 1
            elif var_patho in benign_set:
                bp6 = 1

    return (pp5, bp6)


def execute(
    proband_var_df: object, clinvar_db_dic: dict, clinvar_col2idx: dict
) -> object:
    """_summary_
    Note: ACMG rule 중에서, ps1/pm5/pp5/bp6 에 해당하는 룰을 구현한 모듈이다. 각각의 rule에 대한 설명은 다음과 같다.

        ps1: 기존에 알려진 pathogenic variant 와 같은 아미노산 변화를 유발하는 변이 (같거나/다른 변이 모두 해당).
        pm5: 기존에 pathogenic 하다고 보고된 missense mutation과 다른 아미노산 변화를 일으키는 새로운 missense 변이
        pp5: 신뢰할 만한 reference 에서 pathogenic 으로 보고된 변이
        bp6: 신뢰할 만한 reference 에서 benign 으로 보고된 변이

        본 모듈에서는 1가지의 데이터베이스가 사용되었다.
        - ClinVar: 동일한 유전자를 기준으로 하여 변이들을 정리하였다(RCV).
        {
            "BRCA2": {
                "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}, G245R::G289R]
                "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}, K570N]
            }
        }

        본 모듈의 작동방식은 다음과 같다.
        먼저, 각 variant에 대하여 missense variant인지 아닌지 조사한다(is_missense_var()).
        만약 해당 변이의 전사체가 missense variant라면, 해당 변이의 gene_symbol과 일치하는 ClinVar 데이터들을
        조사하여, clinvar에 보고된 아미노산 변이와 같은 변이를 유발하는지 확인한다(check_amino_acid_change_in_clinvar()).
        만약, 기존에 pathogenic or likely pathogenic 으로 보고된 변이와 같은 아미노산 변이를 유발하면
        ps1을 부여하고, 다른 아미노산 변이를 유발하면 pm5를 부여한다.
        만약, 해당변이의 전사체가 missense variant가 아니라면, 해당 변이가 clinvar에 보고된 변이인지 조사한다.
        (check_same_variant_in_clinvar()). 만약, 해당 변이가 clinvar에 보고되어 있고 pathogenic 하다면
        pp5를, benign 하다면 bp6를 부여한다.

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

            ps1, pm5, pp5, bp6 = 0, 0, 0, 0
            gene_symbol = variant_dic[var_id][var_feature]["var_infos"][
                df_col2idx["symbol"]
            ]
            if is_missense_var(variant_dic[var_id][var_feature], df_col2idx):
                if (gene_symbol != "-") and (clinvar_db_dic.get(gene_symbol)):
                    # gene_symbol이 없는 경우, 비교의미 없음.
                    ps1, pm5 = check_amino_acid_change_in_clinvar(
                        variant_dic[var_id][var_feature]["var_infos"],
                        df_col2idx,
                        clinvar_db_dic[gene_symbol],
                        clinvar_col2idx,
                    )
            else:  # 동일한 염기 서열 변화가 있는지 조사. gene_symbol 과 무관. 염기의 유무만 판단
                pp5, bp6 = check_same_variant_in_clinvar(
                    var_id, gene_symbol, clinvar_db_dic, clinvar_col2idx
                )

            variant_dic[var_id][var_feature]["evidence_score_dic"]["ps1"] = ps1
            variant_dic[var_id][var_feature]["evidence_score_dic"]["pm5"] = pm5
            variant_dic[var_id][var_feature]["evidence_score_dic"]["pp5"] = pp5
            variant_dic[var_id][var_feature]["evidence_score_dic"]["bp6"] = bp6

    proband_var_df.variant_dic = variant_dic

    return proband_var_df
