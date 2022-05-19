# Module for PP3, BP4, (BP7) rule
# Predictive data, computational evidence

import pytest, mock

REVEL_DB = "/data/projects/ACMG/database/revel_with_transcript_ids"


@pytest.fixture
def spliceai_db_dic():

    db_dic = {
        "1-69270-A-G": {"OR4F5": [0.01, 0.08, 0.00, 0.00]},
        "1-69511-A-G": {"OR4F5": [0.00, 0.00, 0.02, 0.00]},
        "1-865738-A-G": {
            "SAMD11": [0.00, 0.00, 0.00, 0.00],
            "AL645608.1": [0.6, 0.00, 0.00, 0.00],
        },
        "7-135347239-G-GC": {},
        "1-111436857-TC-GT": {"CD53": None},
    }

    return db_dic


@pytest.fixture
def df_col2idx():

    df_col2idx: dict = {
        "var_id": 0,  # 1-69270-A-G
        "chrom": 1,
        "location": 2,
        "gene": 3,
        "feature": 4,
        "consequence": 5,
        "cDNA_pos": 6,
        "CDS_pos": 7,
        "protein_pos": 8,
        "AA_change": 9,
        "codon_change": 10,
        "symbol": 11,
        "strand": 12,
        "revel": 13,
        "gnomad_ac": 14,
        "gnomad_an": 15,
        "gnomad_af": 16,
    }

    return df_col2idx


@pytest.fixture
def variant_dic():

    variant_dic = {
        "1-35142-G-A": {
            "ENST00000417324": {
                "var_infos": [
                    "1-35142-G-A",
                    "1",
                    [866319],
                    "-",
                    "ENST00000417324",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            "ENST00000342066": {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "-",
                    "ENST00000342066",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
        },
        "1-69614-C-T": {
            "ENST00000335137": {
                "var_infos": [
                    "1-69614-C-T",
                    "1",
                    [179427536],
                    "-",
                    "ENST00000335137",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "TTN-AS1",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            "NM101": {
                "var_infos": [
                    "1-69614-C-T",
                    "1",
                    [179427536],
                    "-",
                    "NM101",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "TTN-AS1",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
        },
    }

    return variant_dic


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid\n"
        "1,35142,35142,G,A,T,M,0.027,ENST00000417324\n"
        "1,35142,35142,G,C,T,R,0.035,ENST00000417324\n"
        "1,35142,35142,G,T,T,K,0.043,ENST00000417324\n"
        "1,35143,35143,T,A,T,S,0.018,ENST00000417324\n"
        "1,35143,35143,T,C,T,A,0.034,ENST00000417324\n"
        "1,35143,35143,T,G,T,P,0.039,ENST00000417324\n"
        "1,69614,69614,C,T,P,L,0.107,ENST00000534990;ENST00000335137\n"
        "1,69616,69616,A,G,R,G,0.066,ENST00000534990;ENST00000335137\n"
        "1,69616,69616,A,T,R,W,0.091,ENST00000534990;ENST00000335137\n"
    ),
)
def test_add_revel_into_var_infos(mock_open, variant_dic, df_col2idx):

    expected = {
        "1-35142-G-A": {
            "ENST00000417324": {
                "var_infos": [
                    "1-35142-G-A",
                    "1",
                    [866319],
                    "-",
                    "ENST00000417324",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "1",
                    0.027,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            "ENST00000342066": {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "-",
                    "ENST00000342066",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
        },
        "1-69614-C-T": {
            "ENST00000335137": {
                "var_infos": [
                    "1-69614-C-T",
                    "1",
                    [179427536],
                    "-",
                    "ENST00000335137",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "TTN-AS1",
                    "1",
                    0.107,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            "NM101": {
                "var_infos": [
                    "1-69614-C-T",
                    "1",
                    [179427536],
                    "-",
                    "NM101",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "TTN-AS1",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
        },
    }

    assert expected == add_revel_into_var_infos(
        "some_path.txt", variant_dic, df_col2idx
    )


@pytest.mark.parametrize(
    "var_infos_dic, expected",
    [
        (
            {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "SAMD11",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            True,
        ),
        (
            {
                "var_infos": [
                    "2-179427536-T-C",
                    "2",
                    [179427536],
                    "ENSG00000237298",
                    "ENST00000591332",
                    "missense_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "TTN-AS1",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            True,
        ),
        (
            {
                "var_infos": [
                    "7-107838464-G-T",
                    "7",
                    [107838464],
                    "ENSG00000091129",
                    "ENST00000351718",
                    "stop_gained,NMD_transcript_variant",
                    "1697/6228",
                    "1269/3552",
                    "423/1183",
                    "V",
                    "gtC/gtA",
                    "NRCAM",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            False,
        ),
        (
            {
                "var_infos": [
                    "12-57114100-A-G",
                    "12",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "missense_variant,splice_region_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "NACA",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            True,
        ),
        (
            {
                "var_infos": [
                    "18-55268866-C-T",
                    "18",
                    [55268866],
                    "ENSG00000134440",
                    "ENST00000589001",
                    "NMD_transcript_variant,missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "NARS",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            True,
        ),
    ],
)
def test_is_missense_var(var_infos_dic, expected, df_col2idx):

    assert expected == is_missense_var(var_infos_dic, df_col2idx)


@pytest.mark.parametrize(
    "var_infos_dic, expected",
    [
        (
            {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "SAMD11",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            False,
        ),
        (
            {
                "var_infos": [
                    "2-179427536-T-C",
                    "2",
                    [179427536],
                    "ENSG00000237298",
                    "ENST00000591332",
                    "splice_region_variant,synonymous_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "TTN-AS1",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            True,
        ),
        (
            {
                "var_infos": [
                    "7-107838464-G-T",
                    "7",
                    [107838464],
                    "ENSG00000091129",
                    "ENST00000351718",
                    "synonymous_variant,NMD_transcript_variant",
                    "1697/6228",
                    "1269/3552",
                    "423/1183",
                    "V",
                    "gtC/gtA",
                    "NRCAM",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            True,
        ),
        (
            {
                "var_infos": [
                    "12-57114100-A-G",
                    "12",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "missense_variant,splice_region_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "NACA",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            False,
        ),
    ],
)
def test_is_synonymous_var(var_infos_dic, expected, df_col2idx):

    assert expected == is_synonymous_var(var_infos_dic, df_col2idx)


@pytest.mark.parametrize(
    "var_infos_dic, expected",
    [
        (
            {
                "var_infos": [
                    "1-69270-A-G",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "OR4F5",
                    "1",
                    0.52,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (1, 0),
        ),
        (
            {
                "var_infos": [
                    "1-69270-A-G",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "no gene",
                    "1",
                    0.12,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 1),
        ),
        (
            {
                "var_infos": [
                    "1-69270-A-G",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "missense_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "OR4F5",
                    "1",
                    0.11,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 1),
        ),
        (
            {
                "var_infos": [
                    "1-865738-A-G",
                    "1",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "missense_variant,splice_region_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "AL645608.1",
                    "-1",
                    0.1,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (1, 0),
        ),
        (
            {
                "var_infos": [
                    "1-865738-A-G",
                    "1",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "missense_variant,splice_region_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "AL645608.1",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (1, 0),
        ),
        (
            {
                "var_infos": [
                    "1-865738-A-G",
                    "1",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "missense_variant,splice_region_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "SAMD11",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 1),
        ),
        (
            {
                "var_infos": [
                    "1-111436857-TC-GT",
                    "1",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "missense_variant,splice_region_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "CD53",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 0),
        ),
    ],
)
def test_predict_missense_pathogenicity(
    var_infos_dic, expected, df_col2idx, spliceai_db_dic
):

    assert expected == predict_missense_pathogenicity(
        var_infos_dic, df_col2idx, spliceai_db_dic
    )


@pytest.mark.parametrize(
    "var_infos_dic, expected",
    [
        (
            {
                "var_infos": [
                    "1-69511-A-G",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "OR4F5",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 1),
        ),
        (
            {
                "var_infos": [
                    "1-865738-A-G",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "SAMD11",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 1),
        ),
        (
            {
                "var_infos": [
                    "1-865738-A-G",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "AL645608.1",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (1, 0),
        ),
        (
            {
                "var_infos": [
                    "1-8657-A-G",
                    "1",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "AL645",
                    "-1",
                    0.1,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 0),
        ),
        (
            {
                "var_infos": [
                    "7-135347239-G-GC",
                    "1",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "intron_variant,non_coding_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "some gene",
                    "-1",
                    0.1,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            },
            (0, 0),
        ),
    ],
)
def test_ppredict_others_pathogenicity(
    var_infos_dic, expected, df_col2idx, spliceai_db_dic
):

    assert expected == predict_others_pathogenicity(
        var_infos_dic, df_col2idx, spliceai_db_dic
    )


#############################################################


def add_revel_into_var_infos(
    file: str, variant_dic: dict, df_col2idx: dict
) -> dict:
    """_summary_
    Note:
        REVEL db의 용량이 매우 크므로, 따로 저장하지 않고, REVEL score값을 variant_infos
        리스트의 마지막에 append한다.

        ##chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
        # 1,35142,35142,G,A,T,M,0.027,ENST00000417324

        먼저, REVEL db를 한 줄씩 읽어들이면서, 각 라인에 해당하는 variant_id를 환자의 variant 목록에
        있는지 확인한 뒤, 환자 var_id의 각 feature가 REVEL DB의 ensembl_transcipt_id에 해당하는지 확인한다.
        이때, vep annotation된 variant의 feature는 RefSeq과 Ensembl transcriptID가 모두 존재하지만,
        REVEL은 Ensembl_transcript id만 존재한다. 따라서, 별도의 RefSeq-Ensembl의 대응과정을 거치지 않고,
        Ensembl에 해당하는 Variant 들만 REVEL score를 부여하였다. 그 외에는, default 값으로 None을 부여함.

    Args:
        variant_dic (dict): class VariantDF.variant_dic
        df_col2idx (dict): DF column index dictionary

    Returns:
        dict: 각각의 variant_info에 revel값이 update된 variant_dic

    Examples:
        >>> ['12-49500509-C-T', 'ENSG00000139636', 'ENST00000551169',
            'missense_variant,NMD_transcript_variant', '72/771', '73/120',
            '25/39', 'V/I', 'Gta/Ata', 'LMBR1L', None]
            ->
            ['12-49500509-C-T', 'ENSG00000139636', 'ENST00000551169',
            'missense_variant,NMD_transcript_variant', '72/771', '73/120',
            '25/39', 'V/I', 'Gta/Ata', 'LMBR1L', 0.051]
    """

    with open(file) as infile:
        for line in infile:
            if line.startswith("chr"):
                f_col2idx = {
                    val: idx for idx, val in enumerate(line.strip().split(","))
                }
                continue
            else:  # new variation
                row = line.strip().split(",")
                chr = row[f_col2idx["chr"]]
                pos = row[f_col2idx["hg19_pos"]]
                ref = row[f_col2idx["ref"]]
                alt = row[f_col2idx["alt"]]
                var_id = f"{chr}-{pos}-{ref}-{alt}"  # 1-985445-G-GT

            if var_id in variant_dic:
                for var_feature in variant_dic[var_id]:
                    # if has same transcipt_id w/ REVEL
                    if var_feature in row[f_col2idx["Ensembl_transcriptid"]]:
                        variant_dic[var_id][var_feature]["var_infos"][
                            df_col2idx["revel"]
                        ] = float(row[f_col2idx["REVEL"]])
                    else:
                        pass  # default == None

    return variant_dic


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


def predict_missense_pathogenicity(
    var_info_dic: dict, df_col2idx: dict, spliceai_db_dic: dict,
) -> tuple:
    """_summary_
    Note:
        REVEL score, SpliceAI score를 확인하여 pp3, bp4룰 부여하는 함수.

        REVEL score를 우선으로 판정한다. REVEL score가 0.5가 넘으면 pp3를 부여하였으며,
        REVEL_score가 없거나 0.5이하인 경우에는 SpliceAI의 값을 사용하여 score가
        0.5를 넘으면 pp3, 그렇지 않으면 bp4를 부여하였다. 두 값이 모두 없는 경우에는 둘 다 0을 부여하였다.

        SpliceAI dictionary의 경우, 같은 ID라도 gene symbol에 따라 값이 다르므로, 이를 구분하였다.
        **SpliceAI
        {
            ("id (1-1571640-G-T)): { <predict_info_dic>

                ( gene_symbol_1: [DS_AG, DS_AL, DS_DG, DS_DL] )
                    "CDK11B": [0.00, 0.01, 0.00, 0.00]
                    "PREAMEF12": [0.12, 0.11, 0.00, 0.01]
            }
        }
    Args:
        var_infos_dic (dict): variation-feature의 세부 정보 딕셔너리.
        df_col2idx (dict): column index dictionary
        spliceai_db_dic (dict): score_dic_of SpliceAI

    Returns:
        tuple(int, int): (pp3, bp4)

    Examples:
        >>> 1. (revel, spliceAI) = (0.52, [0.1, 0.1, 0.0, 0.0])
                -> pp3 = 1, bp4 = 0
            2. (revel, spliceAI) = (0.12, [0.1, 0.1, 0.0, 0.0])
                -> pp3 = 0, bp4 = 1
            3. (revel, spliceAI) = (0.1, [0.5, 0.6, 0.1, 0.0])
                -> pp3 = 1, bp4 = 0
            4. (revel, spliceAI) = (None, [0.1, 0.1, 0.0, 0.0])
                -> pp3 = 0, bp4 = 1
            5. (revel, spliceAI) = (None, [0.6, 0.1, 0.0, 0.0])
                -> pp3 = 1, bp4 = 0
            6. (revel, spliceAI) = (None, [., ., ., .])
                -> pp3 = 0, bp4 = 0
    """

    pp3, bp4 = 0, 0
    spliceai_score = None

    # 1-69270-A-G
    var_id = var_info_dic["var_infos"][df_col2idx["var_id"]]

    # float or None
    revel_score = var_info_dic["var_infos"][df_col2idx["revel"]]
    if revel_score:
        if revel_score >= 0.5:
            pp3 = 1
            return (pp3, bp4)  # (1, 0)
        else:
            bp4 = 1  # temporaily (0, 1)

    # spliceai, revel이 모두 없는 경우도 있음.
    # Spliceai_score 확인
    if spliceai_db_dic.get(var_id):
        gene_symbol = var_info_dic["var_infos"][df_col2idx["symbol"]]
        # gene_symbol에 따라 값이 다름.
        if gene_symbol in spliceai_db_dic[var_id]:
            try:
                spliceai_score = max(spliceai_db_dic[var_id][gene_symbol])
                if spliceai_score >= 0.5:
                    pp3, bp4 = 1, 0
                else:
                    pp3, bp4 = 0, 1
            except TypeError:  # when spliceai score is Nonetype
                pass

    return (pp3, bp4)


def predict_others_pathogenicity(
    var_info_dic: dict, df_col2idx: dict, spliceai_db_dic: dict,
) -> tuple:
    """_summary_
    Note:
        SpliceAI score를 확인하여 pp3, bp4룰 부여하는 함수. missense variant가 아닌 경우에는,
        SpliceAI score만 사용하여 pp3, bp4를 부여하였다. Score가 없는 경우, pp3, bp4에 0을 부여함.

        SpliceAI dictionary의 경우, 같은 ID라도 gene symbol에 따라 값이 다르므로, 이를 구분하였다.
        **SpliceAI
        {
            ("id (1-1571640-G-T)): { <predict_info_dic>

                ( gene_symbol_1: [DS_AG, DS_AL, DS_DG, DS_DL] )
                    "CDK11B": [0.00, 0.01, 0.00, 0.00]
                    "PREAMEF12": [0.12, 0.11, 0.00, 0.01]
            }
        }
    Args:
        var_infos_dic (dict): variation-feature의 세부 정보 딕셔너리.
        df_col2idx (dict): column index dictionary
        spliceai_db_dic (dict): score_dic_of SpliceAI

    Returns:
        tuple(int, int): (pp3, bp4)

    Examples:
        >>> 1. (spliceAI) = ([0.1, 0.1, 0.0, 0.0])
                -> pp3 = 0, bp4 = 1
            2. (spliceAI) = ([0.5, 0.6, 0.1, 0.0])
                -> pp3 = 1, bp4 = 0
            3. (spliceAI) = ([., ., ., .])
                -> pp3 = 0, bp4 = 0
    """
    pp3, bp4 = 0, 0

    """spliceAI"""
    # 1-69270-A-G
    var_id = var_info_dic["var_infos"][df_col2idx["var_id"]]
    # Spliceai_score 확인
    if spliceai_db_dic.get(var_id):
        gene_symbol = var_info_dic["var_infos"][df_col2idx["symbol"]]
        # gene_symbol에 따라 값이 다름.
        if gene_symbol in spliceai_db_dic[var_id]:
            try:
                spliceai_score = max(spliceai_db_dic[var_id][gene_symbol])
                if spliceai_score >= 0.5:
                    pp3 = 1
                else:
                    bp4 = 1
            except TypeError:  # when spliceai result is Nonetype
                pass

    return (pp3, bp4)


def is_synonymous_var(var_infos_dic: dict, df_col2idx: dict) -> bool:

    """_summary_
    Note:
        해당 변이의 결과에 synonymous가 포함되어 있는지 확인하는 함수.

    Args:
        var_infos_dic (dict): variation-feature의 세부 정보 딕셔너리.
        df_col2idx (dict): column index dictionary

    Returns:
        bool: True or False
    """

    var_infos = var_infos_dic["var_infos"]  # uploadted_id

    if "synonymous_variant" in var_infos[df_col2idx["consequence"]]:
        return True
    else:
        return False


def execute(proband_var_df: object, spliceai_db_dic: dict) -> object:
    """_summary_
    Note: ACMG rule 중에서, pp3/bp4/bp7 에 해당하는 룰을 구현한 모듈이다. 각각의 rule에 대한 설명은 다음과 같다.

        pp3: computational evidence에 의하여 유전자나 생산물에 이상이 있을 것으로 예측되는 변이
        bp4: computational evidence에 의하여 유전자나 생산물에 이상이 없을 것으로 예측되는 변이
        bp7: synonymous variants. computational evidence에 의하여 splicing site에 영향을 주지 않을 것으로
            예측되거나, not conserved nucleotide에서의 변이.

        본 모듈에는 2가지의 예측 프로그램이 사용되었다.
        - REVEL: an ensemble method for predicting the pathogenicity of missense variants.
        - SpliceAI: an mothod for predicting splicing from primary sequence with deep learning.

        SpliceAI는 제공되는 tool을 다운받아, 환자의 VCF file을 input으로 하여 얻은 예측값을 main() 함수에서
        딕셔너리 형태로 가공하였다.
        REVEL은 제공되는 data를 온라인으로 다운받아 사용하였다. 다만, 파일용량이 큰 관계로 SpliceAI와는 달리 별도로
        가공하지 않고, 본 모듈의 add_revel_into_var_infos() 함수를 사용하여 환자의 각 variant_info에
        추가한 뒤, 본 모듈 내부에서 SpliceAI와 함께 pp3/bp4/bp7 rule을 평가하는데 사용되었다.

        본 모듈의 작동방식은 다음과 같다.
        먼저, REVEL database의 정보를 환자의 variant_info에 추가한다(add_revel_into_var_infos()).
        이후, 각 variant에 대하여 missense variant인지 아닌지 조사한다(is_missense_var()). 만약 해당 변이의 전사체가
        missense variant라면, REVEL score와 SpliceAI score를 사용하여 pp3와 bp4를 부여한다(predict_missense_
        pathogenicity()). 만약, 해당 변이의 전사체가 missense variant가 아니라면, SpliceAI Score를 사용하여 pp3,
        bp4를 부여하고(predict_others_pathogenicity()), bp4를 부여 받은 변이가 synonymous variant라면(is_
        synonymous_var()), 추가로 bp7을 부여한다. 이 때, bp7의 조건에는 non-conserved 라는 조건도 존재하지만, 본
        모듈에서는 이를 고려하지 않았다.

        missense variant인 경우에는, REVEL과 SpliceAI를 모두 사용하였다. 다만, REVEL을 우선으로하여 REVEL score가
        0.5가 넘으면 pp3를 부여하였으며, REVEL_score가 없거나 0.5이하인 경우에는 SpliceAI의 값을 사용하여 score가
        0.5를 넘으면 pp3, 그렇지 않으면 bp4를 부여하였다. missense variant가 아닌 경우에는, SpliceAI score만
        사용하여 pp3, bp4를 부여하였다.

        데이터베이스 제공자가 제시하는 기준은 다음과 같다.
        <REVEL>
        > 0.5 (likely disease causing), < 0.5 (likely benign)
        <SpliceAI>
        0.2 (high recall), 0.5 (recommended), 0.8 (high precision)

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        spliceai_db_dic (dict): 환자의 각 variant에 대하여 spliceai로 미리 계산된 예측값 dictionary

    Returns:
        object: proband_variant_DF의 variant_dic["evidnece_score_dic"]이 update된 object.

    References:
        1. https://github.com/Illumina/SpliceAI - SpliceAI cutoff
        2. https://asia.ensembl.org/info/genome/variation/prediction/protein_function.html - REVEL
    """

    df_col2idx = proband_var_df.df_col2idx
    variant_dic = proband_var_df.variant_dic
    variant_dic = add_revel_into_var_infos(variant_dic, df_col2idx)

    for var_id in variant_dic:
        for var_feature in variant_dic[var_id]:

            pp3, bp4, bp7 = 0, 0, 0
            if is_missense_var(variant_dic[var_id][var_feature], df_col2idx):
                # refer REVEL SpliceAI
                pp3, bp4 = predict_missense_pathogenicity(
                    variant_dic[var_id][var_feature],
                    df_col2idx,
                    spliceai_db_dic,
                )
            else:  # refer SpliceAI
                pp3, bp4 = predict_others_pathogenicity(
                    variant_dic[var_id][var_feature],
                    df_col2idx,
                    spliceai_db_dic,
                )
                if (bp4 == 1) and is_synonymous_var(
                    variant_dic[var_id][var_feature], df_col2idx
                ):  # if synonymous variant w/ BP4
                    bp7 = 1

            variant_dic[var_id][var_feature]["evidence_score_dic"]["pp3"] = pp3
            variant_dic[var_id][var_feature]["evidence_score_dic"]["bp4"] = bp4
            variant_dic[var_id][var_feature]["evidence_score_dic"]["bp7"] = bp7

    proband_var_df.variant_dic = variant_dic

    return proband_var_df
