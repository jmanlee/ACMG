# Module for PM2, BA1, BS1
# Population data, vs normal people variants

import pytest, mock


GNOMAD_DB = "/data/projects/ACMG/database/gnomad.exomes.r2.1.1.sites.vcf.gz"


@pytest.fixture
def disease_db_dic():

    disease_db_dic = {
        "INO80": [  # -
            [
                "Combined immunodeficiency with granulomatosis",
                ["-"],
                ["-"],
                ["-"],
                ["-"],
            ]
        ],
        "RAG1": [  # AR
            [
                "Combined immunodeficiency with granulomatosis",
                ["Autosomal recessive", "Autosomal dominant"],
                ["-"],
                ["-"],
                ["-"],
            ],
            [
                "Omenn syndrome",
                ["Autosomal recessive"],
                ["Pediatric"],
                ["HP:0002240", "HP:0001880",],
                [
                    "Hepatomegaly",
                    "Eosinophilia",
                    "Autosomal recessive inheritance",
                ],
            ],
        ],
        "RAG2": [  # AD
            [
                "Combined immunodeficiency with granulomatosis",
                ["Autosomal dominant"],
                ["-"],
                ["-"],
                ["-"],
            ]
        ],
        "VASN": [  # YD
            [
                "some disease name",
                ["Y-linked", "Autosomal recessive"],
                ["-"],
                ["-"],
                ["-"],
            ]
        ],
        "ALX3": [  # XR
            [
                "Frontonasal dysplasia 1",
                ["X-linked recessive", "Sporadic"],
                ["Neonatal"],
                ["HP:0006992", "HP:0001156", "HP:0030084",],
                ["Anterior basal encephalocele", "Brachydactyly",],
            ],
            [
                "Combined immunodeficiency with granulomatosis",
                ["Autosomal dominant"],
                ["-"],
                ["-"],
                ["-"],
            ],
        ],
    }

    return disease_db_dic


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
        "1-138980-G-C": {
            "ENST00000417324": {
                "var_infos": [
                    "1-138980-G-C",
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
                    "1-138980-G-C",
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
            }
        },
        "1-970680-C-T": {
            "ENST00000335137": {
                "var_infos": [
                    "1-970680-C-T",
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
            }
        },
    }

    return variant_dic


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        '##INFO=<ID=controls_nhomalt_popmax,Number=A,Type=Integer,Description="Count of homozygous individuals in the population with the maximum allele frequency in the controls subset">\n'
        "##contig=<ID=Y,length=59373566,assembly=gnomAD_GRCh37>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1	138980	rs796175332	G	C	11050.36	PASS	AC=2;AN=134554;AF=1.48639e-05;rf_tp_probability=2.14504e-01;\n"
        "1	139253	rs1193670589	G	C	1651.52	PASS	AC=1;AN=133560;AF=7.48727e-06;rf_tp_probability=5.21094e-01;FS=1.9411\n"
        "1	955647	rs1194105483	C	T	709.35	PASS	AC=1;AN=83836;AF=1.19281e-05;rf_tp_probability=4.21673e-01;FS=3.60300e+00;InbreedingCo\n"
        "1	970679	rs777369292	T	C	8511.30	PASS	AC=3;AN=250896;AF=1.19571e-05;rf_tp_probability=7.91696e-01;FS=\n"
        "1	970680	rs746480380	C	T	4074.49	PASS	AC=0;AN=0;rf_tp_probability=7.78117e-01;FS=1.27600e+00;Inb\n"
    ),
)
def test_add_gnomad_into_var_infos(mock_open, variant_dic, df_col2idx):

    expected = {
        "1-138980-G-C": {
            "ENST00000417324": {
                "var_infos": [
                    "1-138980-G-C",
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
                    2,
                    134554,
                    1.48639e-05,
                ],
                "evidence_score_dic": {},
            },
            "ENST00000342066": {
                "var_infos": [
                    "1-138980-G-C",
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
                    2,
                    134554,
                    1.48639e-05,
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
            }
        },
        "1-970680-C-T": {
            "ENST00000335137": {
                "var_infos": [
                    "1-970680-C-T",
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
            }
        },
    }

    assert expected == add_gnomad_into_var_infos(
        "some_path.txt", variant_dic, df_col2idx
    )


@pytest.mark.parametrize(
    "an, af, expected",
    [(100, 0.6, 0), (100, 0.03, 0), (1001, 0.0014, 0), (10023, 0.06, 1)],
)
def test_assign_ba1_rule(an: int, af: float, expected):

    assert expected == assign_ba1_rule(an, af)


@pytest.mark.parametrize(
    "an, af, inheritence, expected",
    [
        (100, 0.6, "AD", (0, 0)),
        (10001, 0.00001, "XD", (1, 0)),
        (10001, 0.0001, "XD", (0, 0)),
        (10001, 0.0005, "AD", (0, 1)),
        (1002221, 0.00001, "AR", (1, 0)),
        (10023, 0.0002, "XR", (0, 0)),
        (55212, 0.01, "XR", (0, 1)),
        (122222, 0.00001, "YD", (1, 0)),
        (122222, 0.00001, "Y", (0, 0)),
    ],
)
def test_compare_with_disease_inheritence(an, af, inheritence, expected):

    assert expected == compare_with_disease_inheritence(an, af, inheritence)


@pytest.mark.parametrize(
    "an, af, gene_symbol, expected",
    [
        (100, 0.6, "RAG2", (0, 0)),  # AD
        (10001, 0.00001, "RAG2", (1, 0)),
        (10001, 0.0001, "RAG2", (0, 0)),
        (10001, 0.0005, "RAG2", (0, 1)),
        (1002221, 0.00001, "RAG1", (1, 0)),  # AR
        (10023, 0.00023, "ALX3", (0, 0)),  # XR
        (55212, 0.01, "ALX3", (0, 1)),
        (122222, 0.00001, "VASN", (1, 0)),  # YD
        (122222, 0.00001, "INO80", (0, 0)),  # -
    ],
)
def test_assign_pm2bs1_rule(an, af, gene_symbol, disease_db_dic, expected):

    disease_col2idx: dict = {
        "title": 0,
        "inheritance": 1,
        "onsetAges": 2,
        "symtoms_id": 3,
        "symtoms": 4,
    }

    assert expected == assign_pm2bs1_rule(
        an, af, gene_symbol, disease_db_dic, disease_col2idx
    )


#######################


def add_gnomad_into_var_infos(
    file: str, variant_dic: dict, df_col2idx: dict
) -> dict:
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
    with open(file) as infile:
        for line in infile:
            if line.startswith("#") and (not line.startswith("##")):
                f_col2idx = {
                    val: idx
                    for idx, val in enumerate(
                        line.strip("#").strip().split("\t")
                    )
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
                    gnomad_ac = (
                        row[f_col2idx["INFO"]].split("AC=")[1].split(";")[0]
                    )
                    gnomad_an = (
                        row[f_col2idx["INFO"]].split("AN=")[1].split(";")[0]
                    )
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

    if inheritence in ["AD", "XD", "YD"]:  # dominant disease
        if gnomad_af < 0.00002:  # 0.002%
            pm2 = 1
        elif gnomad_af > 0.0002:  # 0.02%
            bs1 = 1
    elif inheritence in ["AR", "XR"]:  # recessive disease
        if gnomad_af < 0.0001:  # 0.01%
            pm2 = 1
        elif gnomad_af > 0.001:  # 0.1%
            bs1 = 1
    else:
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
        Y_count = 0

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

            Y_count += each_disease_infos[disease_col2idx["inheritance"]].count(
                "Y-linked"
            )

        # 보고된 비율이 더 많은 쪽을 inheritence로 설정.
        # 같은 수인 경우 dominant로 설정하며, X-linked 우선.
        if (XD_count + XR_count) > 0:
            if XD_count >= XR_count:
                inheritence = "XD"
            else:
                inheritence = "XR"
        elif Y_count > 0:
            inheritence = "YD"
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
