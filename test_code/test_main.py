from collections import defaultdict
import pytest, mock


@pytest.fixture
def variant_df():
    variant_df = VariantDF()
    return variant_df


@pytest.fixture
def df_col2idx(variant_df):
    return variant_df.df_col2idx


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        "## SYMBOL_SOURCE : Source of gene symbol\n"
        "## HGNC_ID : Stable identifer of HGNC gene symbol\n"
        "## REFSEQ_MATCH : RefSeq transcript match status\n"
        "## SOURCE : Source of transcript\n"
        "## REFSEQ_OFFSET : HGVS adjustment length required due to mismatch between RefSeq transcript and the reference genome\n"
        "## GIVEN_REF : Reference allele from input\n"
        "## USED_REF : Reference allele as used to get consequences\n"
        "## BAM_EDIT : Indicates success or failure of edit using BAM file\n"
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
        "1-866319-G-A\t1:866319\tA\tENSG00000187634\tENST00000341065\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=G;IMPACT=MODIFIER;STRAND=1;FLAGS=cds_start_NF;SYMBOL=SAMD11;SYMBOL_SOURCE=HGNC;HGNC_ID=28706;SOURCE=Ensembl;GIVEN_REF=G;USED_REF=G\n"
        "1-866319-G-A\t1:866319\tA\tENSG00000187634\tENST00000342066\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=G;IMPACT=MODIFIER;STRAND=1;SYMBOL=SAMD11;SYMBOL_SOURCE=HGNC;HGNC_ID=28706;SOURCE=Ensembl;GIVEN_REF=G;USED_REF=G\n"
        "1-866319-G-A\t1:866319\tA\tENSG00000268179\tENST00000598827\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=G;IMPACT=MODIFIER;STRAND=-1;SYMBOL=AL645608.1;SYMBOL_SOURCE=Clone_based_ensembl_gene;SOURCE=Ensembl;GIVEN_REF=G;USED_REF=G\n"
        "2-179427536-T-C\t2:179427536\tC\tENSG00000237298\tENST00000591332\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=T;IMPACT=MODIFIER;STRAND=1;SYMBOL=TTN-AS1;SYMBOL_SOURCE=HGNC;HGNC_ID=44124;SOURCE=Ensembl;GIVEN_REF=T;USED_REF=T\n"
        "5-127522153-G-GT\t5:127522153-127522154\tT\t6558\tNM_001256461.2\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=-;IMPACT=MODIFIER;STRAND=1;SYMBOL=SLC12A2;SYMBOL_SOURCE=EntrezGene;HGNC_ID=10911;SOURCE=RefSeq;GIVEN_REF=-;USED_REF=-\n"
        "7-107838464-G-T	7:107838464	T	ENSG00000091129	ENST00000351718	Transcript	synonymous_variant	1697/6228	1269/3552	423/1183	V	gtC/gtA	-	REF_ALLELE=G;IMPACT=LOW;STRAND=-1;SYMBOL=NRCAM;SYMBOL_SOURCE=HGNC;HGNC_ID=7994;SOURCE=Ensembl;GIVEN_REF=G;USED_REF=G"
        "9-33799263-C-T\t9:33799263\tT\t5646\tNM_007343.4\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=C;IMPACT=MODIFIER;DISTANCE=34;STRAND=1;SYMBOL=PRSS3;SYMBOL_SOURCE=EntrezGene;HGNC_ID=9486;SOURCE=RefSeq;GIVEN_REF=C;USED_REF=C\n"
        "12-57114100-A-G\t12:57114100\tG\tENSG00000196531\tENST00000552540\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=A;IMPACT=MODIFIER;STRAND=-1;SYMBOL=NACA;SYMBOL_SOURCE=HGNC;HGNC_ID=7629;SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A\n"
        "16-325771-GT-G\t16:325772\t-\tENSG00000076344\tENST00000168869\tTranscript\tintron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=T;IMPACT=MODIFIER;STRAND=-1;SYMBOL=RGS11;SYMBOL_SOURCE=HGNC;HGNC_ID=9993;SOURCE=Ensembl;GIVEN_REF=T;USED_REF=T\n"
        "17-9846521-G-C\t17:9846521\tC\t8522\tNM_201433.2\tTranscript\tsynonymous_variant\t864/8269\t648/1431\t216/476\tT\tacC/acG\t-\tREF_ALLELE=G;IMPACT=LOW;STRAND=-1;SYMBOL=GAS7;SYMBOL_SOURCE=EntrezGene;HGNC_ID=4169;SOURCE=RefSeq;GIVEN_REF=G;USED_REF=G\n"
        "18-55268866-C-T\t18:55268866\tT\tENSG00000134440\tENST00000589001\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tREF_ALLELE=C;IMPACT=MODIFIER;DISTANCE=4300;STRAND=-1;SYMBOL=NARS;SYMBOL_SOURCE=HGNC;HGNC_ID=7643;SOURCE=Ensembl;GIVEN_REF=C;USED_REF=C\n"
        "X-44202890-G-GC\tX:44202890-44202891\tC\tENSG00000183690\tENST00000420999\tTranscript\t5_prime_UTR_variant\t28-29/2464\t-\t-\t-\t-\t-\tREF_ALLELE=-;IMPACT=MODIFIER;STRAND=-1;SYMBOL=EFHC2;SYMBOL_SOURCE=HGNC;HGNC_ID=26233;SOURCE=Ensembl;GIVEN_REF=-;USED_REF=-\n"
    ),
)
def test_parse_variant_file(mock_open, variant_df):
    variant_df.parse_variant_file("some_path.txt")
    expected = {
        "1-866319-G-A": {
            "ENST00000341065": {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "intron_variant",
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
            "ENST00000342066": {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000342066",
                    "intron_variant",
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
            "ENST00000598827": {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "ENSG00000268179",
                    "ENST00000598827",
                    "intron_variant",
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
        },
        "2-179427536-T-C": {
            "ENST00000591332": {
                "var_infos": [
                    "2-179427536-T-C",
                    "2",
                    [179427536],
                    "ENSG00000237298",
                    "ENST00000591332",
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
        "5-127522153-G-GT": {
            "NM_001256461.2": {
                "var_infos": [
                    "5-127522153-G-GT",
                    "5",
                    [127522153, 127522154],
                    "6558",
                    "NM_001256461.2",
                    "intron_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "SLC12A2",
                    "1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            }
        },
        "7-107838464-G-T": {
            "ENST00000351718": {
                "var_infos": [
                    "7-107838464-G-T",
                    "7",
                    [107838464],
                    "ENSG00000091129",
                    "ENST00000351718",
                    "synonymous_variant",
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
            }
        },
        "12-57114100-A-G": {
            "ENST00000552540": {
                "var_infos": [
                    "12-57114100-A-G",
                    "12",
                    [57114100],
                    "ENSG00000196531",
                    "ENST00000552540",
                    "intron_variant",
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
            }
        },
        "16-325771-GT-G": {
            "ENST00000168869": {
                "var_infos": [
                    "16-325771-GT-G",
                    "16",
                    [325772],
                    "ENSG00000076344",
                    "ENST00000168869",
                    "intron_variant,NMD_transcript_variant",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                    "RGS11",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            }
        },
        "17-9846521-G-C": {
            "NM_201433.2": {
                "var_infos": [
                    "17-9846521-G-C",
                    "17",
                    [9846521],
                    "8522",
                    "NM_201433.2",
                    "synonymous_variant",
                    "864/8269",
                    "648/1431",
                    "216/476",
                    "T",
                    "acC/acG",
                    "GAS7",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            }
        },
        "18-55268866-C-T": {
            "ENST00000589001": {
                "var_infos": [
                    "18-55268866-C-T",
                    "18",
                    [55268866],
                    "ENSG00000134440",
                    "ENST00000589001",
                    "downstream_gene_variant",
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
            }
        },
        "X-44202890-G-GC": {
            "ENST00000420999": {
                "var_infos": [
                    "X-44202890-G-GC",
                    "X",
                    [44202890, 44202891],
                    "ENSG00000183690",
                    "ENST00000420999",
                    "5_prime_UTR_variant",
                    "28-29/2464",
                    "-",
                    "-",
                    "-",
                    "-",
                    "EFHC2",
                    "-1",
                    None,
                    0,
                    0,
                    None,
                ],
                "evidence_score_dic": {},
            }
        },
    }
    assert expected == variant_df.variant_dic


@pytest.mark.parametrize(
    "row, expected",
    [
        (
            [
                "1-866319-G-A",
                "1:866319",
                "A",
                "ENSG00000187634",
                "ENST00000341065",
                "Transcript",
                "intron_variant",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "REF_ALLELE=G;IMPACT=MODIFIER;STRAND=1;FLAGS=cds_start_NF;SYMBOL=SAMD11;SYMBOL_SOURCE=HGNC;HGNC_ID=28706;SOURCE=Ensembl;GIVEN_REF=G;USED_REF=G",
            ],
            {
                "var_infos": [
                    "1-866319-G-A",
                    "1",
                    [866319],
                    "ENSG00000187634",
                    "ENST00000341065",
                    "intron_variant",
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
        ),
        (
            [
                "7-107838464-G-T",
                "7:107838464",
                "T",
                "ENSG00000091129",
                "ENST00000351718",
                "Transcript",
                "synonymous_variant",
                "1697/6228",
                "1269/3552",
                "423/1183",
                "V",
                "gtC/gtA",
                "-",
                "REF_ALLELE=G;IMPACT=LOW;STRAND=-1;SYMBOL=NRCAM;SYMBOL_SOURCE=HGNC;HGNC_ID=7994;SOURCE=Ensembl;GIVEN_REF=G;USED_REF=G",
            ],
            {
                "var_infos": [
                    "7-107838464-G-T",
                    "7",
                    [107838464],
                    "ENSG00000091129",
                    "ENST00000351718",
                    "synonymous_variant",
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
        ),
    ],
)
def test_make_variant_info_dic(row, expected, variant_df):

    line = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    f_col2idx = {
        val: idx for idx, val in enumerate(line.strip("#").strip().split("\t"))
    }

    assert expected == variant_df.make_variant_info_dic(row, f_col2idx)


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        "##contig=<ID=X,length=155270560>\n"
        "##contig=<ID=Y,length=59373566>\n"
        "##contig=<ID=MT,length=16569>\n"
        "##reference=file:///data/etc/reference/GRCh37_mainchr/GRCh37.fa\n"
        "##source=SelectVariants\n"
        '##INFO=<ID=OLD_VARIANT,Number=.,Type=String,Description="Original chr:pos:ref:alt encoding">\n'
        '##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">\n'
        '##INFO=<ID=OLD_CLUMPED,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">\n'
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ETD21-RWNYM\n"
        "1	866319	1-866319-G-A	G	A	1225.77	PASS	AC=2;AF=1;AN=2;DP=34;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=34.24;SOR=2.303;set=variant	GT:AD:DP:GQ:PL	1/1:0,34:34:99:1254,102,0\n"
        "1	866511	1-866511-C-CCCCT	C	CCCCT	1633.14	PASS	AC=2;AF=1;AN=2;BaseQRankSum=-0.654;ClippingRankSum=0;DP=85;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60.25;MQRankSum=0.395;QD=30.63;ReadPosRankSum=-0.375;SOR=5.283;set=variant2	GT:AD:DP:GQ:PL	1/1:2,41:43:10:1670,10,0\n"
        "1	871334	1-871334-G-T	G	T	932.77	PASS	AC=1;AF=0.5;AN=2;BaseQRankSum=-4.452;ClippingRankSum=0;DP=74;ExcessHet=3.0103;FS=7.394;MLEAC=1;MLEAF=0.5;MQ=60;MQRankSum=0;QD=12.78;ReadPosRankSum=-0.255;SOR=0.211;set=variant	GT:AD:DP:GQ:PL	0/1:35,38:73:99:961,0,1122\n"
    ),
)
def test_parse_vcf_genotype(mock_open):

    expected = {
        "1-866319-G-A": "1/1",
        "1-866511-C-CCCCT": "1/1",
        "1-871334-G-T": "0/1",
    }

    assert expected == parse_vcf_genotype("some_path.txt", "M")


class VariantDF:
    """_summary_
    Note:
        VEP로 annocation된 vcf 파일의 변환값을 읽어서 "ID-transcipt-infos"로 구성된
    dictionary와 infos[list] 값에 해당하는 index dictionary를 저장하는 dataframe class.
    {id : {feature : { [var_infos], {evidence_score} } } }로 이루어진 3중 dictionary와
    var_infos_list의 index를 지정하기 위한 col2idx dictionary가 있다.

    Examples:
    >>> df_col2idx = {"var_id": 0, "gene": 1, ...}

        variant_dic = {
            "var_key_id": {
                "var_key_feature": {
                    - var_infos: [var_id, gene, ...]
                    - evidence_score_dic: {}}
                }
            }
            "1-69270-A-G": {
                "ENST00000335137": {
                    - var_infos: [var_id, gene, feature, consequence, cDNA_pos, CDS_pos,
                            Protein_pos, AA_change, codon_change, strand, symbol, ...]
                    - evidence_score_dic: dict(bool)
                }
                .
                .
            }
        }"""

    def __init__(self):
        self.df_col2idx: dict = {
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
        self.variant_dic: dict = None

    def parse_variant_file(self, vep_file: str):
        """_summary_
        Note:
            VEP로 annocation된 vcf 파일의 변환값을 읽은 뒤, 각 라인을 파징하여 variant들의 정보
            (ID-transcipt(feature)-infos) 가 담긴 딕셔너리를 VariantDF class의
            변수(self.variant_dic)에 저장한다. 저장된 정보(list)에 대한 index 값은 VariantDF
            class 의 내부 변수 (self.df_col2idx)에 내재되어 있다.

        Args:
            vep_file (str): VEP_file address

        Examples:
        >>> variant_dic = {
                "1-21929260-A-G": {
                    "NR_170901.1": {
                        - var_infos: ["1-21929260-A-G", "1", [1234], "5909",
                            "NR_170901.1","intron_variant,non_coding_transcript_variant",
                                    "-", "-", "-", "-", "-", "-1", "RAP1GAP"]
                        - evidence_score_dic: dict(bool)
                    }
                }
            }"""

        with open(vep_file) as infile:
            # {id : {feature : { [var_infos], {evidence_score} }}}
            variant_dic = defaultdict(dict)
            for line in infile:
                if line.startswith("#") and (not line.startswith("##")):
                    f_col2idx = {
                        val: idx
                        for idx, val in enumerate(
                            line.strip("#").strip().split("\t")
                        )
                    }
                if not line.startswith("#"):  # variation info
                    row = line.strip().split("\t")
                    # "1-69270-A-G", "ENST00000335137"
                    var_key_id = row[f_col2idx["Uploaded_variation"]]
                    var_key_feature = row[f_col2idx["Feature"]]
                    variant_info_dic = self.make_variant_info_dic(
                        row, f_col2idx
                    )
                    variant_dic[var_key_id].update(
                        {var_key_feature: variant_info_dic}
                    )

            self.variant_dic = variant_dic

    def make_variant_info_dic(self, var_row: list, file_col2idx: dict) -> dict:
        """_summary_
        Note:
            ID-feature로 특정지어지는 variant에 대해 {"var_infos": [], "evidence_score_dict": {}}
            로 구성된 <variant_info_dic> 을 만들어 반환한다. "symbol" 의 경우, vep 옵션에서
            추가한 것으로, "Extra" column에 해당하여 따로 추출하였음. Revel의 경우, 추후 revel
            score를 저장하기 위한 용도.

        Args:
            var_row (list): splitted each file line
            file_col2idx (dict): file index dictionary

        Returns:
            dict: {var_infos: [infos], "evidence_score_dict": {}}

        Examples:
        >>> {
                - "var_infos": ["1-21929260-A-G", "1", ["1234"], "5909", "NR_170901.1",
                            "intron_variant,non_coding_transcript_variant",
                            "-", "-", "-", "-", "-", "-1", "RAP1GAP"]
                - "evidence_score_dic": dict(bool)
            }"""

        # variable initialization
        variant_info_dic = {"var_infos": None, "evidence_score_dic": dict()}
        gene_symbol = ""
        strand = ""
        revel = None
        gnomad_ac = 0
        gnomad_an = 0
        gnomad_af = None

        # variable assignment
        chrom, location = var_row[file_col2idx["Location"]].split(":")
        location = list(map(int, location.split("-")))  # [1234] or [1234, 1236]

        if "SYMBOL=" in var_row[file_col2idx["Extra"]]:
            gene_symbol = (
                var_row[file_col2idx["Extra"]].split("SYMBOL=")[1].split(";")[0]
            )
        if "STRAND=" in var_row[file_col2idx["Extra"]]:
            strand = (
                var_row[file_col2idx["Extra"]].split("STRAND=")[1].split(";")[0]
            )

        var_infos = [
            var_row[file_col2idx["Uploaded_variation"]],
            chrom,
            location,
            var_row[file_col2idx["Gene"]],
            var_row[file_col2idx["Feature"]],
            var_row[file_col2idx["Consequence"]],
            var_row[file_col2idx["cDNA_position"]],
            var_row[file_col2idx["CDS_position"]],
            var_row[file_col2idx["Protein_position"]],
            var_row[file_col2idx["Amino_acids"]],
            var_row[file_col2idx["Codons"]],
            gene_symbol,
            strand,
            revel,
            gnomad_ac,
            gnomad_an,
            gnomad_af,
        ]
        variant_info_dic["var_infos"] = var_infos

        return variant_info_dic


def parse_vcf_genotype(vcf_file: str, family="") -> dict:
    """_summary_
    Note:
        variant의 de novo 판정을 위해서, vcf file에 있는 GT 정보를 저장한다.

    Returns:
        dict: {var_id : GT}

    Examples:
    >>> {
            - 1-866319-G-A: '1/1'
            - 1-897325-G-C: '0/1'
            - 1-866511-C-CCCCT: '0/1'
        }
    """
    ## CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ETD21-RWNY
    with open(vcf_file) as infile:

        vcf_genotype_dic = dict()
        # {var_id : GT}
        for line in infile:
            if line.startswith("#") and (not line.startswith("##")):
                f_col2idx = {
                    val: idx
                    for idx, val in enumerate(
                        line.strip("#").strip().split("\t")
                    )
                }
            if not line.startswith("#"):  # variation info
                row = line.strip().split("\t")
                var_id = row[f_col2idx["ID"]]  # 1-985445-G-GT
                var_info = row[
                    f_col2idx[f"ETD21-RWNY{family}"]
                ]  # 1/1:0,7:7:21:226,21,0
                genotypes = var_info.split(":")[0]  # GT:AD:DP:GQ:PL, 1/1
                vcf_genotype_dic[var_id] = genotypes

    return vcf_genotype_dic
