from .helper import *
from .rules import *
from collections import defaultdict

# python -m ACMG
# input file
PROBAND_VEP = "/data/projects/ACMG/input/proband.preprocessed_37.txt"
PROBAND_VCF = "/data/projects/ACMG/input/proband.preprocessed.vcf"
FATHER_VCF = "/data/projects/ACMG/input/proband.father.preprocessed.vcf"
MOTHER_VCF = "/data/projects/ACMG/input/proband.mother.preprocessed.vcf"

# database
CLINVAR_DB = "/data/projects/ACMG/database/clinvar_parsed_single.txt"
DISEASE_DB = "/data/projects/ACMG/database/disease.txt"
GNOMAD_DB = "/data/projects/ACMG/database/gnomad.exomes.r2.1.1.sites.vcf.gz"
REVEL_DB = "/data/projects/ACMG/database/revel_with_transcript_ids"
SPLICEAI_DB = "/data/projects/ACMG/database/proband.spliceai.vcf"
REPEAT_DB = "/data/projects/ACMG/database/hg19.fa.out"


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
                        - var_infos: ["1-21929260-A-G", "1", ["1234"], "5909",
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


def main():
    """_summary_
    : 환자와 부모의 VCF 파일을 VEP로 annotation한 파일을 parsing 하여, 각각의 변이 정보가 담긴 data
    frame을 구성한다. 또한, 필요에 따라 다양한 데이터베이스(clinvar, spliceai 등..)를 파징하여 dictionary
    형식으로 저장하거나 VariantDF.variant_dic에 정보를 추가한다. ACMG rule을 구현한 모듈들에서 각 변이를
    분석하여 rule들을 할당한 뒤, update 된 VariantDF object 를 반환한다.

    최종적으로, 각각의 변이에 할당된 evidence rule 수를 바탕으로 bayesian framework를 이용하여 pathogenicity
    를 판별한 뒤, 이를 내림차순으로 출력한다.
    """

    # VEP annotated VCF 파일을 저장하는 과정
    proband_var_df = VariantDF()
    proband_var_df.parse_variant_file(PROBAND_VEP)

    # de novo 판정을 위한 genotype 정보 저장하는 과정
    proband_genotype_dic = parse_vcf_genotype(PROBAND_VCF)
    father_genotype_dic = parse_vcf_genotype(FATHER_VCF, "F")
    mother_genotype_dic = parse_vcf_genotype(MOTHER_VCF, "M")

    # 데이터베이스 전처리 과정
    spliceai_db_dic = dbparser.parse_spliceai_db(SPLICEAI_DB)
    clinvar_db_dic, clinvar_col2idx = dbparser.parse_clinvar_db(CLINVAR_DB)
    disease_db_dic, disease_col2idx = dbparser.parse_disease_db(DISEASE_DB)
    repeat_db_dic = dbparser.parse_repeatmasker_db(REPEAT_DB)

    # ACMG module
    # revel-5분. pp3 4210, bp4 728559, bp7 70829
    proband_var_df = pp3bp4bp7.execute(proband_var_df, spliceai_db_dic)
    # pp2 2335, bp1 7766
    proband_var_df = pp2bp1.execute(
        proband_var_df, clinvar_db_dic, clinvar_col2idx
    )
    # pvs1 176.
    proband_var_df = pvs1.execute(
        proband_var_df,
        clinvar_db_dic,
        clinvar_col2idx,
        disease_db_dic,
        disease_col2idx,
    )
    # gnomad-1시간. pm2 1894, ba1 525127, bs1 12380
    proband_var_df = pm2ba1bs1.execute(
        proband_var_df, disease_db_dic, disease_col2idx
    )
    # 약 14초. ps1 = 77, pm5 = 310(ps1 과 중복상태), pp5 =  88, bp6 =  170850
    proband_var_df = ps1pm5pp5bp6.execute(
        proband_var_df, clinvar_db_dic, clinvar_col2idx
    )
    # 8초. ps2 33492(vcf = 2851개) de novo
    proband_var_df = ps2.execute(
        proband_var_df,
        proband_genotype_dic,
        father_genotype_dic,
        mother_genotype_dic,
    )
    # 약 20초. pm4 = 1261, bp3 = 766
    proband_var_df = pm4bp3.execute(proband_var_df, repeat_db_dic)
    bayesframe.calculate_acmg(proband_var_df, disease_db_dic)


if __name__ == "__main__":
    # 약 3500초 소요
    main()
