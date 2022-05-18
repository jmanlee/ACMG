# module for dasebase parsing

from collections import defaultdict
import gzip


def read_big_file(filename: str) -> str:
    """_summary_
    Note:
        input file의 각 line 정보를 generator 형식으로, 매 호출시 반환한다.

    Args:
        filename (str): file address

    Yields:
        line (str): 각 line의 정보를 매 호출 시 string으로 반환
    """

    with open(filename) as infile:
        for line in infile:
            if not line:
                break
            yield line


def read_big_gz_file(filename: str) -> str:
    """_summary_
    Note:
        input gz 압축 file의 각 line 정보를 generator 형식으로, 매 호출시 반환한다.

    Args:
        filename (str): file address

    Yields:
        line (str): 각 line의 정보를 매 호출 시 string으로 반환
    """

    with gzip.open(filename, "rb") as infile:
        for line in infile:
            if not line:
                break
            yield line.decode(encoding="utf-8")


def parse_spliceai_db(filename: str) -> dict:
    """_summary_
    Note:
        환자의 vcf 파일을 input으로 하여 얻은 spliceai 예측값을 parsing 하여, 각 변이 id에 대해
        spliceai prediction score를 저장한다. 이 때, 동일한 변이 위치라도 1 가지 이상의 유전자가
        해당될 수 있기 때문에, "gene_symbol"을 key 값 으로 갖는 sub_dictionary를 생성한다.

    Args:
        filename (str): file address

    Returns:
        spliceai_db_dic (dict): { "var_id":
                                    { "gene1": [predict score], "gene2": [PS] }
                                }
    Examples:
        >>> {
                ("id (1-1571640-G-T)): { <predict_info_dic>

                    ( gene_symbol_1: [DS_AG, DS_AL, DS_DG, DS_DL] )
                        "CDK11B": [0.00, 0.01, 0.00, 0.00]
                        "PREAMEF12": [0.12, 0.11, 0.00, 0.01]
                }
            }
    """

    with open(filename) as infile:

        spliceai_db_dic = dict()
        # {id : { gene_symbol : [predict score]}}
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
                predict_info_dic = make_spliceai_var_info_dic(row, f_col2idx)
                spliceai_db_dic[var_id] = predict_info_dic

    return spliceai_db_dic


def make_spliceai_var_info_dic(var_row: list, file_col2idx: dict) -> dict:
    """_summary_
    Note:
        spliceai 결과값이 추가된 vcf 파일의 각 라인을 파징하여, {"gene1": [score], "gene2": [score]}
        로 구성된 <predict_info_dic> 을 만들어 반환한다. 이 때, 각 "id"에 대하여 여러 gene이 해당하는 경우가
        있다. 또한, 값이 없는 경우(.|.|.|.)가 있기 때문에, ValueError 예외처리를 통해 None값을 저장하였다.

    Args:
        var_row (list): splitted each file line
        file_col2idx (dict): file index dictionary

    Returns:
        dict: { "gene1": [predict score], "gene2": [DS_AG, DS_AL, DS_DG, DS_DL] }

    Examples:
        >>> {
                "CDK11B": [0.00, 0.01, 0.00, 0.00]
                "PREAMEF12": [0.12, 0.11, 0.00, 0.01]
            }
    """
    ## CHROM POS    ID           REF    ALT         QUAL    FILTER
    #  1   985445 1-985445-G-GT    G      GT       1488.73    PASS
    ## INFO
    # AC=1;AF=0.5;AN=2;BaseQRankSum=0.063;ClippingRankSum=0;DP=133;ExcessHet=3.0103;
    # FS=2.331;MLEAC=1;MLEAF=0.5;MQ=60;MQRankSum=0;QD=11.91;ReadPosRankSum=0.012;SOR=0.926;set=variant2;
    # (Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">)
    # SpliceAI=T|MMP23B|0.00|0.00|0.00|0.00|-17|-27|-11|-15,T|CDK11B|0.00|0.01|0.00|0.00|-9|-37|-46|10

    predict_info_dic = dict()  # initialization

    if "SpliceAI=" in var_row[file_col2idx["INFO"]]:

        # store with as gene
        predict_values: list = (
            var_row[file_col2idx["INFO"]]
            .split("SpliceAI=")[1]
            .split(";")[0]
            .split(",")
        )
        # multiple genes
        for i in range(len(predict_values)):
            predict_result: list = predict_values[i].split("|")
            gene_symbol = predict_result[1]

            try:
                values = list(map(float, predict_result[2:6]))
                predict_info_dic[gene_symbol] = values
            except ValueError:
                predict_info_dic[gene_symbol] = None

    return predict_info_dic


def parse_revel_db(filename: str) -> tuple():  # ok
    """사용하지 않는 함수. dbparser 대신, iteration으로 각 variant의 info란에 저장하는 형식으로 변경.
    _summary_
    Note:
        가능한 missense variation에 대한 예측 score를 제공하는, revel database(7Gb)를 파징하여 dictionary
        형식으로 반환하는 함수. {"var_id" : [revel_infos]}. transcript_id가 여러 개인 경우 있음.

    Args:
        filename (str): file address

    Returns:
        dict:   {
                    "var_key": [aaref, aaalt, revel_score, [enst_ids]]
                }
    Examples:
        >>> {
                "1-69728-T-C": ["I", "T", 0.035, ["ENST00000534990", "ENST00000335137"]],
            }
    """

    revel_col2idx: dict = {
        "aaref": 0,
        "aaalt": 1,
        "revel_score": 2,
        "enst_id": 3,
    }

    with open(filename) as infile:

        revel_db_dic = dict()
        # {id : [proedict score]}
        for line in infile:
            if line.startswith("chr"):
                f_col2idx = {
                    val: idx for idx, val in enumerate(line.strip().split(","))
                }
            else:  # new info
                row = line.strip().split(",")
                chr = row[f_col2idx["chr"]]
                pos = row[f_col2idx["hg19_pos"]]
                ref = row[f_col2idx["ref"]]
                alt = row[f_col2idx["alt"]]
                var_id = f"{chr}-{pos}-{ref}-{alt}"  # 1-985445-G-GT

                predict_infos: list = make_revel_predict_infos(row, f_col2idx)
                revel_db_dic[var_id] = predict_infos

    return revel_db_dic, revel_col2idx


def make_revel_predict_infos(row_var: list, file_col2idx: dict) -> list:
    """_summary_
    Note:
        revel의 각 라인을 파징하여, 예측값에 대한 정보들을 리스트로 반환한다. 이 때, Ensembl transcript id는
        여러개가 있을 수 있으므로, split 하여 리스트로 저장한다.

    Args:
        var_row (list): splitted each file line
        file_col2idx (dict): file index dictionary

    Returns:
        list: [aaref, aaalt, revel_score, [enst_ids]]

    Examples:
    >>> ["I", "T", 0.035, ["ENST00000534990", "ENST00000335137"]
    """
    ##chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
    # 1,35142,35142,G,A,T,M,0.027,ENST00000417324
    # 1,35142,35142,G,C,T,R,0.035,ENST00000417324

    ensembl_ids: list = row_var[file_col2idx["Ensembl_transcriptid"]].split(";")

    predict_infos: list = [
        row_var[file_col2idx["aaref"]],
        row_var[file_col2idx["aaalt"]],
        float(row_var[file_col2idx["REVEL"]]),
        ensembl_ids,
    ]

    return predict_infos


def parse_clinvar_db(filename: str) -> tuple:
    """_summary_
    Note:
        ClinVar database 파일을 parsing 하여, "gene_symbol"을 key로 하는 dictionary를
        반환한다. 파일의 라인은 각각 특정한 variant_id(e.g. 13-32950906-C-A)를 가지고 있으며,
        변이에 대한 복수의 보고가 정리되어 있다. 이 때, 동일한 variant_id가 존재하는 경우가 간혹 있는데,
        (직접 확인했을 때) 내용 중 대부분 증상에 관한 차이만 있으므로 별도로 저장하지 않고 덮어쓴다. 동일한 유전자에
        해당하는 var_id는 { "var_id": [molecular_consequences, pathogenicity(dict)] } 와 같은
        sub_dictionary 형태로 기존의 value(dict)에 update해서 저장한다.

        * molecular_consquence: Pathogenic, Benign 등..
        * pathogenicity: missense variant, intron variant 등..

    Args:
        filename (str): file address

    Returns:
        dict: { "gene symbol": { "var_id": [pathogenicity, molecular_consequences(dict), aa_change] }}

    Examples:
        >>> {
                "BRCA2": {
                    "13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}, "G245R::G289R"]
                    "13-32950906-C-A": ["Pathogenic", {"missense variant": 2, "intron variant": 1}, "K570N"]
                }
            }
    """

    clinvar_col2idx: dict = {
        "pathogenicity": 0,
        "consequence_dic": 1,
        "aa_change": 2,
    }

    with open(filename) as infile:

        clinvar_db_dic = defaultdict(dict)
        # { symbol: {var_id: [var_infos]}}
        for line in infile:
            if line.startswith("#"):
                f_col2idx = {
                    val: idx
                    for idx, val in enumerate(
                        line.strip("#").strip().split("\t")
                    )
                }
            else:  # clinvar info
                row = line.strip().split("\t")
                gene_symbol = row[f_col2idx["variant:gene:symbol"]]
                gene_symbols = gene_symbol.split("::")
                clinvar_var_info_dic = make_clinvar_gene_info_dic(
                    row, f_col2idx
                )
                for g_symbol in gene_symbols:
                    clinvar_db_dic[g_symbol].update(clinvar_var_info_dic)

    return clinvar_db_dic, clinvar_col2idx


def make_clinvar_gene_info_dic(var_row: list, file_col2idx: dict) -> dict:
    """_summary_
    Note:
        ClinVar database의 각 라인을 파징하여, { symbol: {var_id: [var_infos]}}와 같이 구성된
        <variant_info_dic>을 만들어 반환한다. 이 때, pathogenicity의 경우, 한 id의 rcv에 다수의 보고가 포함되므로
        각 보고의 pathogencitity의 수를 count하여 dictionary 형태로 저장한다.

    Args:
        var_row (list): splitted each file line
        file_col2idx (dict): file index dictionary

    Returns:
        dict: { "var_id": [pathogenicity, molecular_consequences(dict), aa_change] }

    Examples:
        >>> {"13-32921028-CTTTCGG-C": ["Pathogenic", {"splice donor variant": 1}, "G245R::G289R"] }
    """

    clinvar_var_info_dic = dict()

    var_id: str = var_row[file_col2idx["normalized_variant:GRCh37:CPRA"]]
    var_pathogenicity: str = var_row[file_col2idx["pathogenicity"]]
    var_consequences: list = var_row[
        file_col2idx["variant:molecular_consequence"]
    ].split("::")
    var_aa_change: str = var_row[file_col2idx["variant:Protein1LetterCode"]]

    var_consequence_dic = defaultdict(int)

    for cons in var_consequences:
        var_consequence_dic[cons.split(":")[-1]] += 1

    clinvar_var_info_dic[var_id] = [
        var_pathogenicity,
        var_consequence_dic,
        var_aa_change,
    ]

    return clinvar_var_info_dic


def parse_disease_db(filename: str) -> tuple:
    """_summary_
    Note:
        자체적으로 구축한 disease database 파일을 parsing 하여, "gene_symbol"을 key로 하는
        dictionary를 반환한다. 파일의 라인은 각각 특정한 질병에 대해, gene symbol, onset age,
        inheritence 정보 등이 정리되어 있다. 이 때, 동일한 gene symbol이 존재하는 경우가 간혹 있는데,
        동일한 유전자에 해당하는 질병 정보는 [title, inheritence, onsetAge] 와 같은 리스트의 형태로
        value(list에 append해서 저장한다. 최종적으로 생선된 dictionary와 저장된 column의 index
        정보를 반환한다.

    Args:
        filename (str): file address

    Returns:
        dict: { "geneSymbol": [ ["title", "inheritance", "onsetAges", "symtoms_id", "symtoms"] ] }
        disease_col2idx: dict = {"title": 0, "inheritance": 1, "onsetAges": 2, "symtoms_id": 3, "symtoms": 4}

    Examples:
        >>> {
                "TBCE": [
                    ["Ciliary dyskinesia, primary, 14" ,["Autosomal recessive"], ["Pediatric"], symtoms_id, symtoms],
                    ["Encephalopathy, progressive, with amyotrophy and optic atrophy" ,["Autosomal recessive"]
                    , ["Infancy", "Neonatal"], symtoms_id, symtoms]
                ]
            }
            disease_col2idx: dict = {"title": 0, "inheritance": 1, "onsetAges": 2, "symtoms_id": 3, "symtoms": 4}
    """

    disease_col2idx: dict = {
        "title": 0,
        "inheritance": 1,
        "onsetAges": 2,
        "symtoms_id": 3,
        "symtoms": 4,
    }

    with open(filename) as infile:

        disease_db_dic = defaultdict(list)
        # { symbol: [ [disease_infos] ] }
        for line in infile:
            if line.startswith("#"):
                f_col2idx = {
                    val: idx
                    for idx, val in enumerate(
                        line.strip("#").strip().split("\t")
                    )
                }
            else:  # disease info
                row = line.strip().split("\t")
                gene_symbol = row[f_col2idx["geneSymbol"]]
                gene_symbols = gene_symbol.split("::")
                disease_infos: list = make_disease_info_list(row, f_col2idx)

                for g_symbol in gene_symbols:
                    disease_db_dic[g_symbol].append(disease_infos)

    return disease_db_dic, disease_col2idx


def make_disease_info_list(disease_row: list, file_col2idx: dict) -> list:
    """_summary_
    Note:
        Disease database의 각 라인을 파징하여, ["title", "inheritance", "onsetAges"] 와
        같이 구성된 <disease_info_list>를 만들어 반환한다.

    Args:
        disease_row (list): splitted each file line
        file_col2idx (dict): file index dictionary

    Returns:
        list: ["title", "inheritances", "onsetAges", "symtoms_id", "symtoms"]

    Examples:
        >>> ["Ciliary dyskinesia, primary, 14" ,["Autosomal recessive"], ["Pediatric"], symtoms_id, symtoms]
    """

    disease_infos = list()

    # contents
    disease_title: str = disease_row[file_col2idx["title"]]
    disease_inheritances: list = disease_row[
        file_col2idx["inheritances:value"]
    ].split("||")
    disease_onset_age: list = disease_row[
        file_col2idx["onsetAges:value"]
    ].split("||")
    disease_symtoms_id: list = disease_row[file_col2idx["symptoms:id"]].split(
        "||"
    )
    disease_symtoms: list = disease_row[file_col2idx["symptoms:name"]].split(
        "||"
    )

    # append infos
    disease_infos.extend(
        [
            disease_title,
            disease_inheritances,
            disease_onset_age,
            disease_symtoms_id,
            disease_symtoms,
        ]
    )

    return disease_infos


def parse_repeatmasker_db(filename: str) -> dict:
    """_summary_
    Note:
        Repeat region에 대한 database 파일을 parsing 하여, "Chromosome"을 key로 하는 dictionary를
        반환한다. 파일의 라인은 각각 특정한 염색체의 반복서열 범위의 start, end를 포함하고 있다. 동일한 염색체에 해
        당하는 (start, end) 세트를 리스트에 append한 뒤, 염색체를 key로 하는 dictionary에 저장한다.

    Args:
        filename (str): file address

    Returns:
        dict: {"chrom#": [(start1, end1), (start2, end2)]}

    Examples:
        >>> {
                "1": [(10001, 10467), (10469, 11447), ...]
                "2": [(444214, 444451), (445112, 445234), ...]
            } (start, end)는 오름차순 정렬되어 있음.
    """

    with open(filename) as infile:

        repeatmasker_db_dic = defaultdict(list)
        # { chr# : [(start, end), ..]}

        # 첫 세줄이 column 명 관련. 따로 파징하기 어려워서 3줄 readline.
        a = infile.readline()
        b = infile.readline()
        c = infile.readline()

        # append position infos
        for line in infile:
            row = line.strip().split()
            chrom = row[4]  # chr1
            start_pos = int(row[5])
            end_pos = int(row[6])

            repeatmasker_db_dic[chrom[3:]].append((start_pos, end_pos))

    return repeatmasker_db_dic
