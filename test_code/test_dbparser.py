# module for dasebase parsing test

from collections import defaultdict
import gzip
import pytest, mock


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        '##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">\n'
        '##INFO=<ID=OLD_CLUMPED,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">\n'
        '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">\n'
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ETD21-RWNY\n"
        "1	69270	1-69270-A-G	A	G	197.8	SNP_FILTER	AC=2;AF=1;AN=2;DP=7;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=30.06;QD=28.26;SOR=4.174;set=FilteredInAll;SpliceAI=G|OR4F5|0.01|0.08|0.00|0.00|-10|26|-28|-25	GT:AD:DP:GQ:PL	1/1:0,7:7:21:226,21,0\n"
        "1	69511	1-69511-A-G	A	G	2934.77	PASS	AC=2;AF=1;AN=2;DP=97;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=43.29;QD=30.26;SOR=0.992;set=variant;SpliceAI=G|OR4F5|0.00|0.00|0.02|0.00|26|32|26|34	GT:AD:DP:GQ:PL	1/1:0,97:97:99:2963,291,0\n"
        "1	865738	1-865738-A-G	A	G	2970.77	PASS	AC=1;AF=0.5;AN=2;BaseQRankSum=0.777;ClippingRankSum=0;DP=237;ExcessHet=3.0103;FS=3.944;MLEAC=1;MLEAF=0.5;MQ=60;MQRankSum=0;QD=12.59;ReadPosRankSum=-0.25;SOR=0.971;set=variant;SpliceAI=G|SAMD11|0.00|0.00|0.00|0.00|-8|2|-22|-50,G|AL645608.1|0.00|0.00|0.00|0.00|-7|-19|-17|-41	GT:AD:DP:GQ:PL	0/1:131,105:236:99:2999,0,3791\n"
        "7	135347239	7-135347239-G-GC	G	GC	2201.73	PASS	AC=1;AF=0.5;AN=2;DP=89;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.5;MQ=60;QD=26.21;SOR=1.139;set=variant2;OLD_MULTIALLELIC=7:135347239:GC/G/GCC;OLD_CLUMPED=7:135347239:GC/GCC	GT:AD:DP:GQ:PL	./1:0,53:84:99:2239,924,975\n"
        "1	111436857	1-111436857-TC-GT	TC	GT	1495.77	PASS	AC=2;AF=1;AN=2;DP=34;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=23.93;SOR=2.584;set=variant;SpliceAI=GT|CD53|.|.|.|.|.|.|.|.	GT:AD:DP:GQ:PGT:PID:PL	1/1:0,34:34:99:1|1:111436857_T_G:1524,102,0\n"
    ),
)
def test_parse_spliceai_db(mock_open):

    expected = {
        "1-69270-A-G": {"OR4F5": [0.01, 0.08, 0.00, 0.00]},
        "1-69511-A-G": {"OR4F5": [0.00, 0.00, 0.02, 0.00]},
        "1-865738-A-G": {
            "SAMD11": [0.00, 0.00, 0.00, 0.00],
            "AL645608.1": [0.00, 0.00, 0.00, 0.00],
        },
        "7-135347239-G-GC": {},
        "1-111436857-TC-GT": {"CD53": None},
    }

    assert expected == parse_spliceai_db("spliceai.vcf")


@pytest.mark.parametrize(
    "row, expected",
    [
        (
            [
                "1",
                "69270",
                "1-69270-A-G",
                "A",
                "G",
                "197.8",
                "SNP_FILTER",
                "AC=2;AF=1;AN=2;DP=7;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=30.06;QD=28.26;SOR=4.174;set=FilteredInAll;SpliceAI=G|OR4F5|0.01|0.08|0.00|0.00|-10|26|-28|-25",
                "GT:AD:DP:GQ:PL",
                "1/1:0,7:7:21:226,21,0",
            ],
            {"OR4F5": [0.01, 0.08, 0.00, 0.00]},
        ),
        (
            [
                "1",
                "69511",
                "1-69511-A-G",
                "A",
                "G",
                "2934.77",
                "PASS",
                "AC=2;AF=1;AN=2;DP=97;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=43.29;QD=30.26;SOR=0.992;set=variant;SpliceAI=G|OR4F5|0.00|0.00|0.02|0.00|26|32|26|34",
                "GT:AD:DP:GQ:PL",
                "1/1:0,97:97:99:2963,291,0",
            ],
            {"OR4F5": [0.00, 0.00, 0.02, 0.00]},
        ),
        (
            [
                "1",
                "865738",
                "1-865738-A-G",
                "A",
                "G",
                "2970.77",
                "PASS",
                "AC=1;AF=0.5;AN=2;BaseQRankSum=0.777;ClippingRankSum=0;DP=237;ExcessHet=3.0103;FS=3.944;MLEAC=1;MLEAF=0.5;MQ=60;MQRankSum=0;QD=12.59;ReadPosRankSum=-0.25;SOR=0.971;set=variant;SpliceAI=G|SAMD11|0.00|0.00|0.00|0.00|-8|2|-22|-50,G|AL645608.1|0.00|0.00|0.00|0.00|-7|-19|-17|-41",
                "GT:AD:DP:GQ:PL",
                "0/1:131,105:236:99:2999,0,3791",
            ],
            {
                "SAMD11": [0.00, 0.00, 0.00, 0.00],
                "AL645608.1": [0.00, 0.00, 0.00, 0.00],
            },
        ),
        (
            [
                "7",
                "135347239",
                "7-135347239-G-GC",
                "G",
                "GC",
                "2201.73",
                "PASS",
                "AC=1;AF=0.5;AN=2;DP=89;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.5;MQ=60;QD=26.21;SOR=1.139;set=variant2;OLD_MULTIALLELIC=7:135347239:GC/G/GCC;OLD_CLUMPED=7:135347239:GC/GCC",
                "GT:AD:DP:GQ:PL",
                "./1:0,53:84:99:2239,924,975",
            ],
            {},
        ),
        (
            [
                "1",
                "111436857",
                "1-111436857-TC-GT",
                "TC",
                "GT",
                "1495.77",
                "PASS",
                "AC=2;AF=1;AN=2;DP=34;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=23.93;SOR=2.584;set=variant;SpliceAI=GT|CD53|.|.|.|.|.|.|.|.",
                "GT:AD:DP:GQ:PGT:PID:PL",
                "1/1:0,34:34:99:1|1:111436857_T_G:1524,102,0",
            ],
            {"CD53": None},
        ),
    ],
)
def test_make_spliceai_var_info_dic(row, expected):

    line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tETD21-RWNY"
    f_col2idx = {
        val: idx for idx, val in enumerate(line.strip("#").strip().split("\t"))
    }

    assert expected == make_spliceai_var_info_dic(row, f_col2idx)


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        "#clinvar_set_id\ttitle\trcv\tcreation_date\tlast_updated\tlast_evaluated\treview_status\tstar\tpathogenicity\tinheritance\torigin\tindividuals\tPMIDs\tgenotype_type\tgenotype_vcv\tvariant:variation_type\tvariant:vcv\tvariant:variant_type\tvariant:allele_id\tvariant:molecular_consequence\tvariant:HGVSc\tvariant:HGVSp\tvariant:HGVSg\tvariant:HGVS_previous\tvariant:Protein1LetterCode\tvariant:GRCh38:Assembly\tvariant:GRCh38:Chr\tvariant:GRCh38:positionVCF\tvariant:GRCh38:start\tvariant:GRCh38:stop\tvariant:GRCh38:referenceAllele\tvariant:GRCh38:alternateAllele\tvariant:GRCh37:Assembly\tvariant:GRCh37:Chr\tvariant:GRCh37:positionVCF\tvariant:GRCh37:start\tvariant:GRCh37:stop\tvariant:GRCh37:referenceAllele\tvariant:GRCh37:alternateAllele\tvariant:NCBI36:Assembly\tvariant:NCBI36:Chr\tvariant:NCBI36:positionVCF\tvariant:NCBI36:start\tvariant:NCBI36:stop\tvariant:NCBI36:referenceAllele\tvariant:NCBI36:alternateAllele\tvariant:gene:id\tvariant:gene:symbol\tvariant:uniprot\tvariant:dbsnp\tvariant:omim_gene_id\tdisease_id\tdisease_name\tdisease_mechanism\tscv:scv\tscv:pathogenicity\tscv:review_status\tscv:pmid\tscv:submitter\tscv:last_evaluated\tscv:submit_date\tscv:gene_symbol\tscv_pathogenicity\trcv_pathogenicity\tnormalized_variant:GRCh38:Assembly\tnormalized_variant:GRCh38:Chr\tnormalized_variant:GRCh38:positionVCF\tnormalized_variant:GRCh38:referenceAllele\tnormalized_variant:GRCh38:alternateAllele\tnormalized_variant:GRCh37:Assembly\tnormalized_variant:GRCh37:Chr\tnormalized_variant:GRCh37:positionVCF\tnormalized_variant:GRCh37:referenceAllele\tnormalized_variant:GRCh37:alternateAllele\tnormalized_variant:GRCh38:CPRA\tnormalized_variant:GRCh37:CPRA\n"
        "74627773	NM_000769.4(CYP2C19):c.-806C>A AND clopidogrel response - Dosage, Efficacy, Toxicity/ADR	RCV000211201	2016-05-18	2021-09-29	2016-06-14	reviewed by expert panel	3	drug response	-	germline	-	22955794::23364775::20801498::22704413::23922007::21392617::21716271::20826260::23726091::22190063::20083681::22228204::22990067::24019397::19463375::22028352::23809542::22462746::22992668::20492469	-	-	Variant	VCV000225946	single nucleotide variant	227770	-	-	-	NC_000010.11:g.94761900C>A::NG_008384.3:g.4220C>A::NG_055436.1:g.1260C>A	NC_000010.10:g.96521657C>A	-	GRCh38	10	94761900	94761900	94761900	C	A	GRCh37	10	96521657	96521657	96521657	C	A	-	-	-	-	-	-	-	1557::110599570	CYP2C19::LOC110599570	-	rs12248560	-	MedGen:CN236507	clopidogrel response - Dosage, Efficacy, Toxicity/ADR	-	SCV000268179	drug response	reviewed by expert panel	22955794::23364775::20801498::22704413::23922007::21392617::20826260::23726091::22190063::20083681::22228204::22990067::24019397::19463375::22028352::23809542::22462746::22992668::20492469	PharmGKB	2016-06-14	2018-06-18	CYP2C19	-	-	GRCh38	10	94761900	C	A	GRCh37	10	96521657	C	A	10-94761900-C-A	10-96521657-C-A\n"
        "85343921	NM_000274.4(OAT):c.1311G>C (p.Leu437Phe) AND Ornithine aminotransferase deficiency	RCV000000181	2012-08-13	2022-02-19	1992-02-15	no assertion criteria provided	0	Benign	-	germline	-	1737786	-	-	Variant	VCV000000158	single nucleotide variant	15197	NM_001322968.2:c.1311G>C:nonsense::NM_001322974.2:c.711G>C:frameshift variant::NM_001322966.2:c.1311G>C:missense variant::NM_001322965.2:c.1311G>C:missense variant::NM_001171814.2:c.897G>C:missense variant::NM_001322969.2:c.1311G>C:missense variant::NM_001322970.2:c.1311G>C:missense variant::NM_001322967.2:c.1311G>C:missense variant::NM_000274.4:c.1311G>C:missense variant::NM_001322971.2:c.990G>C:missense variant	NM_001322965.2:c.1311G>C::NM_001171814.2:c.897G>C::NM_001322968.2:c.1311G>C::NM_001322970.2:c.1311G>C::NM_001322967.2:c.1311G>C::NM_001322974.2:c.711G>C::NM_001322969.2:c.1311G>C::NM_001322971.2:c.990G>C::NM_001322966.2:c.1311G>C::NM_000274.4:c.1311G>C	NP_001309898.1:p.Leu437Phe::NP_001165285.1:p.Leu299Phe::NP_001309903.1:p.Leu237Phe::NP_000265.1:p.Leu437Phe::NP_001309899.1:p.Leu437Phe::P04181:p.Leu437Phe::NP_001309900.1:p.Leu330Phe::NP_001309894.1:p.Leu437Phe::NP_001309895.1:p.Leu437Phe::NP_001309896.1:p.Leu437Phe::NP_001309897.1:p.Leu437Phe	NC_000010.11:g.124397951C>G::LRG_685:g.26000G>C::NG_008861.1:g.26000G>C	NC_000010.10:g.126086520C>G	L299F::L237F::L437F::L330F	GRCh38	10	124397951	124397951	124397951	C	G	GRCh37	10	126086520	126086520	126086520	C	G	-	-	-	-	-	-	-	4942	OAT	P04181#VAR_000586	rs1800456	613349.0013	Genetic Alliance:Ornithine+aminotransferase+deficiency/5411::OMIM:258870::MedGen:C0018425::Orphanet:414::SNOMED CT:276426004::MONDO:MONDO:0009796	Ornithine aminotransferase deficiency	-	SCV000020324	Benign	no assertion criteria provided	1737786	OMIM	1992-02-15	2018-02-01	OAT	benign	-	GRCh38	10	124397951	C	G	GRCh37	10	126086520	C	G	10-124397951-C-G	10-126086520-C-G\n"
        "85346099	NM_016341.4(PLCE1):c.3346C>T (p.Arg1116Ter) AND Nephrotic syndrome, type 3	RCV000002437	2012-08-13	2022-02-19	2006-12-01	no assertion criteria provided	0	Pathogenic	-	germline	-	34436835::17086182	-	-	Variant	VCV000002346	single nucleotide variant	17385	NM_016341.4:c.3346C>T:nonsense::NM_001165979.2:c.2422C>T:nonsense::NM_001288989.2:c.3346C>T:nonsense	NM_016341.4:c.3346C>T::NM_001288989.2:c.3346C>T::NM_001165979.2:c.2422C>T	NP_057425.3:p.Arg1116Ter::NP_001159451.1:p.Arg808Ter::NP_001275918.1:p.Arg1116Ter	NC_000010.11:g.94254256C>T::NG_015799.1:g.265268C>T	NC_000010.10:g.96014013C>T	R808*::R1116*	GRCh38	10	94254256	94254256	94254256	C	T	GRCh37	10	96014013	96014013	96014013	C	T	-	-	-	-	-	-	-	51196	PLCE1	-	rs121912602	608414.0003	Orphanet:656::MedGen:C1853124::Genetic Alliance:Nephrotic+syndrome%2C+type+3/8988::OMIM:610725::MONDO:MONDO:0012546	Nephrotic syndrome, type 3	-	SCV000022595	Pathogenic	no assertion criteria provided	17086182	OMIM	2006-12-01	2016-04-26	PLCE1	pathogenic	-	GRCh38	10	94254256	C	T	GRCh37	10	96014013	C	T	10-94254256-C-T	10-96014013-C-T\n"
        "85346098	NM_016341.4(PLCE1):c.1477C>T (p.Arg493Ter) AND Nephrotic syndrome, type 3	RCV000002436	2012-08-13	2022-02-19	2006-12-01	no assertion criteria provided	0	Pathogenic	-	germline	-	34436835::17086182	-	-	Variant	VCV000002345	single nucleotide variant	17384	NM_001165979.2:c.553C>T:nonsense::NM_001288989.2:c.1477C>T:nonsense::NM_016341.4:c.1477C>T:nonsense	NM_001165979.2:c.553C>T::NM_001288989.2:c.1477C>T::NM_016341.4:c.1477C>T	NP_001275918.1:p.Arg493Ter::NP_057425.3:p.Arg493Ter::NP_001159451.1:p.Arg185Ter	NC_000010.11:g.94132444C>T::NG_015799.1:g.143456C>T	NC_000010.10:g.95892201C>T	R493*::R185*	GRCh38	10	94132444	94132444	94132444	C	T	GRCh37	10	95892201	95892201	95892201	C	T	-	-	-	-	-	-	-	51196	PLCE1	-	rs121912601	608414.0002	Orphanet:656::MedGen:C1853124::Genetic Alliance:Nephrotic+syndrome%2C+type+3/8988::OMIM:610725::MONDO:MONDO:0012546	Nephrotic syndrome, type 3	-	SCV000022594	Pathogenic	no assertion criteria provided	17086182	OMIM	2006-12-01	2016-04-26	PLCE1	pathogenic	-	GRCh38	10	94132444	C	T	GRCh37	10	95892201	C	T	10-94132444-C-T	10-95892201-C-T\n"
        "85401547	NM_004826.4(ECEL1):c.2278T>C (p.Cys760Arg) AND Distal arthrogryposis type 5D	RCV000087050	2014-02-24	2022-02-19	2019-07-12	no assertion criteria provided	0	Pathogenic/Likely pathogenic	-	germline	-	23236030::2323603	-	-	Variant	VCV000100650	single nucleotide variant	106522	NM_001290787.2:c.2272T>C:missense variant::NM_004826.4:c.2278T>C:missense variant	NM_004826.4:c.2278T>C::NM_001290787.2:c.2272T>C	O95672:p.Cys760Arg::NP_001277716.1:p.Cys758Arg::NP_004817.2:p.Cys760Arg	NG_034065.1:g.12657T>C::NC_000002.12:g.232480203A>G	NC_000002.11:g.233344913A>G::NM_004826.2:c.2278T>C	C758R::C760R	GRCh38	2	232480203	232480203	232480203	A	G	GRCh37	2	233344913	233344913	233344913	A	G	-	-	-	-	-	-	-	9427	ECEL1	O95672#VAR_069995	rs587777129	605896.0007	OMIM:615065::MedGen:C3554415::MONDO:MONDO:0014028::Office of Rare Diseases:13059::Orphanet:329457	Distal arthrogryposis type 5D	-	SCV000119864||SCV002024449	Pathogenic||Likely pathogenic	no assertion criteria provided||no assertion criteria provided	23236030||-	OMIM||PerkinElmer Genomics	2013-04-15||2019-07-12	2015-06-11||2021-11-19	ECEL1||-	pathogenic::likely pathogenic	-	GRCh38	2	232480203	A	G	GRCh37	2	233344913	A	G	2-232480203-A-G	2-233344913-A-G\n"
    ),
)
def test_parse_clinvar_db(mock_open):

    clinvar_col2idx: dict = {
        "pathogenicity": 0,
        "consequence_dic": 1,
        "aa_change": 2,
    }

    expected = (
        {
            "CYP2C19": {"10-96521657-C-A": ["drug response", {"-": 1}, "-"]},
            "LOC110599570": {
                "10-96521657-C-A": ["drug response", {"-": 1}, "-"]
            },
            "OAT": {
                "10-126086520-C-G": [
                    "Benign",
                    {
                        "missense variant": 8,
                        "nonsense": 1,
                        "frameshift variant": 1,
                    },
                    "L299F::L237F::L437F::L330F",
                ]
            },
            "PLCE1": {
                "10-96014013-C-T": [
                    "Pathogenic",
                    {"nonsense": 3},
                    "R808*::R1116*",
                ],
                "10-95892201-C-T": [
                    "Pathogenic",
                    {"nonsense": 3},
                    "R493*::R185*",
                ],
            },
            "ECEL1": {
                "2-233344913-A-G": [
                    "Pathogenic/Likely pathogenic",
                    {"missense variant": 2},
                    "C758R::C760R",
                ]
            },
        },
        clinvar_col2idx,
    )

    assert expected == parse_clinvar_db("some_path.txt")


@pytest.mark.parametrize(
    "row, expected",
    [
        (
            [
                "74627773",
                "NM_000769.4(CYP2C19):c.-806C>A AND clopidogrel response - Dosage, Efficacy, Toxicity/ADR",
                "RCV000211201",
                "2016-05-18",
                "2021-09-29",
                "2016-06-14",
                "reviewed by expert panel",
                "3",
                "drug response",
                "-",
                "germline",
                "-",
                "22955794::23364775::20801498::22704413::23922007::21392617::21716271::20826260::23726091::22190063::20083681::22228204::22990067::24019397::19463375::22028352::23809542::22462746::22992668::20492469",
                "-",
                "-",
                "Variant",
                "VCV000225946",
                "single nucleotide variant",
                "227770",
                "-",
                "-",
                "-",
                "NC_000010.11:g.94761900C>A::NG_008384.3:g.4220C>A::NG_055436.1:g.1260C>A",
                "NC_000010.10:g.96521657C>A",
                "-",
                "GRCh38",
                "10",
                "94761900",
                "94761900",
                "94761900",
                "C",
                "A",
                "GRCh37",
                "10",
                "96521657",
                "96521657",
                "96521657",
                "C",
                "A",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "1557::110599570",
                "CYP2C19::LOC110599570",
                "-",
                "rs12248560",
                "-",
                "MedGen:CN236507",
                "clopidogrel response - Dosage, Efficacy, Toxicity/ADR",
                "-",
                "SCV000268179",
                "drug response",
                "reviewed by expert panel",
                "22955794::23364775::20801498::22704413::23922007::21392617::20826260::23726091::22190063::20083681::22228204::22990067::24019397::19463375::22028352::23809542::22462746::22992668::20492469",
                "PharmGKB",
                "2016-06-14",
                "2018-06-18",
                "CYP2C19",
                "-",
                "-",
                "GRCh38",
                "10",
                "94761900",
                "C",
                "A",
                "GRCh37",
                "10",
                "96521657",
                "C",
                "A",
                "10-94761900-C-A",
                "10-96521657-C-A",
            ],
            {"10-96521657-C-A": ["drug response", {"-": 1}, "-"]},
        ),
        (
            [
                "85343921",
                "NM_000274.4(OAT):c.1311G>C (p.Leu437Phe) AND Ornithine aminotransferase deficiency",
                "RCV000000181",
                "2012-08-13",
                "2022-02-19",
                "1992-02-15",
                "no assertion criteria provided",
                "0",
                "Benign",
                "-",
                "germline",
                "-",
                "1737786",
                "-",
                "-",
                "Variant",
                "VCV000000158",
                "single nucleotide variant",
                "15197",
                "NM_001322968.2:c.1311G>C:nonsense::NM_001322974.2:c.711G>C:frameshift variant::NM_001322966.2:c.1311G>C:missense variant::NM_001322965.2:c.1311G>C:missense variant::NM_001171814.2:c.897G>C:missense variant::NM_001322969.2:c.1311G>C:missense variant::NM_001322970.2:c.1311G>C:missense variant::NM_001322967.2:c.1311G>C:missense variant::NM_000274.4:c.1311G>C:missense variant::NM_001322971.2:c.990G>C:missense variant",
                "NM_001322965.2:c.1311G>C::NM_001171814.2:c.897G>C::NM_001322968.2:c.1311G>C::NM_001322970.2:c.1311G>C::NM_001322967.2:c.1311G>C::NM_001322974.2:c.711G>C::NM_001322969.2:c.1311G>C::NM_001322971.2:c.990G>C::NM_001322966.2:c.1311G>C::NM_000274.4:c.1311G>C",
                "NP_001309898.1:p.Leu437Phe::NP_001165285.1:p.Leu299Phe::NP_001309903.1:p.Leu237Phe::NP_000265.1:p.Leu437Phe::NP_001309899.1:p.Leu437Phe::P04181:p.Leu437Phe::NP_001309900.1:p.Leu330Phe::NP_001309894.1:p.Leu437Phe::NP_001309895.1:p.Leu437Phe::NP_001309896.1:p.Leu437Phe::NP_001309897.1:p.Leu437Phe",
                "NC_000010.11:g.124397951C>G::LRG_685:g.26000G>C::NG_008861.1:g.26000G>C",
                "NC_000010.10:g.126086520C>G",
                "L299F::L237F::L437F::L330F",
                "GRCh38",
                "10",
                "124397951",
                "124397951",
                "124397951",
                "C",
                "G",
                "GRCh37",
                "10",
                "126086520",
                "126086520",
                "126086520",
                "C",
                "G",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "4942",
                "OAT",
                "P04181#VAR_000586",
                "rs1800456",
                "613349.0013",
                "Genetic Alliance:Ornithine+aminotransferase+deficiency/5411::OMIM:258870::MedGen:C0018425::Orphanet:414::SNOMED CT:276426004::MONDO:MONDO:0009796",
                "Ornithine aminotransferase deficiency",
                "-",
                "SCV000020324",
                "Benign",
                "no assertion criteria provided",
                "1737786",
                "OMIM",
                "1992-02-15",
                "2018-02-01",
                "OAT",
                "benign",
                "-",
                "GRCh38",
                "10",
                "124397951",
                "C",
                "G",
                "GRCh37",
                "10",
                "126086520",
                "C",
                "G",
                "10-124397951-C-G",
                "10-126086520-C-G",
            ],
            {
                "10-126086520-C-G": [
                    "Benign",
                    {
                        "missense variant": 8,
                        "nonsense": 1,
                        "frameshift variant": 1,
                    },
                    "L299F::L237F::L437F::L330F",
                ]
            },
        ),
    ],
)
def test_make_clinvar_gene_info_dic(row, expected):

    line = "#clinvar_set_id\ttitle\trcv\tcreation_date\tlast_updated\tlast_evaluated\treview_status\tstar\tpathogenicity\tinheritance\torigin\tindividuals\tPMIDs\tgenotype_type\tgenotype_vcv\tvariant:variation_type\tvariant:vcv\tvariant:variant_type\tvariant:allele_id\tvariant:molecular_consequence\tvariant:HGVSc\tvariant:HGVSp\tvariant:HGVSg\tvariant:HGVS_previous\tvariant:Protein1LetterCode\tvariant:GRCh38:Assembly\tvariant:GRCh38:Chr\tvariant:GRCh38:positionVCF\tvariant:GRCh38:start\tvariant:GRCh38:stop\tvariant:GRCh38:referenceAllele\tvariant:GRCh38:alternateAllele\tvariant:GRCh37:Assembly\tvariant:GRCh37:Chr\tvariant:GRCh37:positionVCF\tvariant:GRCh37:start\tvariant:GRCh37:stop\tvariant:GRCh37:referenceAllele\tvariant:GRCh37:alternateAllele\tvariant:NCBI36:Assembly\tvariant:NCBI36:Chr\tvariant:NCBI36:positionVCF\tvariant:NCBI36:start\tvariant:NCBI36:stop\tvariant:NCBI36:referenceAllele\tvariant:NCBI36:alternateAllele\tvariant:gene:id\tvariant:gene:symbol\tvariant:uniprot\tvariant:dbsnp\tvariant:omim_gene_id\tdisease_id\tdisease_name\tdisease_mechanism\tscv:scv\tscv:pathogenicity\tscv:review_status\tscv:pmid\tscv:submitter\tscv:last_evaluated\tscv:submit_date\tscv:gene_symbol\tscv_pathogenicity\trcv_pathogenicity\tnormalized_variant:GRCh38:Assembly\tnormalized_variant:GRCh38:Chr\tnormalized_variant:GRCh38:positionVCF\tnormalized_variant:GRCh38:referenceAllele\tnormalized_variant:GRCh38:alternateAllele\tnormalized_variant:GRCh37:Assembly\tnormalized_variant:GRCh37:Chr\tnormalized_variant:GRCh37:positionVCF\tnormalized_variant:GRCh37:referenceAllele\tnormalized_variant:GRCh37:alternateAllele\tnormalized_variant:GRCh38:CPRA\tnormalized_variant:GRCh37:CPRA\n"
    f_col2idx = {
        val: idx for idx, val in enumerate(line.strip("#").strip().split("\t"))
    }

    assert expected == make_clinvar_gene_info_dic(row, f_col2idx)


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        "#omimPhenoId\tomimGeneId\torphaId\torphaGeneSymbol\tncbiGeneId\tensemblGeneId\thgncId\ttitle\tpreferredTitle\tgeneSymbol\tinheritances:value\tinheritances:source\tonsetAges:value\tonsetAges:source\tprevalences:value\tprevalences:source\tsymptoms:id\tsymptoms:name\tsymptoms:source\tomimFlag\torphaExtraInfo\tomimPMID\talternativeGeneSymbols\tpheno_geno_source\n"
        '-	-	ORPHA:157949	INO80::RAG1::RAG2	-	-	-	Combined immunodeficiency with granulomatosis	-	INO80::RAG1::RAG2	-	-	-	-	-	-	-	-	-	-	{"name": "Combined immunodeficiency with granulomatosis", "disorder_type": "Disease", "external_references": {"ICD-10": [{"reference": "D81.1", "label": "NTBT (ORPHA codes Narrower Term maps to a Broader Term)"}], "OMIM": [{"reference": "233650", "label": "E (Exact mapping: the two concepts are equivalent)"}], "UMLS": [{"reference": "C2673536", "label": "E (Exact mapping: the two concepts are equivalent)"}]}, "genes": {"INO80": {"source": "25312759[PMID]", "name": "INO80 complex ATPase subunit", "synonyms": ["INO80 complex subunit A", "INO80A", "KIAA1259", "hINO80"], "references": {"Genatlas": "INO80", "Ensembl": "ENSG00000128908", "HGNC": "26956", "Reactome": "Q9ULG1", "OMIM": "610169", "SwissProt": "Q9ULG1"}}, "RAG1": {"source": "18463379[PMID]", "name": "recombination activating 1", "synonyms": ["MGC43321", "RING finger protein 74", "RNF74", "V(D)J recombination-activating protein 1", "recombination activating protein 1"], "references": {"Ensembl": "ENSG00000166349", "Genatlas": "RAG1", "HGNC": "9831", "OMIM": "179615", "Reactome": "P15918", "SwissProt": "P15918"}}, "RAG2": {"source": "18463379[PMID]", "name": "recombination activating 2", "synonyms": [], "references": {"Ensembl": "ENSG00000175097", "Genatlas": "RAG2", "HGNC": "9832", "OMIM": "179616", "Reactome": "P55895", "SwissProt": "P55895"}}}}	-	INO80::RAG1::RAG2	ORPHA\n'
        "603554	179615	-	-	5896	ENSG00000166349	9831	Omenn syndrome	OMENN SYNDROME	RAG1	Autosomal recessive	OMIM::CGD::HPO	Pediatric	CGD	-	-	HP:0002240||HP:0001880||HP:0000007||HP:0001873||HP:0001744||HP:0001903||HP:0002716||HP:0001596||HP:0002293||HP:0001019||HP:0002014||HP:0003075||HP:0003212||HP:0001508||HP:0002090||HP:0001072||HP:0005365||HP:0004324||HP:0002028||HP:0200034||HP:0001425||HP:0005374||HP:0011123||HP:0008940||HP:0002665||HP:0000100||HP:0100806||HP:0001945||HP:0031273||HP:0005871||HP:0004385||HP:0001433||HP:0008873||HP:0000988||HP:0010976||HP:0002718||HP:0000778||HP:0004429||HP:0002841	Hepatomegaly||Eosinophilia||Autosomal recessive inheritance||Thrombocytopenia||Splenomegaly||Anemia||Lymphadenopathy||Alopecia||Alopecia of scalp||Erythroderma||Diarrhea||Hypoproteinemia||Increased circulating IgE level||Failure to thrive||Pneumonia||Thickened skin||Severe B lymphocytopenia||Increased body weight||Chronic diarrhea||Papule||Heterogeneous||Cellular immunodeficiency||Inflammatory abnormality of the skin||Generalized lymphadenopathy||Lymphoma||Nephrotic syndrome||Sepsis||Fever||Shock||Metaphyseal chondrodysplasia||Protracted diarrhea||Hepatosplenomegaly||Disproportionate short-limb short stature||Skin rash||B lymphocytopenia||Recurrent bacterial infections||Hypoplasia of the thymus||Recurrent viral infections||Recurrent fungal infections	OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO||OMIM::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::MANUAL::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||HPO||HPO::JAY_O||HPO::JAY_O||HPO::JAY_O||HPO::JAY_O	#	-	2248640::4115568::15696198::8450050::6243521::11313270::2010548::15731174::5809843::6326896::7608815::3879354::3679186::17476359::17476358::14328107::7050708::8776375::1986108::9630231	RAG1	OMIM::MANUAL::CGD\n"
        "-	608843	-	-	114990	ENSG00000168140	18517	-	-	VASN	-	-	-	-	-	-	-	-	-	*	-	-	VASN	OMIM\n"
        "136760	606014	ORPHA:391474	ALX3	257	ENSG00000156150	449	Frontonasal dysplasia 1	FRONTONASAL DYSPLASIA 1; FND1	ALX3	Autosomal recessive||Sporadic	OMIM::ORPHA::CGD::HPO||HPO	Neonatal	ORPHA	Unknown__Worldwide	ORPHA	HP:0006992||HP:0001156||HP:0030084||HP:0000369||HP:0040019||HP:0000286||HP:0006931||HP:0000431||HP:0007541||HP:0001274||HP:0000568||HP:0000518||HP:0005258||HP:0030024||HP:0000384||HP:0000589||HP:0009099||HP:0000327||HP:0001249||HP:0001636||HP:0002738||HP:0000349||HP:0000007||HP:0000508||HP:0000405||HP:0012385||HP:0000316||HP:0000161||HP:0004423||HP:0000175||HP:0000368||HP:0000486||HP:0000612||HP:0002084||HP:0002650||HP:0002938||HP:0004112||HP:0007370||HP:0008591||HP:0010297||HP:0011817||HP:0025247||HP:0100490||HP:0000873||HP:0040075||HP:0000846||HP:0000365||HP:0001627||HP:0000445||HP:0000463||HP:0002000||HP:0002006||HP:0025514||HP:0010761||HP:0012032||HP:0002282||HP:0001425||HP:0001426||HP:0001012||HP:0001566||HP:0011803||HP:0003745||HP:0000006||HP:0000969||HP:0000456||HP:0000455||HP:0100259||HP:0000821||HP:0000824||HP:0000202||HP:0001004||HP:0030854||HP:0000343||HP:0000527||HP:0100629||HP:0009466||HP:0001162||HP:0009473	Anterior basal encephalocele||Brachydactyly||Clinodactyly||Low-set ears||Finger clinodactyly||Epicanthus||Pericallosal lipoma||Wide nasal bridge||Frontal cutaneous lipoma||Agenesis of corpus callosum||Microphthalmia||Cataract||Pectoral muscle hypoplasia/aplasia||Pretragal ectopia||Preauricular skin tag||Coloboma||Median cleft palate||Hypoplasia of the maxilla||Intellectual disability||Tetralogy of Fallot||Hypoplastic frontal sinuses||Widow's peak||Autosomal recessive inheritance||Ptosis||Conductive hearing impairment||Camptodactyly||Hypertelorism||Median cleft lip||Cranium bifidum occultum||Cleft palate||Low-set, posteriorly rotated ears||Strabismus||Iris coloboma||Encephalocele||Scoliosis||Lumbar hyperlordosis||Midline nasal groove||Aplasia/Hypoplasia of the corpus callosum||Congenital conductive hearing impairment||Bifid tongue||Basal encephalocele||Dermoid cyst||Camptodactyly of finger||Diabetes insipidus||Hypopituitarism||Adrenal insufficiency||Hearing impairment||Abnormal heart morphology||Wide nose||Anteverted nares||Short columella||Facial cleft||Morning glory anomaly||Broad columella||Lipoma||Gray matter heterotopia||Heterogeneous||Multifactorial inheritance||Multiple lipomas||Widely-spaced maxillary central incisors||Bifid nose||Sporadic||Autosomal dominant inheritance||Edema||Bifid nasal tip||Broad nasal tip||Postaxial polydactyly||Hypothyroidism||Decreased response to growth hormone stimulation test||Oral cleft||Lymphedema||Scleral staphyloma||Long philtrum||Long eyelashes||Midline facial cleft||Radial deviation of finger||Postaxial hand polydactyly||Joint contracture of the hand	OMIM::HPO::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::ORPHA::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::ORPHA::HPO::JAY_O||OMIM::ORPHA::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO||OMIM::ORPHA::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::ORPHA::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::ORPHA::MANUAL::HPO::JAY_O||ORPHA||ORPHA||ORPHA||ORPHA||ORPHA::MANUAL::JAY_O||ORPHA||ORPHA||ORPHA||ORPHA||ORPHA||ORPHA||ORPHA||ORPHA||ORPHA||ORPHA::MANUAL::JAY_O||ORPHA||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::HPO::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL||MANUAL||MANUAL::JAY_O||MANUAL::HPO::JAY_O||MANUAL::HPO::JAY_O||MANUAL::HPO||MANUAL||MANUAL::JAY_O||MANUAL::HPO::JAY_O||MANUAL::HPO::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||HPO::JAY_O||HPO::JAY_O||HPO::JAY_O	#	-	20106874::3560167::13206986::8362915::19365836::9689987::17963218::8741108::15127764::2738904::15384079::7363499::10564879::5444583::2840620::7762593::4003439::19409524::17955515	ALX3::FND1	OMIM::ORPHA::MANUAL::CGD\n"
    ),
)
def test_parse_disease_db(mock_open):

    disease_col2idx: dict = {
        "title": 0,
        "inheritance": 1,
        "onsetAges": 2,
        "symtoms_id": 3,
        "symtoms": 4,
    }

    expected = (
        {
            "INO80": [
                [
                    "Combined immunodeficiency with granulomatosis",
                    ["-"],
                    ["-"],
                    ["-"],
                    ["-"],
                ]
            ],
            "RAG1": [
                [
                    "Combined immunodeficiency with granulomatosis",
                    ["-"],
                    ["-"],
                    ["-"],
                    ["-"],
                ],
                [
                    "Omenn syndrome",
                    ["Autosomal recessive"],
                    ["Pediatric"],
                    [
                        "HP:0002240",
                        "HP:0001880",
                        "HP:0000007",
                        "HP:0001873",
                        "HP:0001744",
                        "HP:0001903",
                        "HP:0002716",
                        "HP:0001596",
                        "HP:0002293",
                        "HP:0001019",
                        "HP:0002014",
                        "HP:0003075",
                        "HP:0003212",
                        "HP:0001508",
                        "HP:0002090",
                        "HP:0001072",
                        "HP:0005365",
                        "HP:0004324",
                        "HP:0002028",
                        "HP:0200034",
                        "HP:0001425",
                        "HP:0005374",
                        "HP:0011123",
                        "HP:0008940",
                        "HP:0002665",
                        "HP:0000100",
                        "HP:0100806",
                        "HP:0001945",
                        "HP:0031273",
                        "HP:0005871",
                        "HP:0004385",
                        "HP:0001433",
                        "HP:0008873",
                        "HP:0000988",
                        "HP:0010976",
                        "HP:0002718",
                        "HP:0000778",
                        "HP:0004429",
                        "HP:0002841",
                    ],
                    [
                        "Hepatomegaly",
                        "Eosinophilia",
                        "Autosomal recessive inheritance",
                        "Thrombocytopenia",
                        "Splenomegaly",
                        "Anemia",
                        "Lymphadenopathy",
                        "Alopecia",
                        "Alopecia of scalp",
                        "Erythroderma",
                        "Diarrhea",
                        "Hypoproteinemia",
                        "Increased circulating IgE level",
                        "Failure to thrive",
                        "Pneumonia",
                        "Thickened skin",
                        "Severe B lymphocytopenia",
                        "Increased body weight",
                        "Chronic diarrhea",
                        "Papule",
                        "Heterogeneous",
                        "Cellular immunodeficiency",
                        "Inflammatory abnormality of the skin",
                        "Generalized lymphadenopathy",
                        "Lymphoma",
                        "Nephrotic syndrome",
                        "Sepsis",
                        "Fever",
                        "Shock",
                        "Metaphyseal chondrodysplasia",
                        "Protracted diarrhea",
                        "Hepatosplenomegaly",
                        "Disproportionate short-limb short stature",
                        "Skin rash",
                        "B lymphocytopenia",
                        "Recurrent bacterial infections",
                        "Hypoplasia of the thymus",
                        "Recurrent viral infections",
                        "Recurrent fungal infections",
                    ],
                ],
            ],
            "RAG2": [
                [
                    "Combined immunodeficiency with granulomatosis",
                    ["-"],
                    ["-"],
                    ["-"],
                    ["-"],
                ]
            ],
            "VASN": [["-", ["-"], ["-"], ["-"], ["-"]]],
            "ALX3": [
                [
                    "Frontonasal dysplasia 1",
                    ["Autosomal recessive", "Sporadic"],
                    ["Neonatal"],
                    [
                        "HP:0006992",
                        "HP:0001156",
                        "HP:0030084",
                        "HP:0000369",
                        "HP:0040019",
                        "HP:0000286",
                        "HP:0006931",
                        "HP:0000431",
                        "HP:0007541",
                        "HP:0001274",
                        "HP:0000568",
                        "HP:0000518",
                        "HP:0005258",
                        "HP:0030024",
                        "HP:0000384",
                        "HP:0000589",
                        "HP:0009099",
                        "HP:0000327",
                        "HP:0001249",
                        "HP:0001636",
                        "HP:0002738",
                        "HP:0000349",
                        "HP:0000007",
                        "HP:0000508",
                        "HP:0000405",
                        "HP:0012385",
                        "HP:0000316",
                        "HP:0000161",
                        "HP:0004423",
                        "HP:0000175",
                        "HP:0000368",
                        "HP:0000486",
                        "HP:0000612",
                        "HP:0002084",
                        "HP:0002650",
                        "HP:0002938",
                        "HP:0004112",
                        "HP:0007370",
                        "HP:0008591",
                        "HP:0010297",
                        "HP:0011817",
                        "HP:0025247",
                        "HP:0100490",
                        "HP:0000873",
                        "HP:0040075",
                        "HP:0000846",
                        "HP:0000365",
                        "HP:0001627",
                        "HP:0000445",
                        "HP:0000463",
                        "HP:0002000",
                        "HP:0002006",
                        "HP:0025514",
                        "HP:0010761",
                        "HP:0012032",
                        "HP:0002282",
                        "HP:0001425",
                        "HP:0001426",
                        "HP:0001012",
                        "HP:0001566",
                        "HP:0011803",
                        "HP:0003745",
                        "HP:0000006",
                        "HP:0000969",
                        "HP:0000456",
                        "HP:0000455",
                        "HP:0100259",
                        "HP:0000821",
                        "HP:0000824",
                        "HP:0000202",
                        "HP:0001004",
                        "HP:0030854",
                        "HP:0000343",
                        "HP:0000527",
                        "HP:0100629",
                        "HP:0009466",
                        "HP:0001162",
                        "HP:0009473",
                    ],
                    [
                        "Anterior basal encephalocele",
                        "Brachydactyly",
                        "Clinodactyly",
                        "Low-set ears",
                        "Finger clinodactyly",
                        "Epicanthus",
                        "Pericallosal lipoma",
                        "Wide nasal bridge",
                        "Frontal cutaneous lipoma",
                        "Agenesis of corpus callosum",
                        "Microphthalmia",
                        "Cataract",
                        "Pectoral muscle hypoplasia/aplasia",
                        "Pretragal ectopia",
                        "Preauricular skin tag",
                        "Coloboma",
                        "Median cleft palate",
                        "Hypoplasia of the maxilla",
                        "Intellectual disability",
                        "Tetralogy of Fallot",
                        "Hypoplastic frontal sinuses",
                        "Widow's peak",
                        "Autosomal recessive inheritance",
                        "Ptosis",
                        "Conductive hearing impairment",
                        "Camptodactyly",
                        "Hypertelorism",
                        "Median cleft lip",
                        "Cranium bifidum occultum",
                        "Cleft palate",
                        "Low-set, posteriorly rotated ears",
                        "Strabismus",
                        "Iris coloboma",
                        "Encephalocele",
                        "Scoliosis",
                        "Lumbar hyperlordosis",
                        "Midline nasal groove",
                        "Aplasia/Hypoplasia of the corpus callosum",
                        "Congenital conductive hearing impairment",
                        "Bifid tongue",
                        "Basal encephalocele",
                        "Dermoid cyst",
                        "Camptodactyly of finger",
                        "Diabetes insipidus",
                        "Hypopituitarism",
                        "Adrenal insufficiency",
                        "Hearing impairment",
                        "Abnormal heart morphology",
                        "Wide nose",
                        "Anteverted nares",
                        "Short columella",
                        "Facial cleft",
                        "Morning glory anomaly",
                        "Broad columella",
                        "Lipoma",
                        "Gray matter heterotopia",
                        "Heterogeneous",
                        "Multifactorial inheritance",
                        "Multiple lipomas",
                        "Widely-spaced maxillary central incisors",
                        "Bifid nose",
                        "Sporadic",
                        "Autosomal dominant inheritance",
                        "Edema",
                        "Bifid nasal tip",
                        "Broad nasal tip",
                        "Postaxial polydactyly",
                        "Hypothyroidism",
                        "Decreased response to growth hormone stimulation test",
                        "Oral cleft",
                        "Lymphedema",
                        "Scleral staphyloma",
                        "Long philtrum",
                        "Long eyelashes",
                        "Midline facial cleft",
                        "Radial deviation of finger",
                        "Postaxial hand polydactyly",
                        "Joint contracture of the hand",
                    ],
                ]
            ],
        },
        disease_col2idx,
    )

    assert expected == parse_disease_db("some_path.txt")


@pytest.mark.parametrize(
    "row, expected",
    [
        (
            [
                "-",
                "608843",
                "-",
                "-",
                "114990",
                "ENSG00000168140",
                "18517",
                "-",
                "-",
                "VASN",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "*",
                "-",
                "-",
                "VASN",
                "OMIM",
            ],
            ["-", ["-"], ["-"], ["-"], ["-"]],
        ),
        (
            [
                "603554",
                "179615",
                "-",
                "-",
                "5896",
                "ENSG00000166349",
                "9831",
                "Omenn syndrome",
                "OMENN SYNDROME",
                "RAG1",
                "Autosomal recessive",
                "OMIM::CGD::HPO",
                "Pediatric",
                "CGD",
                "-",
                "-",
                "HP:0002240||HP:0001880||HP:0000007||HP:0001873||HP:0001744||HP:0001903||HP:0002716||HP:0001596||HP:0002293||HP:0001019||HP:0002014||HP:0003075||HP:0003212||HP:0001508||HP:0002090||HP:0001072||HP:0005365||HP:0004324||HP:0002028||HP:0200034||HP:0001425||HP:0005374||HP:0011123||HP:0008940||HP:0002665||HP:0000100||HP:0100806||HP:0001945||HP:0031273||HP:0005871||HP:0004385||HP:0001433||HP:0008873||HP:0000988||HP:0010976||HP:0002718||HP:0000778||HP:0004429||HP:0002841",
                "Hepatomegaly||Eosinophilia||Autosomal recessive inheritance||Thrombocytopenia||Splenomegaly||Anemia||Lymphadenopathy||Alopecia||Alopecia of scalp||Erythroderma||Diarrhea||Hypoproteinemia||Increased circulating IgE level||Failure to thrive||Pneumonia||Thickened skin||Severe B lymphocytopenia||Increased body weight||Chronic diarrhea||Papule||Heterogeneous||Cellular immunodeficiency||Inflammatory abnormality of the skin||Generalized lymphadenopathy||Lymphoma||Nephrotic syndrome||Sepsis||Fever||Shock||Metaphyseal chondrodysplasia||Protracted diarrhea||Hepatosplenomegaly||Disproportionate short-limb short stature||Skin rash||B lymphocytopenia||Recurrent bacterial infections||Hypoplasia of the thymus||Recurrent viral infections||Recurrent fungal infections",
                "OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO||OMIM::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::MANUAL::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||OMIM::MANUAL::HPO::JAY_O||OMIM::HPO::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||MANUAL::JAY_O||HPO||HPO::JAY_O||HPO::JAY_O||HPO::JAY_O||HPO::JAY_O",
                "#",
                "-",
                "2248640::4115568::15696198::8450050::6243521::11313270::2010548::15731174::5809843::6326896::7608815::3879354::3679186::17476359::17476358::14328107::7050708::8776375::1986108::9630231",
                "RAG1",
                "OMIM::MANUAL::CGD",
            ],
            [
                "Omenn syndrome",
                ["Autosomal recessive"],
                ["Pediatric"],
                [
                    "HP:0002240",
                    "HP:0001880",
                    "HP:0000007",
                    "HP:0001873",
                    "HP:0001744",
                    "HP:0001903",
                    "HP:0002716",
                    "HP:0001596",
                    "HP:0002293",
                    "HP:0001019",
                    "HP:0002014",
                    "HP:0003075",
                    "HP:0003212",
                    "HP:0001508",
                    "HP:0002090",
                    "HP:0001072",
                    "HP:0005365",
                    "HP:0004324",
                    "HP:0002028",
                    "HP:0200034",
                    "HP:0001425",
                    "HP:0005374",
                    "HP:0011123",
                    "HP:0008940",
                    "HP:0002665",
                    "HP:0000100",
                    "HP:0100806",
                    "HP:0001945",
                    "HP:0031273",
                    "HP:0005871",
                    "HP:0004385",
                    "HP:0001433",
                    "HP:0008873",
                    "HP:0000988",
                    "HP:0010976",
                    "HP:0002718",
                    "HP:0000778",
                    "HP:0004429",
                    "HP:0002841",
                ],
                [
                    "Hepatomegaly",
                    "Eosinophilia",
                    "Autosomal recessive inheritance",
                    "Thrombocytopenia",
                    "Splenomegaly",
                    "Anemia",
                    "Lymphadenopathy",
                    "Alopecia",
                    "Alopecia of scalp",
                    "Erythroderma",
                    "Diarrhea",
                    "Hypoproteinemia",
                    "Increased circulating IgE level",
                    "Failure to thrive",
                    "Pneumonia",
                    "Thickened skin",
                    "Severe B lymphocytopenia",
                    "Increased body weight",
                    "Chronic diarrhea",
                    "Papule",
                    "Heterogeneous",
                    "Cellular immunodeficiency",
                    "Inflammatory abnormality of the skin",
                    "Generalized lymphadenopathy",
                    "Lymphoma",
                    "Nephrotic syndrome",
                    "Sepsis",
                    "Fever",
                    "Shock",
                    "Metaphyseal chondrodysplasia",
                    "Protracted diarrhea",
                    "Hepatosplenomegaly",
                    "Disproportionate short-limb short stature",
                    "Skin rash",
                    "B lymphocytopenia",
                    "Recurrent bacterial infections",
                    "Hypoplasia of the thymus",
                    "Recurrent viral infections",
                    "Recurrent fungal infections",
                ],
            ],
        ),
    ],
)
def test_make_disease_info_list(row, expected):

    line = "#omimPhenoId\tomimGeneId\torphaId\torphaGeneSymbol\tncbiGeneId\tensemblGeneId\thgncId\ttitle\tpreferredTitle\tgeneSymbol\tinheritances:value\tinheritances:source\tonsetAges:value\tonsetAges:source\tprevalences:value\tprevalences:source\tsymptoms:id\tsymptoms:name\tsymptoms:source\tomimFlag\torphaExtraInfo\tomimPMID\talternativeGeneSymbols\tpheno_geno_source\n"
    f_col2idx = {
        val: idx for idx, val in enumerate(line.strip("#").strip().split("\t"))
    }

    assert expected == make_disease_info_list(row, f_col2idx)


@mock.patch(
    "builtins.open",
    new_callable=mock.mock_open,
    read_data=(
        "   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat\n"
        "score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID\n"
        "\n"
        "  463   1.3  0.6  1.7  chr1        10001   10468 (249240153) +  (TAACCC)n      Simple_repeat            1  463    (0)      1\n"
        " 3612  11.4 21.5  1.3  chr1        10469   11447 (249239174) C  TAR1           Satellite/telo       (399) 1712    483      2\n"
        " 2597  18.4  2.7  3.3  chr4      184488284 184488398 (6665878) C  MLT2B2         LTR/ERVL               (0)  515    402 1297781\n"
    ),
)
def test_parse_repeatmasker_db(mock_open):

    expected = {
        "1": [(10001, 10468), (10469, 11447)],
        "4": [(184488284, 184488398)],
    }

    assert expected == parse_repeatmasker_db("some_path.fa.out")


######################################################################


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
