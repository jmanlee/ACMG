# module for compare EVIDENCE results.

from collections import defaultdict

EVIDENCE_3B = "/data/projects/ACMG/output/ETD21-RWNY.acmg_assigned.tsv"
EVIDENCE_JM = "/data/projects/ACMG/output/ACMG_result_0520_2.txt"


def parse_evidence(filename: str) -> dict:
    """_summary_
    Note:
        variant_id - feature: [EVIDENCE RULE] 
        의 형태로 데이터를 가공한다. 
    """

    with open(filename) as infile:

        evidence_result_dic = defaultdict(dict)

        for line in infile:
            if line.startswith("#"):
                f_col2idx = {
                    val: idx
                    for idx, val in enumerate(
                        line.strip("#").strip().split("\t")
                    )
                }
            if not line.startswith("#"):
                row = line.strip().split("\t")
                var_id = row[f_col2idx["Uploaded_variation"]]
                var_feature = row[f_col2idx["Feature"]]
                acmg_rule = row[f_col2idx["ACMG_RULE"]]
                acmg_rules: list = acmg_rule.split("||")

                evidence_result_dic[var_id].update({var_feature: acmg_rules})

    return evidence_result_dic


def parse_my_evidence(filename: str) -> dict:
    """_summary_
    Note:
        variant_id - feature: [My RULE] 
        의 형태로 데이터를 가공한다. 
    """

    with open(filename) as infile:

        my_evidence_result_dic = defaultdict(dict)

        for line in infile:
            if line.startswith("#"):
                f_col2idx = {
                    val: idx
                    for idx, val in enumerate(
                        line.strip("#").strip().split("\t")
                    )
                }
            if not line.startswith("#"):
                row = line.strip().split("\t")
                var_id = row[f_col2idx["Variant_ID"]]
                var_feature = row[f_col2idx["Feature"]]
                acmg_rule = row[f_col2idx["ACMG_rule"]]
                acmg_rules: list = acmg_rule.split("||")

                my_evidence_result_dic[var_id].update({var_feature: acmg_rules})

    return my_evidence_result_dic


def compare_rules(evidence_result_dic: dict(), my_result_dic: dict()):
    """_summary_
    Note:
        My ACMG vs. EVIDENCE rule을 비교하여 비율을 계산한다. 
    
    Print:
        1. My rule / EVIDENCE
        2. EVIDENCE / My rule
    """

    ## pvs1, ps1, ps2, pm2, pm4, pm5, pp2, pp3, pp5 9개
    ## ba1, bs1, bp1, bp3, bp4, bp6, bp7 7개

    # [count, total count]
    rules_dic_1 = {
        "PVS1": [0, 0],
        "PS1": [0, 0],
        "PS2": [0, 0],
        "PM2": [0, 0],
        "PM4": [0, 0],
        "PM5": [0, 0],
        "PP2": [0, 0],
        "PP3": [0, 0],
        "PP5": [0, 0],
        "BA1": [0, 0],
        "BS1": [0, 0],
        "BP1": [0, 0],
        "BP3": [0, 0],
        "BP4": [0, 0],
        "BP6": [0, 0],
        "BP7": [0, 0],
    }

    for var_id in evidence_result_dic:
        for var_feature in evidence_result_dic[var_id]:
            evidences: list = evidence_result_dic[var_id][var_feature]

            # EVIDENCE가 기준인 경우, my_variant에서 key가 없더라도 total count에 포함.
            try:
                my_evidences: list = my_result_dic[var_id][var_feature]
            except KeyError:
                for acmg_rule, acmg_rule_counts in rules_dic_1.items():
                    for rule_3b in evidences:
                        if rule_3b.startswith(acmg_rule) and (
                            not rule_3b.endswith("-")
                        ):
                            acmg_rule_counts[1] += 1
                continue

            # EVIDENCE의 rule을 세고, 그 rule이 My-rule에도 있는지 counting
            for acmg_rule, acmg_rule_counts in rules_dic_1.items():
                for rule_3b in evidences:
                    if rule_3b.startswith(acmg_rule) and (
                        not rule_3b.endswith("-")
                    ):
                        acmg_rule_counts[1] += 1
                        for rule_my in my_evidences:
                            if rule_my.startswith(acmg_rule):
                                acmg_rule_counts[0] += 1

    for acmg_rule, acmg_rule_counts in rules_dic_1.items():
        try:
            print(
                f"{acmg_rule}: {acmg_rule_counts[0]/acmg_rule_counts[1]:.4f}, {acmg_rule_counts}"
            )
        except ZeroDivisionError:
            print(f"{acmg_rule}: None, {acmg_rule_counts}")

    # [count, total count]
    rules_dic_2 = {
        "PVS1": [0, 0],
        "PS1": [0, 0],
        "PS2": [0, 0],
        "PM2": [0, 0],
        "PM4": [0, 0],
        "PM5": [0, 0],
        "PP2": [0, 0],
        "PP3": [0, 0],
        "PP5": [0, 0],
        "BA1": [0, 0],
        "BS1": [0, 0],
        "BP1": [0, 0],
        "BP3": [0, 0],
        "BP4": [0, 0],
        "BP6": [0, 0],
        "BP7": [0, 0],
    }

    for var_id in my_result_dic:
        for var_feature in my_result_dic[var_id]:
            my_evidences: list = my_result_dic[var_id][var_feature]

            # My-rule이 기준인 경우, EVIENCE_variant에서 key가 없으면 제외.
            try:
                evidences: list = evidence_result_dic[var_id][var_feature]
            except KeyError:
                continue

            # My-rule을 세고, 그 rule이 EVIDENCE 결과에도 있는지 counting
            for acmg_rule, acmg_rule_counts in rules_dic_2.items():
                for rule_my in my_evidences:
                    if rule_my.startswith(acmg_rule):
                        acmg_rule_counts[1] += 1

                        for rule_3b in evidences:
                            if rule_3b.startswith(acmg_rule) and (
                                not rule_3b.endswith("-")
                            ):
                                acmg_rule_counts[0] += 1

    for acmg_rule, acmg_rule_counts in rules_dic_2.items():
        try:
            print(
                f"{acmg_rule}: {acmg_rule_counts[0]/acmg_rule_counts[1]:.4f}, {acmg_rule_counts}"
            )
        except ZeroDivisionError:
            print(f"{acmg_rule}: None, {acmg_rule_counts}")


def main():
    """_summary_
    Note:
        내가 구현한 ACMG rule 할당과, EVIDENCE의 rule 할당을 비교하는 코드.

    Print:
        1. My rule / EVIDENCE
        2. EVIDENCE / My rule
    """

    evidence_result_dic = parse_evidence(EVIDENCE_3B)
    my_evidence_result_dic = parse_my_evidence(EVIDENCE_JM)

    compare_rules(evidence_result_dic, my_evidence_result_dic)


if __name__ == "__main__":

    main()
