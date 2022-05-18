# Module for PM4, BP3
# repetitive region, in-frame deletion/insertion, stop-loss


def is_inframe_indel_change_var(var_infos_dic: dict, df_col2idx: dict) -> bool:
    """_summary_
    Note:
        해당 변이의 결과가 inframe deletion/insertion 인지 확인하는 함수.

    Args:
        var_infos_dic (dict): variation-feature의 세부 정보 딕셔너리.
        df_col2idx (dict): column index dictionary

    Returns:
        bool: True or False
    """

    var_infos = var_infos_dic["var_infos"]  # uploaded_id

    if "inframe" in var_infos[df_col2idx["consequence"]]:
        return True
    else:
        return False


def is_stop_lost_var(var_infos_dic: dict, df_col2idx: dict) -> bool:
    """_summary_
    Note:
        해당 변이의 결과가 stop_lost 인지 확인하는 함수.

    Args:
        var_infos_dic (dict): variation-feature의 세부 정보 딕셔너리.
        df_col2idx (dict): column index dictionary

    Returns:
        bool: True or False
    """

    var_infos = var_infos_dic["var_infos"]  # uploaded_id

    if "stop_lost" in var_infos[df_col2idx["consequence"]]:
        return True
    else:
        return False


def match_with_repeat_region(
    chrom: str, var_start: int, var_end: int, repeat_db_dic: dict
) -> tuple:
    """_summary_
    Note:
        해당 inframe_indel이 repetitive region에 포함되는지 확인하는 함수.
        variation의 '절반이상'이 반복서열에 포함되는 경우를 기준으로 하였다.

        repeat_db_dic = {
                    "1": [(10001, 10467), (10469, 11447), ...]
                    "2": [(444214, 444451), (445112, 445234), ...]
                } (start, end)은 오름차순 정렬되어 있음.

        repeat_db_dic의 구간은 오름차순으로 정렬되어 있으므로, 원하는 구간을 빠르게 찾기위해 binary search를
        이용할 수 있다. '절반이상'이 반복서열 구간에 포함되기 위해서는 변이의 중간위치(center_pos)가 반복서열
        구간 내에 반드시 포함되어야 한다.

                    (  *  ) <- Varirant
        ...|----|....|--------|.........|----|... <- repeat regions
                          (  *  )

        따라서, 반복서열의 start postion을 기준으로 해당 변이의 중간지점이 포함되는 구간을 찾고, center_pos가
        해당 반복서열의 (start, end)에 속하는지 확인한다. 반복서열에 속하면 bp3를, 아니면 pm4를 부여한다.

    Args:
        Chrom (str): 해당 변이가 속하는 염색체
        var_start (int): 변이의 시작지점
        var_end (int): 변이의 끝지점 (insertion의 경우 삽입된 위치 1개가 표시)
        repeat_db_dic (dict): {Chr: [(start, end), ..]} 로 정리된 Dictionary

    Returns:
        tuple: (pm3, bp4)
    """

    pm4, bp3 = 0, 0

    # chrom  -> [ (start,end), (10001, 10467), (10469, 11447), ...]
    repetitive_regions: list = repeat_db_dic[chrom]

    # Binary search
    first, last = 0, len(repetitive_regions) - 1
    center_pos = round((var_start + var_end) / 2)

    while first <= last:
        mid = (first + last) // 2
        if (
            repetitive_regions[mid][0]
            <= center_pos
            < repetitive_regions[mid + 1][0]
        ):
            break
        if repetitive_regions[mid][0] < var_start:
            first = mid + 1
        else:
            last = mid - 1

    # repetitive region에 절반이상 포함되면 bp3 부여
    if center_pos <= repetitive_regions[mid][1]:
        bp3 = 1
    else:
        pm4 = 1

    return (pm4, bp3)


def execute(proband_var_df: object, repeat_db_dic: dict,) -> object:
    """_summary_
    Note:
        ACMG rule 중에서, pm4/bp3 에 해당하는 룰을 구현한 모듈이다. 각각의 rule에 대한 설명은 다음과 같다.

        pm4: non-repeat region에서 inframe indel 또는 stop-lost 변이로 인한 단백질의 길이 변화
        bp3: repetitive region에서 발생한 inframe indel로 인한 단백질의 길이 변화

        본 모듈에서는 1가지의 데이터베이스가 사용되었다.
        - repeatMasker: SINE, LINE 등 여러 class의 반복서열 위치가 정리되어 있다. 가공한 정보에는 위치정보만 포함됨.
        {
            "1": [(10001, 10467), (10469, 11447), ...]
            "2": [(444214, 444451), (445112, 445234), ...]
        } (start, end)는 오름차순 정렬되어 있음.

        본 모듈의 작동방식은 다음과 같다.
        먼저, 각 variant에 대하여 inframe indel인지 stop lost인지 조사한다(is_inframe_indel_change_var(),
        is_stop_lost_var()). 만약 해당 변이가 inframe indel 이라면, 해당 구역이 repetitive region인지 확인
        한 뒤, 반복서열 위치이면 bp3를, 아니면 pm4 rule을 부여하였다(match_with_repeat_region()). 이 때, 반복
        서열 구간에 범위가 겹칠 수가 있는데, 이 경우 해당 범주에 절반이상 포함되는 것을 기준으로 하였다. 이 때, repeat
        DB는 구간이 오름차순 정렬되어 있으므로, binary search를 통해 빠르게 구간을 특정할 수 있다.
        Insertion의 경우, 구간의 표현이 (pos, pos+1)과 같이 포현되는데, 반복서열의 바로 옆 서열인 경우도 bp3를
        부여하였다. 만약 해당 변이의 전사체가 Stop_lost라면, pm4 rule 부여하였다.

    Args:
        proband_var_df (object): Class VariantDF <- VEP file
        repeat_db_dic (dict): repeatMasker 데이터베이스를 미리 가공한 dictionary

    Returns:
        object: proband_variant_DF의 variant_dic["evidnece_score_dic"]이 update된 object.

    Examples:
    >>> (repeat start, var_start, var_end, repeat end) -> rule
    (66766337, 66766356 66766357, 66766356) -> bp3
    (44520238, 44520238 44520240, 44520263) -> bp3
    (131241030, 131241029 131241030, 131241054) -> bp3
    """

    df_col2idx = proband_var_df.df_col2idx
    variant_dic = proband_var_df.variant_dic

    for var_id in variant_dic:
        for var_feature in variant_dic[var_id]:

            pm4, bp3 = 0, 0

            if is_inframe_indel_change_var(
                variant_dic[var_id][var_feature], df_col2idx
            ):
                # repetitive region을 확인하기 위한 정보 할당.
                chrom = variant_dic[var_id][var_feature]["var_infos"][
                    df_col2idx["chrom"]
                ]
                location: list = variant_dic[var_id][var_feature]["var_infos"][
                    df_col2idx["location"]
                ]
                # (66766356, 66766357)
                var_start = location[0]
                var_end = location[-1]

                # repeat region 검색 및 Rule 할당.
                pm4, bp3 = match_with_repeat_region(
                    chrom, var_start, var_end, repeat_db_dic
                )
            elif is_stop_lost_var(variant_dic[var_id][var_feature], df_col2idx):
                pm4 = 1
            else:
                continue

            variant_dic[var_id][var_feature]["evidence_score_dic"]["pm4"] = pm4
            variant_dic[var_id][var_feature]["evidence_score_dic"]["bp3"] = bp3

    proband_var_df.variant_dic = variant_dic

    return proband_var_df
