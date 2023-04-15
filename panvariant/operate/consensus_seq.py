def get_consensus_seq(aln_db):
    consensus_seq = ""
    smp_list = sorted(aln_db.keys())
    for _ in range(len(aln_db[smp_list[0]])):
        cnt_db = {}
        for smp in smp_list:
            base = aln_db[smp][_]
            if base not in cnt_db:
                cnt_db[base] = 0
            cnt_db[base] += 1
        max_base = ""
        max_cnt = 0
        for base in cnt_db:
            if cnt_db[base] > max_cnt:
                max_cnt = cnt_db[base]
                max_base = base
        if max_cnt == len(smp_list):
            consensus_seq += max_base.upper()
        else:
            consensus_seq += max_base.lower()
    return consensus_seq
