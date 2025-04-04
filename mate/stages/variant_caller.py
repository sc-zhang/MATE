from mate.io.file_operate import FastaIO, VariantIO, AlignIO
from mate.io.message import Message as Msg
from mate.base.consensus_seq import get_consensus_seq
from pathos.multiprocessing import Pool
from os import listdir, path


def __variant_caller_for_single_file(aln_file, var_file, cleanup_aln_file, variant_filter):
    Msg.info("\tLoading %s" % aln_file)
    fasta_io = FastaIO(aln_file)
    fasta_io.read_aln()

    Msg.info("\tGenerating consensus sequence")
    consensus_seq = get_consensus_seq(fasta_io.fasta_db)
    seq_len = len(consensus_seq)

    Msg.info("\tDropping low quality sequences")

    # drop samples with low support kmers
    support_db = {}
    for smp in fasta_io.fasta_db:
        support_db[smp] = []

    seq_cnt = len(fasta_io.fasta_db)
    kmer_length, lower_threshold, missing_threshold = variant_filter.split(':')
    kmer_length = int(kmer_length)
    lower_threshold = float(lower_threshold)

    # remove base if more than lower_threshold ratio of samples with lower supported kmer at this position
    remove_pos = set()
    for i in range(seq_len-kmer_length+1):
        cnt_db = {}
        for smp in fasta_io.fasta_db:
            kmer = fasta_io.fasta_db[smp][i: i+kmer_length]
            if kmer not in cnt_db:
                cnt_db[kmer] = 0
            cnt_db[kmer] += 1
        for kmer in cnt_db:
            if cnt_db[kmer]*1./seq_cnt < lower_threshold:
                remove_pos.add(i)

    # remove base if more than missing_threshold ratio of samples with '-'
    missing_threshold = float(missing_threshold) * seq_cnt
    for pos in range(seq_len):
        cnt = 0
        for smp in fasta_io.fasta_db:
            if fasta_io.fasta_db[smp][pos] != '-':
                cnt += 1
        if cnt >= missing_threshold:
            remove_pos.add(pos)

    cleanup_aln_db = {}
    aln_db = {}
    for smp in fasta_io.fasta_db:
        cleanup_aln_db[smp] = []
        for pos in range(seq_len):
            if pos not in remove_pos:
                cleanup_aln_db[smp].append(fasta_io.fasta_db[smp][pos])
            else:
                cleanup_aln_db[smp].append('-')
    for smp in cleanup_aln_db:
        aln_db[smp] = ''.join(cleanup_aln_db[smp])

    Msg.info("\tRegenerating consensus sequence")
    consensus_seq = get_consensus_seq(aln_db)
    seq_len = len(consensus_seq)

    Msg.info("\tChecking each site")
    full_info = []
    for i in range(seq_len):
        info = []
        ref = consensus_seq[i].upper()
        alt = {}
        for smp in sorted(aln_db):
            base = aln_db[smp][i]
            cur_type = 0
            if base != ref:
                if base not in alt:
                    alt[base] = len(alt)+1
                cur_type = alt[base]
            info.append(cur_type)
        if alt:
            full_info.append([i, ref, alt, info])

    Msg.info("\tWriting results")
    var_io = VariantIO()
    var_io.write_file(var_file, sorted(aln_db.keys()), full_info)

    aln_io = AlignIO()
    aln_io.write_file(cleanup_aln_file, aln_db)
    Msg.info("\tFinished")


def variant_caller(aln_dir, var_dir, cleanup_aln_dir, variant_filter, thread):
    pool = Pool(processes=thread)
    Msg.info("Variant calling")

    res = []
    for fn in listdir(aln_dir):
        Msg.info("\tCalling %s" % fn)
        aln_file = path.join(aln_dir, fn)
        var_file = path.join(var_dir, fn.replace('.aln', '.var'))
        cleanup_aln_file = path.join(cleanup_aln_dir, fn)
        res.append([fn, pool.apply_async(__variant_caller_for_single_file,
                                         (aln_file, var_file, cleanup_aln_file, variant_filter,))])
    pool.close()
    pool.join()

    # If subprocess failed, the error will be caught.
    for aln_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}".format(aln_fn, e))

    Msg.info("Variant called")
