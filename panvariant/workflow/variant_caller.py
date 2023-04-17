from panvariant.io.file_operate import FastaIO, VariantIO
from panvariant.io.message import Message as Msg
from panvariant.operate.consensus_seq import get_consensus_seq
from pathos.multiprocessing import Pool
from os import listdir, path


def __variant_caller_for_single_file(aln_file, var_file):
    Msg.info("\tLoading %s" % aln_file)
    fasta_io = FastaIO(aln_file)
    fasta_io.read_aln()

    Msg.info("\tGenerating consensus sequence")
    consensus_seq = get_consensus_seq(fasta_io.fasta_db)
    seq_len = len(consensus_seq)

    Msg.info("\tDropping low quality sequences")
    # drop samples with missing rate larger than 70%
    aln_db = {}
    for smp in fasta_io.fasta_db:
        cnt = 0
        for i in range(seq_len):
            if fasta_io.fasta_db[smp][i] == '-':
                cnt += 1
        if cnt <= seq_len*.7:
            aln_db[smp] = fasta_io.fasta_db[smp]

    # if only single sample retained, return.
    retain_sample_cnt = len(aln_db)
    if retain_sample_cnt <= 1:
        Msg.info("\tToo few samples, Abort")
        return

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
    Msg.info("\tFinished")


def variant_caller(aln_dir, var_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Variant calling")
    for fn in listdir(aln_dir):
        Msg.info("\tCalling %s" % fn)
        aln_file = path.join(aln_dir, fn)
        var_file = path.join(var_dir, fn.replace('.aln', '.var'))
        pool.apply_async(__variant_caller_for_single_file, (aln_file, var_file, ))
    pool.close()
    pool.join()
    Msg.info("Variant called")
