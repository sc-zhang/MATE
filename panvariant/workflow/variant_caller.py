from panvariant.io.file_operate import FastaIO
from panvariant.io.message import Message as Msg
from panvariant.operate.consensus_seq import get_consensus_seq
from pathos.multiprocessing import Pool


def variant_caller_for_single_file(aln_file):
    Msg.info("\tLoading %s" % aln_file)
    fasta_io = FastaIO(aln_file)
    fasta_io.read_aln()

    Msg.info("\tGenerating consensus sequence")
    consensus_seq = get_consensus_seq(fasta_io.fasta_db)
    seq_len = len(consensus_seq)

    full_info = []
    for i in range(seq_len):
        info = []
        ref = consensus_seq[i].upper()
        alt = {}
        for smp in sorted(fasta_io.fasta_db):
            base = fasta_io.fasta_db[smp][i]
            cur_type = 0
            if base != ref:
                if base not in alt:
                    alt[base] = len(alt)+1
                cur_type = alt[base]
            info.append(cur_type)
        full_info.append([i, ref, alt, info])

    print(full_info)
