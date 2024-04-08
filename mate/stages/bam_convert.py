import pysam
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool
from os import path, listdir
from mate.io.file_operate import FastaIO


class BAM2CDS:
    def __init__(self, in_bam_dir, ref_cds, out_cds_dir, thread):
        self.__in_bam_dir = in_bam_dir
        self.__ref_cds = ref_cds
        self.__out_cds_dir = out_cds_dir
        self.__thread = thread
    
    @staticmethod
    def __bam2cds(in_bam, in_ref, out_cds):
        Msg.info("\tLoading reference")
        fasta_io = FastaIO(in_ref)
        fasta_io.read_fasta()
    
        Msg.info("\tLoading bam")
        failed_flags = 4 | 256 | 512 | 1024 | 2048
        '''
        FLAG:
        1 0x1 template having multiple segments in sequencing
        2 0x2 each segment properly aligned according to the aligner
        4 0x4 segment unmapped
        8 0x8 next segment in the template unmapped
        16 0x10 SEQ being reverse complemented
        32 0x20 SEQ of the next segment in the template being reverse complemented
        64 0x40 the first segment in the template
        128 0x80 the last segment in the template
        256 0x100 secondary alignment
        512 0x200 not passing filters, such as platform/vendor quality controls
        1024 0x400 PCR or optical duplicate
        2048 0x800 supplementary alignment
        '''
    
        match_db = {gid: [{} for _ in range(len(fasta_io.fasta_db[gid]))] for gid in fasta_io.fasta_db}
    
        '''
        CIGAR:
        Op BAM Description Consumes_query Consumes_refernece
        M 0 alignment match (can be a sequence match or mismatch) yes yes
        I 1 insertion to the reference yes no
        D 2 deletion from the reference no yes
        N 3 skipped region from the reference no yes
        S 4 soft clipping (clipped sequences present in SEQ) yes no
        H 5 hard clipping (clipped sequences NOT present in SEQ) no no
        P 6 padding (silent deletion from padded reference) no no
        = 7 sequence match yes yes
        X 8 sequence mismatch yes yes
        '''
        qry_move_set = set([1, 4])
        ref_move_set = set([2, 3])
        seq_length = 0
        with pysam.AlignmentFile(in_bam, 'rb') as fin:
            for line in fin:
                if not (line.flag & failed_flags):
                    ref_name = line.reference_name
                    ref_start = line.reference_start
                    qry_start = line.query_alignment_start
                    qry_seq = line.query_sequence
                    ref_offset = ref_start
                    qry_offset = 0
                    last_ref_offset = ref_offset
                    qry_base = []
                    seq_length += line.query_alignment_length
                    for aln_type, base_cnt in line.cigar:
                        for _ in range(base_cnt):
                            if aln_type in ref_move_set:
                                # qry_base.append('-')
                                qry_sub_seq = ''.join(qry_base)
                                if qry_sub_seq not in match_db[ref_name][last_ref_offset]:
                                    match_db[ref_name][last_ref_offset][qry_sub_seq] = 0
                                match_db[ref_name][last_ref_offset][qry_sub_seq] += 1
                                qry_base = []
                                ref_offset += 1
                                last_ref_offset = ref_offset
                                if ref_offset >= len(fasta_io.fasta_db[ref_name]):
                                    break
                            elif aln_type in qry_move_set:
                                if last_ref_offset != 0:
                                    qry_base.append(qry_seq[qry_offset])
                                qry_offset += 1
                            elif aln_type == 0:
                                qry_base.append(qry_seq[qry_offset])
                                qry_sub_seq = ''.join(qry_base)
                                if qry_sub_seq not in match_db[ref_name][last_ref_offset]:
                                    match_db[ref_name][last_ref_offset][qry_sub_seq] = 0
                                match_db[ref_name][last_ref_offset][qry_sub_seq] += 1
                                qry_base = []
                                ref_offset += 1
                                qry_offset += 1
                                last_ref_offset = ref_offset
    
        ref_length = 0
        for ref_name in fasta_io.fasta_db:
            ref_length += len(fasta_io.fasta_db[ref_name])
    
        seq_cov = seq_length * 1. / ref_length
        
        Msg.info("\tWriting cds")
        with open(out_cds, 'w') as fout:
            for ref_name in sorted(match_db):
                seq = []
                for i in range(len(match_db[ref_name])):
                    info = match_db[ref_name][i]
                    if len(info) == 0:
                        best_base = ''
                    else:
                        best_base = sorted(info, key=lambda x: info[x], reverse=True)[0]
                        if info[best_base] < seq_cov / 5.:
                            best_base = fasta_io.fasta_db[ref_name][i]
                    seq.append(best_base)
                seq = ''.join(seq)
                if len(seq) != 0:
                    fout.write(">%s\n" % ref_name)
                    fout.write("%s\n" % seq)
        Msg.info("\tBam converted")
        
    def convert(self):
        Msg.info("Converting")
        pool = Pool(processes=self.__thread)
        
        res = []
        for fn in listdir(self.__in_bam_dir):
            if not fn.endswith(".bam"):
                continue
            sample_id = '.'.join(fn.split('.')[:-1])
            bam_file = path.join(self.__in_bam_dir, fn)
            out_cds = path.join(self.__out_cds_dir, "%s.cds" % sample_id)
            res.append([fn, pool.apply_async(self.__bam2cds, (bam_file, self.__ref_cds, out_cds, ))])
        pool.close()
        pool.join()

        # If subprocess failed, the error will be caught.
        for bam_fn, r in res:
            try:
                r.get()
            except Exception as e:
                Msg.warn("\tException caught with {}: {}".format(bam_fn, e))

        Msg.info("Extract finished")
