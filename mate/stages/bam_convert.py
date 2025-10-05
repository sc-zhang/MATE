import pysam
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool
from os import path, listdir
from mate.io.file_operate import BedIO, FastaIO


class BAM2CDS:
    def __init__(self, in_bam_dir, ref_file, out_cds_dir, ref_type, thread):
        self.__in_bam_dir = in_bam_dir
        self.__ref_type = ref_type
        self.__ref_file = ref_file
        self.__out_cds_dir = out_cds_dir
        self.__thread = thread

    @staticmethod
    def __reverse_seq(seq):
        base_db = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return ''.join([base_db[_] if _ in base_db else _ for _ in seq.upper()[::-1]])

    def __bam2cds(self, in_bam, in_ref, ref_type, out_cds):
        Msg.info("\tLoading reference")
        ref_io = None
        if ref_type == 'cds':
            ref_io = FastaIO(in_ref)
            ref_io.read_fasta()
        else:
            ref_io = BedIO(in_ref)
            ref_io.read_bed()

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

        if ref_type == 'cds':
            match_db = {gid: [{"-": 0} for _ in range(len(ref_io.fasta_db[gid]))] for gid in ref_io.fasta_db}
        else:
            match_db = {gid: [{"-": 0} for _ in range(ref_io.bed_db[gid][2] - ref_io.bed_db[gid][1])] for gid in
                        ref_io.bed_db}

        '''
        CIGAR:
        Op BAM Description Consumes_query Consumes_reference
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
        qry_move_set = {1, 4}
        ref_move_set = {2, 3}
        seq_length = 0
        with pysam.AlignmentFile(in_bam, 'rb') as fin:
            for line in fin:
                if (not (line.flag & failed_flags)) and line.mapq >= 1:
                    ref_name = line.reference_name
                    if ref_type == 'bed' and ref_name not in ref_io.chrs:
                        continue
                    ref_start = line.reference_start
                    qry_start = line.query_alignment_start
                    qry_seq = line.query_sequence
                    ref_offset = ref_start
                    qry_offset = 0
                    last_ref_offset = ref_offset
                    qry_base = []
                    seq_length += line.query_alignment_length
                    if ref_type == 'cds':
                        gene_name = ref_name
                    else:
                        gene_name = ref_io.get_gid(ref_name, ref_start)
                    if gene_name not in match_db:
                        continue
                    for aln_type, base_cnt in line.cigar:
                        for _ in range(base_cnt):
                            if ref_type == 'cds':
                                if ref_offset >= len(ref_io.fasta_db[gene_name]):
                                    break
                            else:
                                if ref_offset >= ref_io.bed_db[gene_name][2]:
                                    break
                            if aln_type in ref_move_set:
                                # qry_base.append('-')
                                qry_sub_seq = ''.join(qry_base)
                                if ref_type == 'cds':
                                    if qry_sub_seq not in match_db[gene_name][last_ref_offset]:
                                        match_db[gene_name][last_ref_offset][qry_sub_seq] = 0
                                    match_db[gene_name][last_ref_offset][qry_sub_seq] += 1
                                else:
                                    if qry_sub_seq not in match_db[gene_name][last_ref_offset -
                                                                              ref_io.bed_db[gene_name][1]]:
                                        match_db[gene_name][last_ref_offset -
                                                            ref_io.bed_db[gene_name][1]][qry_sub_seq] = 0
                                    match_db[gene_name][last_ref_offset -
                                                        ref_io.bed_db[gene_name][1]][qry_sub_seq] += 1
                                qry_base = []
                                ref_offset += 1
                                last_ref_offset = ref_offset
                            elif aln_type in qry_move_set:
                                # if last_ref_offset != 0:
                                qry_base.append(qry_seq[qry_offset])
                                qry_offset += 1
                            elif aln_type == 0:
                                qry_base.append(qry_seq[qry_offset])
                                qry_sub_seq = ''.join(qry_base)
                                if ref_type == 'cds':
                                    if qry_sub_seq not in match_db[gene_name][last_ref_offset]:
                                        match_db[gene_name][last_ref_offset][qry_sub_seq] = 0
                                    match_db[gene_name][last_ref_offset][qry_sub_seq] += 1
                                else:
                                    if qry_sub_seq not in match_db[gene_name][last_ref_offset -
                                                                              ref_io.bed_db[gene_name][1]]:
                                        match_db[gene_name][last_ref_offset -
                                                            ref_io.bed_db[gene_name][1]][qry_sub_seq] = 0
                                    match_db[gene_name][last_ref_offset -
                                                        ref_io.bed_db[gene_name][1]][qry_sub_seq] += 1
                                qry_base = []
                                ref_offset += 1
                                qry_offset += 1
                                last_ref_offset = ref_offset

        Msg.info("\tWriting cds")
        with open(out_cds, 'w') as fout:
            for gene_name in sorted(match_db):
                seq = []
                for i in range(len(match_db[gene_name])):
                    info = match_db[gene_name][i]
                    if len(info) == 1:
                        best_base = sorted(info)[0]
                    else:
                        sorted_info = sorted(info, key=lambda x: info[x], reverse=True)
                        # best_base = sorted(info, key=lambda x: info[x], reverse=True)[0]
                        if info[sorted_info[1]] > 0.5 * info[sorted_info[0]]:
                            best_base = '-'
                        else:
                            best_base = sorted_info[0]
                    seq.append(best_base)
                seq = ''.join(seq)

                if ref_type == 'bed' and ref_io.bed_db[gene_name][3] == '-':
                    seq = self.__reverse_seq(seq)
                if len(seq) != 0:
                    fout.write(">%s\n" % gene_name)
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
            self.__bam2cds(bam_file, self.__ref_file, self.__ref_type, out_cds)
            # res.append([fn, pool.apply_async(self.__bam2cds, (bam_file, self.__ref_file, self.__ref_type, out_cds,))])
        pool.close()
        pool.join()

        # If subprocess failed, the error will be caught.
        for bam_fn, r in res:
            try:
                r.get()
            except Exception as e:
                Msg.warn("\tException caught with {}: {}".format(bam_fn, e))

        Msg.info("Extract finished")
