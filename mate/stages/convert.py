from os import path, listdir
from mate.io.message import Message as Msg
from mate.io.file_operate import FastaIO


# for running mafft, we need collect all cds with same id in different cds file to one file,
# and the file is named to the cds id
# for example:
# cds files in cds_dir
# cds_file1.cds    cds_file2.cds    cds_file3.cds    cds_file4.cds    ...
# >gid1            >gid1            >gid1            >gid1
# ATCG...          ATGG...          ACCT...          ACCG...
# output split file
# gid.fasta
# >cds_file1
# ATCG...
# >cds_file2
# ATGG...
# >cds_file3
# ACCT...
# >cds_file4
# ACCG...
def convert_cds_files_for_mafft(cds_dir, split_dir):
    Msg.info("Loading cds files")
    cds_db = {}
    for fn in listdir(cds_dir):
        Msg.info("\tLoading %s" % fn)
        sample_id = fn.replace('.cds', '')
        cds_file = path.join(cds_dir, fn)
        fasta_io = FastaIO(cds_file)
        fasta_io.read_fasta()
        for gid in fasta_io.fasta_db:
            if gid not in cds_db:
                cds_db[gid] = {}
            cds_db[gid][sample_id] = fasta_io.fasta_db[gid]
    Msg.info("Loading finished")

    Msg.info("Writing split fasta files")
    for gid in cds_db:
        split_fn = path.join(split_dir, gid+'.fa')
        Msg.info("\tWriting %s" % gid)
        with open(split_fn, 'w') as fout:
            for sample_id in cds_db[gid]:
                fout.write(">%s\n%s\n" % (sample_id, cds_db[gid][sample_id]))
    Msg.info("Writing Finished")
