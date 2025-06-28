from enum import Enum


class CDS_Type(Enum):
    VALID = 1
    LENGTH_ERROR = 2
    MISSING_START_CODON = 3
    MISSING_STOP_CODON = 4
    EARLY_STOP_CODON = 5


class CDS_Check:
    @staticmethod
    def check_cds(cds_seq):
        start_codon = {"ATG"}
        stop_codon = {"TAG", "TAA", "TGA"}
        cds_db = {}

        cds_type = CDS_Type.VALID
        if len(cds_seq) % 3 != 0:
            cds_type = CDS_Type.LENGTH_ERROR
        else:
            if cds_seq[:3] not in start_codon:
                cds_type = CDS_Type.MISSING_START_CODON
            elif cds_seq[-3:] not in stop_codon:
                cds_type = CDS_Type.MISSING_STOP_CODON
            else:
                for _ in range(3, len(cds_seq) - 3, 3):
                    if cds_seq[_: _ + 3] in stop_codon:
                        cds_type = CDS_Type.EARLY_STOP_CODON

        return cds_type
