```python
#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Calling of hyperedited reads, sites, and clusters
#
# Description:
#   Recalculates MD and NM tags against the original reference genome, identifies
#   substitution-class mismatches, applies hyperediting criteria, and writes BED
#   files containing hyperedited reads, individual sites, and editing clusters.
#
# Usage:
#   python3 05_calling_per_sites.py \
#     OUTDIR SAMPLE CLASS REF.fa classes.tsv
#
# Requirements:
#   Python 3
#   SAMtools 1.20
# -----------------------------------------------------------------------------

import collections
import math
import os
import subprocess
import sys


OUTDIR, SAMPLE, CLASS, REF_ORIG, CLASSES = sys.argv[1:6]

# Load the nucleotide substitution associated with the selected class.
class_map = {}

with open(CLASSES) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue

        label, source, destination = line.split()
        class_map[label] = (source, destination)

SRC, DST = class_map[CLASS]

BAM_DIR = f"{OUTDIR}/04_masked_align/{CLASS}/bam"
R1_BAM = f"{BAM_DIR}/{SAMPLE}_R1.{CLASS}.bam"
R2_BAM = f"{BAM_DIR}/{SAMPLE}_R2.{CLASS}.bam"

CALL_DIR = f"{OUTDIR}/05_calls/{CLASS}"


def sam_stream_with_md(bam):
    """Recalculate MD and NM tags against the original reference."""

    command = (
        f"samtools calmd -b {bam} {REF_ORIG} "
        f"| samtools view -h -F 0x900"
    )

    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    for line in process.stdout:
        yield line

    process.stdout.close()
    process.wait()


def parse_cigar(cigar):
    """Convert a CIGAR string into operation-length pairs."""

    numbers = ""
    operations = []

    for character in cigar:
        if character.isdigit():
            numbers += character
        else:
            operations.append((character, int(numbers)))
            numbers = ""

    return operations


def quality_scores(quality):
    """Convert Phred+33 quality characters into numerical scores."""

    return [ord(character) - 33 for character in quality.rstrip()]


def class_match_unstranded(reference_base, read_base, source, destination):
    """Match the selected substitution or its reverse complement."""

    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }

    reference_base = reference_base.upper()
    read_base = read_base.upper()

    forward_match = (
        reference_base == source
        and read_base == destination
    )

    complementary_match = (
        reference_base == complement[source]
        and read_base == complement[destination]
    )

    return forward_match or complementary_match


def expand_mismatches(cigar, md, read_sequence, position, flag):
    """
    Return mismatches as:
    reference position, read position, reference base, read base, and base quality.
    """

    strand = "-" if int(flag) & 16 else "+"

    operations = parse_cigar(cigar)

    read_position = 0
    reference_position = position
    md_position = 0
    mismatches = []

    def next_md_token(md_string, index):

        if index >= len(md_string):
            return None, index

        if md_string[index].isdigit():
            end = index

            while end < len(md_string) and md_string[end].isdigit():
                end += 1

            return ("=", int(md_string[index:end])), end

        if md_string[index] == "^":
            end = index + 1

            while end < len(md_string) and md_string[end].isalpha():
                end += 1

            return ("D", md_string[index + 1:end]), end

        return ("X", md_string[index].upper()), index + 1

    while operations:

        operation, length = operations[0]

        if length == 0:
            operations.pop(0)
            continue

        if operation in ("M", "X", "="):

            token, md_position = next_md_token(md, md_position)

            if token is None:
                read_position += length
                reference_position += length
                operations[0] = (operation, 0)
                continue

            token_type, value = token

            if token_type == "=":

                step = min(value, length)

                read_position += step
                reference_position += step
                length -= step

                operations[0] = (operation, length)

                if value > step:
                    md = str(value - step) + md[md_position:]
                    md_position = 0

            elif token_type == "D":

                reference_position += len(value)

            else:

                if length > 0 and read_position < len(read_sequence):

                    read_base = read_sequence[read_position][0]
                    base_quality = read_sequence[read_position][1]

                    mismatches.append(
                        (
                            reference_position,
                            read_position,
                            value,
                            read_base,
                            base_quality,
                        )
                    )

                    read_position += 1
                    reference_position += 1
                    length -= 1

                    operations[0] = (operation, length)

        elif operation == "I":
            read_position += length
            operations[0] = (operation, 0)

        elif operation in ("D", "N"):
            reference_position += length
            operations[0] = (operation, 0)

        elif operation in ("S", "H"):
            if operation == "S":
                read_position += length

            operations[0] = (operation, 0)

        else:
            operations[0] = (operation, 0)

    read_length = sum(
        length
        for operation, length in parse_cigar(cigar)
        if operation in ("M", "=", "X", "I", "S")
    )

    return mismatches, read_length, strand


def parse_read_end(bam):
    """
    Collect alignments by read name and retain the alignment with the highest
    substitution-class fraction when it exceeds the next-best alignment by 0.10.
    """

    alignments_by_read = collections.defaultdict(list)

    for line in sam_stream_with_md(bam):

        if line.startswith("@"):
            continue

        fields = line.rstrip("\n").split("\t")

        qname = fields[0]
        flag = int(fields[1])
        chromosome = fields[2]
        position = int(fields[3])
        mapq = int(fields[4])
        cigar = fields[5]
        sequence = fields[9]
        quality = fields[10]

        if chromosome == "*" or cigar == "*":
            continue

        md = None

        for tag in fields[11:]:
            if tag.startswith("MD:Z:"):
                md = tag[5:]
                break

        if md is None:
            continue

        read_sequence = list(zip(sequence, quality_scores(quality)))

        mismatches, read_length, strand = expand_mismatches(
            cigar,
            md,
            read_sequence,
            position - 1,
            flag,
        )

        # Retain only mismatches with Phred quality >= 30.
        mismatches = [
            mismatch
            for mismatch in mismatches
            if mismatch[4] >= 30
        ]

        total_mismatches = len(mismatches)

        edited_positions = [
            mismatch[0]
            for mismatch in mismatches
            if class_match_unstranded(
                mismatch[2],
                mismatch[3],
                SRC,
                DST,
            )
        ]

        number_of_edits = len(edited_positions)

        editing_fraction = (
            number_of_edits / total_mismatches
            if total_mismatches > 0
            else 0.0
        )

        reference_span = sum(
            length
            for operation, length in parse_cigar(cigar)
            if operation in ("M", "=", "X", "D", "N")
        )

        start = position - 1
        end = start + reference_span

        alignments_by_read[qname].append(
            {
                "rname": chromosome,
                "start": start,
                "end": end,
                "strand": strand,
                "mapq": mapq,
                "read_len": read_length,
                "mis_total": total_mismatches,
                "edit_refpos": edited_positions,
                "frac": editing_fraction,
                "flag": flag,
            }
        )

    retained_alignments = {}

    for qname, alignments in alignments_by_read.items():

        if len(alignments) > 50:
            continue

        if not alignments:
            continue

        ranked = sorted(
            alignments,
            key=lambda alignment: alignment["frac"],
            reverse=True,
        )

        if (
            len(ranked) == 1
            or ranked[0]["frac"] >= ranked[1]["frac"] + 0.10
        ):
            retained_alignments[qname] = ranked[0]

    return retained_alignments


def passes_hyperediting_filters(record):
    """Apply the hyperediting thresholds used in the analysis."""

    read_length = record["read_len"]
    edited_positions = record["edit_refpos"]
    total_mismatches = record["mis_total"]

    if read_length <= 0:
        return False

    if len(edited_positions) < max(
        1,
        math.ceil(0.05 * read_length),
    ):
        return False

    if (
        total_mismatches == 0
        or len(edited_positions) / total_mismatches < 0.60
    ):
        return False

    first_edit = min(edited_positions)
    last_edit = max(edited_positions)

    cluster_length = last_edit - first_edit + 1

    if cluster_length < 0.10 * read_length:
        return False

    left_boundary = record["start"] + 0.20 * read_length
    right_boundary = record["start"] + 0.80 * read_length

    # Reference offsets are used as an approximation of read-relative windows.
    if last_edit <= left_boundary and first_edit <= left_boundary:
        return False

    if first_edit >= right_boundary and last_edit >= right_boundary:
        return False

    return True


def paired_end_sanity(read1, read2):
    """Check chromosome, distance, and orientation of paired alignments."""

    if read1["rname"] != read2["rname"]:
        return False

    first_start = min(read1["start"], read2["start"])
    second_start = max(read1["start"], read2["start"])

    if second_start - first_start > 500_000:
        return False

    return read1["strand"] != read2["strand"]


# Parse both read ends.
read1_alignments = parse_read_end(R1_BAM)
read2_alignments = parse_read_end(R2_BAM)

read1_pass = {
    qname: record
    for qname, record in read1_alignments.items()
    if record["edit_refpos"] and passes_hyperediting_filters(record)
}

read2_pass = {
    qname: record
    for qname, record in read2_alignments.items()
    if record["edit_refpos"] and passes_hyperediting_filters(record)
}

# Join paired ends when both are available.
retained_reads = []

for qname in set(read1_pass) | set(read2_pass):

    read1 = read1_pass.get(qname)
    read2 = read2_pass.get(qname)

    if read1 and read2:

        if paired_end_sanity(read1, read2):
            retained_reads.extend(
                [
                    (qname, read1),
                    (qname, read2),
                ]
            )

    elif read1:
        retained_reads.append((qname, read1))

    elif read2:
        retained_reads.append((qname, read2))


os.makedirs(CALL_DIR, exist_ok=True)

hyperreads_bed = (
    f"{CALL_DIR}/{SAMPLE}.{CLASS}.hyperreads.bed"
)

sites_bed = (
    f"{CALL_DIR}/{SAMPLE}.{CLASS}.sites.bed"
)

clusters_bed = (
    f"{CALL_DIR}/{SAMPLE}.{CLASS}.clusters.bed"
)

unique_sites = set()

with open(hyperreads_bed, "w") as reads_output, \
     open(clusters_bed, "w") as clusters_output:

    for qname, record in retained_reads:

        reads_output.write(
            f"{record['rname']}\t"
            f"{record['start']}\t"
            f"{record['end']}\t"
            f"{qname}\t"
            f"0\t"
            f"{record['strand']}\n"
        )

        first_edit = min(record["edit_refpos"])
        last_edit = max(record["edit_refpos"])

        clusters_output.write(
            f"{record['rname']}\t"
            f"{first_edit}\t"
            f"{last_edit + 1}\n"
        )

        for position in record["edit_refpos"]:
            unique_sites.add(
                (
                    record["rname"],
                    position,
                )
            )

with open(sites_bed, "w") as sites_output:

    for chromosome, position in sorted(unique_sites):
        sites_output.write(
            f"{chromosome}\t"
            f"{position}\t"
            f"{position + 1}\n"
        )
```

