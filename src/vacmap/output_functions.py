import pysam
import sys
def parse_sam_line(sam_line, header):
    """
    Parse a SAM alignment line (str) into a pysam.AlignedSegment object.
    This function is a pure Python replacement for pysam.AlignedSegment.fromstring.

    Args:
        sam_line (str): One alignment line in SAM format.
        header (dict): pysam header dictionary.

    Returns:
        pysam.AlignedSegment
    """
    fields = sam_line.split('\t')
    if len(fields) < 11:
        raise ValueError(f"Invalid SAM line (too few fields): {sam_line}")

    # Core fields
    qname = fields[0]
    flag = int(fields[1])
    rname = fields[2]
    pos = int(fields[3]) - 1  # pysam uses 0-based positions
    mapq = int(fields[4])
    cigar = fields[5]
    rnext = fields[6]
    pnext = int(fields[7]) - 1 if fields[7] != "0" else -1
    tlen = int(fields[8])
    seq = fields[9]
    qual = fields[10]

    seg = pysam.AlignedSegment(header)
    seg.query_name = qname
    seg.flag = flag
    seg.reference_id = seg.header.get_tid(rname) if rname != '*' else -1
    seg.reference_start = pos if rname != '*' else -1
    seg.mapping_quality = mapq
    seg.cigarstring = cigar if cigar != '*' else None
    seg.next_reference_id = seg.header.get_tid(rnext) if rnext != '*' else -1
    seg.next_reference_start = pnext if rnext != '*' else -1
    seg.template_length = tlen
    seg.query_sequence = seq if seq != '*' else None
    seg.query_qualities = pysam.qualitystring_to_array(qual) if qual != '*' else None

    # Optional fields
    for opt in fields[11:]:
        if ":" in opt:
            tag, type_, value = opt.split(":", 2)
            # Try to convert value to int/float if possible
            if type_ == 'i':
                value = int(value)
            elif type_ == 'f':
                value = float(value)
            seg.set_tag(tag, value, value_type=type_)

    return seg

def sam_bam_writer(cooked_queue, header, out_path="-"):
    """
    Write alignments to SAM/BAM file or stdout.
    Args:
        cooked_queue: multiprocessing queue or similar, yielding lists of alignments (SAM string).
        header: SAM/BAM header dict or pysam.AlignmentHeader.
        out_path: Output file path ('.sam', '.bam'), or '-' for stdout.
    """
    # Convert dict header to pysam.AlignmentHeader if needed
    if isinstance(header, dict):
        header = pysam.AlignmentHeader.from_dict(header)

    if out_path == "-":
        out_stream = sys.stdout
        # Write SAM header
        for line in str(header).split("\n"):
            if line.strip():
                out_stream.write(line + "\n")
        out_stream.flush()
        # Stream alignments
        while True:
            a_list = cooked_queue.get()
            if isinstance(a_list, int):
                break
            out_stream.write("\n".join(a_list) + "\n")
        out_stream.flush()
    elif out_path.endswith(".sam"):
        # SAM: write to file using plain python
        with open(out_path, "w") as outf:
            # Write SAM header
            for line in str(header).split("\n"):
                if line.strip():
                    outf.write(line + "\n")
            while True:
                a_list = cooked_queue.get()
                if isinstance(a_list, int):
                    break
                for a in a_list:
                    outf.write(a + "\n")
    elif out_path.endswith(".bam"):
        # BAM: use pysam
        with pysam.AlignmentFile(out_path, "wb", header=header, threads=16) as outf:
            while True:
                a_list = cooked_queue.get()
                if isinstance(a_list, int):
                    break
                for a in a_list:
                    try:
                        aligned = parse_sam_line(a, header)
                        outf.write(aligned)
                    except Exception as e:
                        print(f"Error writing alignment: {e}\nSAM line: {a}", file=sys.stderr)
    else:
        raise ValueError("Output path must end with .sam, .bam, or be '-' for stdout.")
def sam_bam_writer_asm(cooked_queue, header, out_path="-"):

    """
    Write alignments to SAM/BAM file or stdout.
    Args:
        cooked_queue: multiprocessing queue or similar, yielding lists of alignments (SAM string).
        header: SAM/BAM header dict or pysam.AlignmentHeader.
        out_path: Output file path ('.sam', '.bam'), or '-' for stdout.
    """
    # Convert dict header to pysam.AlignmentHeader if needed
    if isinstance(header, dict):
        header = pysam.AlignmentHeader.from_dict(header)

    if out_path == "-":
        out_stream = sys.stdout
        for line in str(header).split("\n"):
            if line.strip():
                out_stream.write(line + "\n")
        out_stream.flush()
        while True:
            a = cooked_queue.get()
            if isinstance(a, int):
                break
            out_stream.write(a + "\n")
        out_stream.flush()

            
    elif out_path.endswith(".sam"):
        # SAM: write to file using plain python
        with open(out_path, "w") as outf:
            # Write SAM header
            for line in str(header).split("\n"):
                if line.strip():
                    outf.write(line + "\n")
            while True:
                a = cooked_queue.get()
                if isinstance(a, int):
                    break
                outf.write(a + "\n")
                
    elif out_path.endswith(".bam"):
        # BAM: use pysam
        with pysam.AlignmentFile(out_path, "wb", header=header, threads=16) as outf:
            while True:
                a = cooked_queue.get()
                if isinstance(a, int):
                    break
                try:
                    aligned = parse_sam_line(a, header)
                    outf.write(aligned)
                except Exception as e:
                    print(f"Error writing alignment: {e}\nSAM line: {a}", file=sys.stderr)
    else:
        raise ValueError("Output path must end with .sam, .bam, or be '-' for stdout.")
        
import sys
import os
import pysam
import subprocess

def sam_bam_writer(cooked_queue, header, out_path="-"):
    """
    Write alignments to SAM/BAM file or stdout using samtools.
    If 'sorted.bam' is requested, it sorts and then indexes the output.
    """
    # Ensure header is a pysam object so we can easily stringify it
    if isinstance(header, dict):
        header = pysam.AlignmentHeader.from_dict(header)

    # --- CASE 1: Write SAM (Text) or Standard Output ---
    if out_path == "-" or out_path.endswith(".sam"):
        out_stream = sys.stdout if out_path == "-" else open(out_path, "w")
        try:
            out_stream.write(str(header))
            while True:
                a_list = cooked_queue.get()
                if isinstance(a_list, int):
                    break
                out_stream.write("\n".join(a_list) + "\n")
        finally:
            if out_path != "-":
                out_stream.close()

    # --- CASE 2: Write BAM (Using samtools) ---
    elif out_path.endswith(".bam"):
        is_sorted = out_path.endswith("sorted.bam")
        
        # Determine command: Sort or View
        if is_sorted:
            # Stream SAM -> samtools sort -> BAM File
            cmd = ["samtools", "sort", "-@", "8", '--write-index', "-o", out_path, "-"]
        else:
            # Stream SAM -> samtools view -> BAM File
            cmd = ["samtools", "view", "-b", "-@", "8", "-o", out_path, "-"]
        
        # Open subprocess with a pipe for stdin
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, encoding='utf-8', bufsize=64*1024)
        
        try:
            # 1. Write Header
            process.stdin.write(str(header))
            
            # 2. Write Alignments
            while True:
                a_list = cooked_queue.get()
                if isinstance(a_list, int):
                    break
                process.stdin.write("\n".join(a_list) + "\n")
                
            process.stdin.close()
            retcode = process.wait()
            
            if retcode != 0:
                print(f"Error: samtools exited with code {retcode}", file=sys.stderr)


        except BrokenPipeError:
            print("Error: Broken pipe while writing to samtools.", file=sys.stderr)
        except Exception as e:
            print(f"Error writing to samtools: {e}", file=sys.stderr)
            if process.poll() is None:
                process.kill()
    else:
        raise ValueError("Output path must end with .sam, .bam, sorted.bam or be '-' for stdout.")


def sam_bam_writer_asm(cooked_queue, header, out_path="-"):
    """
    Write alignments to SAM/BAM file or stdout (Assembly Mode).
    """
    if isinstance(header, dict):
        header = pysam.AlignmentHeader.from_dict(header)

    # --- CASE 1: Write SAM (Text) or Standard Output ---
    if out_path == "-" or out_path.endswith(".sam"):
        out_stream = sys.stdout if out_path == "-" else open(out_path, "w")
        try:
            out_stream.write(str(header))
            while True:
                a = cooked_queue.get()
                if isinstance(a, int):
                    break
                out_stream.write(a + "\n")
        finally:
            if out_path != "-":
                out_stream.close()

    # --- CASE 2: Write BAM (Using samtools) ---
    elif out_path.endswith(".bam"):
        is_sorted = out_path.endswith("sorted.bam")

        if is_sorted:
            cmd = ["samtools", "sort", "-@", "8", '--write-index', "-o", out_path, "-"]
        else:
            cmd = ["samtools", "view", "-b", "-@", "8", "-o", out_path, "-"]
            
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, encoding='utf-8', bufsize=64*1024)
        
        try:
            process.stdin.write(str(header))
            
            while True:
                a = cooked_queue.get()
                if isinstance(a, int):
                    break
                process.stdin.write(a + "\n")
                
            process.stdin.close()
            retcode = process.wait()
            
            if retcode != 0:
                print(f"Error: samtools exited with code {retcode}", file=sys.stderr)


        except BrokenPipeError:
            print("Error: Broken pipe while writing to samtools.", file=sys.stderr)
        except Exception as e:
            print(f"Error writing to samtools: {e}", file=sys.stderr)
            if process.poll() is None:
                process.kill()
    else:
        raise ValueError("Output path must end with .sam, .bam, sorted.bam or be '-' for stdout.")