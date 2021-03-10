import sys
from io import BytesIO
from typing import Dict, Optional, Tuple, List
import zlib
import os
import struct
from collections import defaultdict, namedtuple


def _read_int(x: BytesIO, n_bytes: int, signed: bool=True) -> int:
    return int.from_bytes(x.read(n_bytes), byteorder="little", signed=signed)

def read_uint8(x: BytesIO) -> int:
    return _read_int(x, 1, False)

def read_int8(x: BytesIO) -> int:
    return _read_int(x, 1, True)

def read_uint16(x: BytesIO) -> int:
    return _read_int(x, 2, False)

def read_int16(x: BytesIO) -> int:
    return _read_int(x, 2, True)

def read_uint32(x: BytesIO) -> int:
    return _read_int(x, 4, False)

def read_int32(x: BytesIO) -> int:
    return _read_int(x, 4, True)

def read_int64(x: BytesIO) -> int:
    return _read_int(x, 8, True)

def read_uint64(x: BytesIO) -> int:
    return _read_int(x, 8, False)


def read_string(x: BytesIO, n_bytes: Optional[int]=None):
    if n_bytes:
        return "".join(chr(_) for _ in x.read(n_bytes))
    else:
        # Assume NUL-terminated string
        s = ""
        while True:
            c = read_int8(x)
            if c == 0:
                break
            s += chr(c)
        return s

BASE_CODES="=ACMGRSVTWYHKDBN"
def read_int4(x: BytesIO) -> Tuple[int, int]:
    z = read_int8(x)
    return ((z >> 4) & 15, z & 15) 


def read_float(x: BytesIO) -> float:
    return struct.unpack("f", x.read(4))[0]


def read_seq(x: BytesIO, l_seq: int) -> str:
    return "".join([BASE_CODES[v] for v in read_int4(x) for _ in range(int((l_seq + 1)/2))])

def read_qual(x: BytesIO, l_seq: int) -> List[int]: # TODO: Format properly, ref spec
    qual_raw = [read_uint8(x) for _ in range(l_seq)]
    return qual_raw


def read_header(data: bytes) -> Tuple[str, Dict[int, str]]:
    total_size = len(data)
    block = BytesIO(data)
    magic = read_string(block, 4)
    assert magic == "BAM\1"
    l_text = read_uint32(block)
    header_text = read_string(block, l_text)
    n_ref = read_uint32(block)
    ref_ids = {-1: "Unmapped"}
    for i in range(n_ref):
        l_name = read_uint32(block)
        name = read_string(block, l_name)
        l_ref = read_uint32(block)
        ref_ids[i] = name


    remaining_size = total_size - block.tell()
    assert remaining_size == 0
    return header_text, ref_ids

Alignment = namedtuple("Alignment", ["ref", "refId", "read_name", "mapq", "bin"] )
def read_alignments(data: bytes, ref_ids: Optional[Dict[int, str]]=None) -> List[Alignment]:
    alignments = []
    total_size = len(data)
    block = BytesIO(data)
    remaining_size = total_size
    while remaining_size > 0:
        block_start = block.tell()
        block_length = read_uint32(block)
        refId = read_int32(block)
        if ref_ids:
            ref = ref_ids[refId]
        else:
            ref = "N/A"
        
        pos = read_int32(block)
        
        l_read_name = read_uint8(block)
        mapq = read_uint8(block)
        bin = read_uint16(block)
        n_cigar_op = read_uint16(block)
        flag = read_uint16(block)
        l_seq = read_uint32(block)
        next_refID = read_int32(block)
        next_pos = read_int32(block)
        tlen = read_int32(block)
        read_name = read_string(block, l_read_name)
        cigar = [read_uint32(block) for _ in range(n_cigar_op)]

        al = Alignment(ref, refId, read_name, mapq, bin)
        alignments.append(al)

        # Skip the rest of the block, as the rest is complex to get right
        # See section 4.2.4 of the SAM v1 spec for details
        block.seek(block_start + block_length + 4)
        remaining_size = total_size - block.tell()
        continue

        # TODO: The following is incomplete
        seq = read_seq(block, l_seq)
        qual = read_qual(block, l_seq)
        tag = read_string(block, 2)

        val_type = read_string(block, 1)
        val_type_switch = {
            "A": lambda: read_string(block, 1),
            "c": lambda: read_int8(block),
            "C": lambda: read_uint8(block),
            "s": lambda: read_int16(block),
            "S": lambda: read_uint16(block),
            "i": lambda: read_int32(block),
            "I": lambda: read_uint32(block),
            "f": lambda: read_float(block),
            "Z": lambda: read_string(block)
        }

        value = val_type_switch[val_type]()
    return alignments




def inspect_bam(filename: str, fast: bool=False) -> None:
    total_size = os.stat(filename).st_size
    header = None
    ref_ids = None
    current_block = -1
    with open(filename, 'rb') as bam:
        while bam.tell() < total_size:
            block_start = bam.tell()
            current_block += 1

            if fast:
                # Skip reading unnecessary fields
                bam.seek(block_start + 10)
                xlen = read_uint16(bam)
                bam.seek(block_start + 16)
                bsize = read_uint16(bam)
                cdata_len = bsize - xlen - 19
                bam.seek(block_start + cdata_len + 22)
                cdata = None
                isize = read_uint32(bam)
            else:
                # Read header fields
                id1 = read_uint8(bam)
                assert id1 == 31
                id2 = read_uint8(bam)
                assert id2 == 139
                cm = read_uint8(bam)
                assert cm == 8
                flg = read_uint8(bam)
                assert flg == 4
                mtime = read_uint32(bam)
                xfl = read_uint8(bam)
                _os = read_uint8(bam)
                xlen = read_uint16(bam)
                si1 = read_uint8(bam)
                assert si1 == 66
                si2 = read_uint8(bam)
                assert si2 == 67
                slen = read_uint16(bam)
                bsize = read_uint16(bam)
                cdata_len = bsize - xlen - 19
                # Read compressed data
                cdata = bam.read(cdata_len)
                # Read checksum and expected size of uncompressed data
                crc32 = read_uint32(bam) # TODO: Check checksum
                isize = read_uint32(bam)

            if cdata is not None:
                # Decompress data
                data = zlib.decompress(cdata, -15)

                # First block is header, subsequent blocks are alignments
                if block_start == 0:
                    header, ref_ids = read_header(data)
                else:
                    alignments = read_alignments(data, ref_ids)
                
                assert len(data) == isize, f"Uncompressed data does not match expected size ({len(data)} != {isize})"

            block_end = bam.tell()

            assert block_end - block_start == bsize +1

            print(f"Block #{current_block} -- Block start: {block_start} Block end: {block_end} raw size: {cdata_len} data size: {isize}" )

def make_virtual_offset(coffset, uoffset):
    return (coffset << 16) | uoffset
    pass

def split_virtual_offset(virtual_offset):
    coffset = virtual_offset >> 16
    uoffset = virtual_offset ^ (coffset << 16)
    return coffset, uoffset

def inspect_bai(filename):
    total_size = os.stat(filename).st_size
    with open(filename, "rb") as bai:
        magic = read_int32(bai)
        n_ref = read_uint32(bai)
        for i in range(n_ref):
            n_bins = read_uint32(bai) 
            #print("n_bins: ", n_bins)
            for j in range(n_bins):
                bin = read_uint32(bai)
                n_chunks = read_uint32(bai)
                for k in range(n_chunks):
                    chunk_beg = read_uint64(bai)
                    chunk_end = read_uint64(bai)
                    assert make_virtual_offset(*split_virtual_offset(chunk_beg)) == chunk_beg
                    assert make_virtual_offset(*split_virtual_offset(chunk_end)) == chunk_end
                    # coffset<<16|uoffset
                    print(f"ref_id: {i} bin: {bin} chunk_beg: {split_virtual_offset(chunk_beg)} chunk_end: {split_virtual_offset(chunk_end)}")
            n_intv = read_uint32(bai)
            for k in range(n_intv):
                ioffset = read_uint64(bai)
        if bai.tell() != total_size:
            n_no_coor = read_uint64(bai)
        else:
            n_no_coor = "N/A"
        print(f"Num unplaced unmapped reads: {n_no_coor}")
        assert bai.tell() == total_size



if __name__ == "__main__":
    f = sys.argv[1]
    if f.endswith(".bam"):
        if len(sys.argv) > 2:
            fast = sys.argv[2] == "--fast"
        else:
            fast = False
        inspect_bam(f, fast)
    elif f.endswith(".bai"):
        inspect_bai(f)
    else:
        raise ValueError(f"Unknown filetype for {f}")
