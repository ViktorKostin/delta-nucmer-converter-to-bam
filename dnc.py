from Bio.Seq import Seq
import pysam, argparse

def delta_to_dict(delta, ref_name):
    result = []
    delta_new = []
    for contig in delta.split('>{} '.format(ref_name))[1:]:
        delta_new.append(contig.split('\n')[:-1])
    for contig in delta_new[:-1]:
        contig_d = {}
        contig_d['chunks'] = []
        mutations = []
        for elem in contig:
            if(elem.find('Contig') > -1):
                contig_d['name'] = elem.split(' ')[0]
                contig_d['length'] = to_int(elem.split(' ')[2])
            else:
                if(elem.find(' ') > -1):
                    contig_almt = {}
                    contig_almt['ref_start'], contig_almt['ref_end'],\
                    contig_almt['contig_start'], contig_almt['contig_end'], _, _, _ = to_int(elem.split(' '))
                    contig_d['chunks'].append(contig_almt)
                    contig_almt['mutations'] = []
                else:
                    if(elem != '0'):
                        mutations.append(elem)
                    else:
                        contig_almt['mutations'] = to_int(mutations)
                        mutations = []
        result.append(contig_d)
    return result

def to_int(elem):
    if(type(elem) == str):
        return int(elem)
    if(type(elem) == list):
        return list(map(lambda x: int(x), elem))

def get_data(file):
    result = file.read()
    if(file.name.find("reference.fasta") > -1):
        ref_name = result.split('\n')[0].replace('>', '')
        result = result.replace('\n', '').replace(ref_name, '')
        return result, ref_name
    return result

def mutations_to_cigar(mutations, contig_length, chunk_size):
    stack = []
    cigar = []
    deletion = 0
    insertion = 0
    if(len(mutations) == 0):
        return [(0, chunk_size)]
    for mutation in mutations:
        if(mutation == 1):
            cigar.append((2,1))
        elif(mutation == -1):
            cigar.append((1,1))
        elif(mutation > 1):
            cigar.append((0,mutation-1))
            cigar.append((2,1))
        elif(mutation < 1):
            cigar.append((0,abs(mutation)-1))
            cigar.append((1,1))
    return cigar


def cigar_compress(cigar):
    cigar_new = []
    for c in cigar:
        if(len(cigar_new) == 0):
            cigar_new.append(c)
        elif(cigar_new[-1][0] == c[0]):
            cigar_new[-1] = (c[0], cigar_new[-1][1] + c[1])
        else:
            cigar_new.append(c)
            pass
    return cigar_new

def cigar_size(cigar):
    result = 0
    for c in cigar:
        result += c[1]
    return result

def cigar_count_dels(cigar):
    result = 0
    for c in cigar:
        if(c[0] == 2):
            result += c[1]
    return result

def delta_prepare(delta, assembly, ref_name):
    delta_d = delta_to_dict(delta, ref_name)
    for contig in delta_d:
        for chunk in contig['chunks']:
            start = min(chunk['contig_start'], chunk['contig_end'])
            end = max(chunk['contig_start'], chunk['contig_end'])
            chunk_size = abs(end - start)+1
            cigar = mutations_to_cigar(chunk['mutations'], contig['length'], chunk_size)
            cigar = cigar_compress(cigar)
            cigar_num_dels = cigar_count_dels(cigar)
            match_end = (0, chunk_size-cigar_size(cigar)+cigar_count_dels(cigar))
            if(match_end[1] > 0):
                cigar.append(match_end)
            if(chunk['contig_start'] > chunk['contig_end']):
                seq = assembly.split(contig['name']+'\n')[1].split('\n>Contig')[0]\
                            [chunk['contig_end']-1:chunk['contig_start']]
                seq = Seq(seq).reverse_complement().__str__()
            else:
                seq = assembly.split(contig['name']+'\n')[1].split('\n>Contig')[0]\
                            [chunk['contig_start']-1:chunk['contig_end']]
            chunk['cigar'] = cigar
            chunk['sequence'] = seq
    return delta_d

def delta_d_to_chunks(delta, assembly, ref_name):
    chunks = []
    for contig in delta_prepare(delta, assembly, ref_name):
        for idx,chunk in enumerate(contig['chunks']):
            chunk['name'] = contig['name']
            chunk['length'] = contig['length']
            chunks.append(chunk)
    chunks = sorted(chunks, key=lambda chunk: chunk['ref_start'])
    return chunks


def create_bam(fname, assembly, delta, reference_len, ref_name):
    header = { 'HD': {'VN': '1.0'},
                'SQ': [{'LN': reference_len, 'SN': ref_name}] }
    with pysam.AlignmentFile('./{}'.format(fname), "wb", header=header) as outf:
        for chunk in delta_d_to_chunks(delta, assembly, ref_name):
            a = pysam.AlignedSegment()
            a.query_name = chunk['name']
            a.query_sequence=chunk['sequence']
            a.flag = 16
            a.reference_id = 0
            a.reference_start = chunk['ref_start']-1
            a.cigar = chunk['cigar']
            a.mapping_quality = 255
            outf.write(a)

def main(args):
    assembly = get_data(args.assembly)
    delta = get_data(args.delta)
    reference, ref_name = get_data(args.reference)
    ref_len = len(reference)
    fname = args.fname
    create_bam(fname, assembly, delta, ref_len, ref_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', metavar='reference.fasta', type=argparse.FileType('r', encoding='UTF-8'), help='reference genome file')
    parser.add_argument('assembly', metavar='assembly.fasta', type=argparse.FileType('r', encoding='UTF-8'), help='assembly file with relative contigs to reference genome')
    parser.add_argument('delta', metavar='out.delta', type=argparse.FileType('r', encoding='UTF-8'), help='delta generated from nucmer')
    parser.add_argument('fname', metavar='out.bam', type=str, nargs='?', default='out.bam', help='filename for generated alignment')
    args = parser.parse_args()
    try:
        main(args)
    except Exception as e:
        print(e)
