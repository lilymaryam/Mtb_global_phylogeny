import argparse
#from termios import FF1

#makes the assumption that there is no header
#this file assumes 2 groups compared to each other (p sure jaccard needs that anyway)

parser = argparse.ArgumentParser()
parser.add_argument('--tsv1', '-t1', type=str, action='store', required=True, help='path to tsv of sample neighborhoods for first dataset')
parser.add_argument('--tsv2', '-t2', type=str, action='store', required=True, help='path to tsv of sample neighborhoods for second dataset')
parser.add_argument('--output_name', '-o', type=str, action='store', default='jaccard_index_scores.tsv', help='filepath and filename to the output file of all jaccard scores')
args = parser.parse_args()

tsv1 = args.tsv1
tsv2 = args.tsv2
out = args.output_name

def fix_file(file):
    sampls = {}
    with open(file) as f:
        for line in f:
            if not line.startswith('ref'):
                line = line.strip().split()
                if line[0] not in sampls:
                    sampls[line[0]] = [line[1]]
                else:
                    sampls[line[0]].append(line[1])
                #print(line)
    print('sampls', sampls)
    return sampls


def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(set(list1)) + len(set(list2))) - intersection
    return float(intersection) / union

def check_overlap(list1, list2):
    #results = {}
    print(list1)
    print(list2)
    truthlen = len(list1)
    overlap = 0
    list1ind = 0
    list2ind = 0
    if list1[list1ind] == None:
        truthlen = 0
        return None, truthlen
    if list1[list1ind] == None and list2[list2ind] == None:
        truthlen = 0
        return None, truthlen
    if list2[list2ind] == None:
        #print('no neighbs', list1, list2[list2ind])
        return int(0), truthlen

    while list1ind < len(list1) and list2ind<len(list2):
        if list1[list1ind] == list2[list2ind]:
            overlap += 1
            list1ind += 1
            list2ind += 1
        elif list1[list1ind] < list2[list2ind]:
            print('list1 less', list1[list1ind], list2[list2ind])
            list1ind += 1
        elif list1[list1ind] > list2[list2ind]:
            print('list2 less', list1[list1ind], list2[list2ind])
            list2ind += 1
    return overlap/truthlen, truthlen

    '''    
    for l1 in list1:
        for l2 in list2:
        if l2 == l1:
            print(l2, l1, 'equal')
        if l2>l1:
            print(l2,l1, 'l2 greater')
        if l2<l1:
            print(l2, l1, 'l2 less')
            '''





def parse_file(files):
    samples = {}
    for i in range(len(files)):
        file = files[i]
        with open(file) as f:
            for line in f:
                if line.startswith('\t'):
                    #print(line)
                    line = line.strip().split()
                    line.insert(0, '')
                    #print(line)
                    #print('line starts w tab')
                elif line == '\n':
                    #print(line)
                    line = line.strip().split()
                    #print(line)
                    line.insert(0, '')
                    #print(line)
                    #print('line starts with ""')
                else:
                    line = line.strip().split()
                    #print(line)

                    #move out of conditional 
                #print('line', line)
                if i == 0:
                    
                    if line[0] not in samples:
                        if len(line) > 1:
                            #print(line[1])
                            samples[line[0]] = [line[1].split(',')]
                        else:
                            samples[line[0]] = [[None]]
                        #if  not line[0].startswith('SRR'):
                        #    print('after',samples[line[0]])
                    else:
                        print('weird', line)
                
                if i == 1:
                
                    if line[0] in samples:
                    #if  not line[0].startswith('SRR'):
                            #print('before',samples[line[0]])
                        if len(line) > 1:
                            #print(line[1])
                            samples[line[0]].append(line[1].split(','))
                        else:
                            samples[line[0]].append([None])
                        #if  not line[0].startswith('SRR'):
                        #    print('after',samples[line[0]])    
                    else:
                        print('weird', line)
    #for s in samples:
    #    print(s, samples[s])
    return samples

#put in main function?

#x nearest
'''
one = fix_file(tsv1)
two = fix_file(tsv2)
keys1 = sorted(list(one.keys()))
keys2 = sorted(list(two.keys()))
assert keys1==keys2
keys=keys1
'''

#within x
files = [tsv1, tsv2]

samples = parse_file(files)




with open(out, 'w') as o:
    #within x
    keys = sorted(list(samples.keys()))
    for k in keys:
        #if you want the neighborhoods in a file, write the following to a file
        #neighborhood.write(k, samples[k][0], samples[k][1], sep = '\t')
        #print(k, samples[k][0], samples[k][1], sep = '\t')
        print(k)
        #within x
        result, truthlen = check_overlap(sorted(samples[k][0]), sorted(samples[k][1]))
        #jsim = jaccard_similarity(samples[k][0], samples[k][1])
        
        # x nearest
        #jsim = jaccard_similarity(one[k],two[k])
        
        #print(round(jsim,4))
        o.write(f'{k}\t{truthlen}\t{result}\n')
        
#what should we consider two empty lists? they are correct in that they are both empty 
#if i make the lists empty i get an error, if i fill them with none, they will give full value 
