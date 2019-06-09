from Bio import AlignIO

allID = set()
alignment_count = 0
seq_count = 0
for alignment_block in AlignIO.parse("/Users/liuyuxuan/Downloads/208dataset/dog_chrX.maf", "maf"):
    # print("printing a new multiple alignment")

    # print("with " + str(alignment_block.get_alignment_length()) + " sites")
    alignment_count += 1
    # print(alignment_count)
    for seqrec in alignment_block:
        # print("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" %
        #       (seqrec.annotations["start"],
        #        seqrec.annotations["strand"],
        #        seqrec.annotations["srcSize"],
        #        seqrec.annotations["size"]))
        allID.add(seqrec.id)
        seq_count += 1
        # print(seqrec.seq)
        # if not bool1:
        #     bool1 = True
        # else:
        #     break

    # print(alignment_block.format("fasta"))
    # if count < 0:
    #     count += 1
    # else:
    #     break
print(allID)
print(len(allID))
print(alignment_count)
print(seq_count)