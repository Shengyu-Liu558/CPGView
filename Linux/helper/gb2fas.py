import sys
from Bio import SeqIO

gbk_filename = sys.argv[1] 
faa_filename = sys.argv[2]
#gbk_filename = "NC_005213.gbk"
#faa_filename = "NC_005213_converted.fna"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print ("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()
