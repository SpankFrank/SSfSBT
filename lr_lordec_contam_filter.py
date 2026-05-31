from argparse import ArgumentParser
from typing   import Iterable, Generator
from math     import inf

from file_services.utils import get_read_reader

class MyArgumentParser(ArgumentParser):

    prog        =   "lr_lordec_contam_filter"

    description =   """
                    A simple script facilitating the filtering contaminations from
                    long reads or transcripts corrected with LoRDEC using Kraken2-filtered
                    short reads.
                    """
    
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)

        self.add_argument("longreads")
        self.add_argument("-v",
                          metavar="",
                          help="Report progress every v processed long reads [default: 100.000]",
                          default=100_000)
        
        fo = self.add_argument_group("Filtering options")
        
        fo.add_argument("-mbc","--min_bases_corrected",
                        type=int,
                        default=21)
        fo.add_argument("mbu","--max_bases_uncorrected",
                        type=int,
                        default=inf)
        fo.add_argument("-mfc","--min_fraction_corrected",
                        type=float,
                        default=0.5)

def main():

    args     = MyArgumentParser().parse_args()
    fs       = get_read_reader(args.longreads)
    accepted = []
    rejected = []
    
    bases_total                     = 0
    bases_accepted                  = 0
    bases_accepted_corrected        = 0
    bases_accepted_uncorrected      = 0
    bases_rejected                  = 0
    
    print("Reading and processing ...")

    for i, longread in enumerate(fs.read(args.longreads)):
        
        if i % args.v == 0:
            print(f"{i:>8}")
            
        bases_total += longread["length"]
        
        read_bases_corrected    = len([c for c in longread["sequence"] if c.isupper()])
        read_bases_uncorrected  = longread["length"] - read_bases_corrected
        read_fraction_corrected = read_bases_corrected / longread["length"]
        
        accept = ((read_bases_corrected > args.min_bases_corrected) and
                  (read_bases_uncorrected < args.max_bases_uncorrected) and
                  (read_fraction_corrected > args.min_fraction_corrected))
        
        if accept:
            accepted.append(longread)
            
            bases_accepted              += longread["length"]
            bases_accepted_corrected    += read_bases_corrected
            bases_accepted_uncorrected  += read_bases_uncorrected
        else:
            rejected.append(longread)
            
            bases_rejected += longread["length"]
    
    print(f"{i:>8}")
    print("Writing ...")
    fs.write(args.longreads+".accepted", accepted)
    fs.write(args.longreads+".rejected", rejected)
    
    reads_accepted = len(accepted)
    reads_rejected = len(rejected)
    reads_total    = reads_accepted + reads_rejected

    print("##################################################################################")
    print(f"Reads accepted:     {reads_accepted:>10} ({(reads_accepted/reads_total)*100:.3f}%)")
    print(f"Reads rejected:     {reads_rejected:>10} ({(reads_rejected/reads_total)*100:.3f}%)")
    print()
    print(f"Bases accepted:     {bases_accepted:>10} ({(bases_accepted/bases_total)*100:.3f}%)")
    print(f"Bases rejected:     {bases_rejected:>10} ({(bases_rejected/bases_total)*100:.3f}%)")
    print()
    print("In accepted reads:")
    print(f"Bases corrected:    {bases_accepted_corrected:>10} ({(bases_accepted_corrected/bases_accepted)*100:.3f}%)")
    print(f"Bases uncorrected:  {bases_accepted_uncorrected:>10} ({(bases_accepted_uncorrected/bases_accepted)*100:.3f}%)")
    print()
    print("##############################################")
    print("#    Simon says: Thanks for using SSfSBT!    #")
    print("##############################################")

if __name__ == "__main_":

    main()
