from argparse   import ArgumentParser
from os         import path
from polars     import read_csv

import logging
import sys

class MyArgumentParser(ArgumentParser):

    prog        =   "kallisto2nanosim"

    description =   """
                    Converting Kallisto transcript abundance files for NanoSim.
                    """
    
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)

        self.add_argument("kallisto_file")
        self.add_argument("nanosim_file")
        self.add_argument("-r","--remove_underscore",
                          action='store_true')

    def parse_args(self):

        self.args = super().parse_args()

        if not path.isfile(self.args.kallisto_file):
            raise Exception("Kallisto file does not exist")
        
        if path.isfile(self.args.nanosim_file):
            raise Exception("NanoSim file exists")
        
        return self.args
        
def check_target_ids(dataframe):

    for id in dataframe["target_id"]:
        if "_" in id:
            w = ["Target IDs contain \"_\" as character.",
                 "NanoSim has/had a bug that might cause it to break because of this.",
                 "Consider using the --remove-underscore parameter,",
                 "or check if the bug has been removed"]
            logging.warning("\n".join(w))
        
def main():

    logging.basicConfig(level    = logging.INFO,
                        format   = "%(asctime)s %(levelname)s %(message)s",
                        datefmt  = "%d-%m-%Y %H:%M:%S",
                        handlers = [logging.StreamHandler(stream=sys.stdout)]
                        )

    args = MyArgumentParser().parse_args()

    df = read_csv(args.kallisto_file, separator="\t")
    df = df.select(["target_id", "est_counts", "tpm"])

    if args.remove_underscore != None:
        df = df.with_columns(df["target_id"].str.replace("_", "").alias("target_id"))
    else:
        check_target_ids(df)

    df.write_csv(args.nanosim_file, separator="\t")

    print("##############################################")
    print("#    Simon says: Thanks for using SSfSBT!    #")
    print("##############################################")

if __name__ == "__main__":
    main()

