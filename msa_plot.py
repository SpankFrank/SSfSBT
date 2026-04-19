from argparse       import ArgumentParser
from Bio.SeqIO      import parse
from Bio.SeqRecord  import SeqRecord
from collections    import Counter
from json           import loads
from matplotlib     import pyplot, axes
from numpy          import array, nan
from pprint         import pprint
from scipy.stats    import entropy
from seaborn        import heatmap

class MyArgumentParser(ArgumentParser):

    prog        =   "msa_plot"

    description =   """
                    Vizualisation of multiple sequence alignments
                    """
    
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)

        self.add_argument("msa",
                          type      = str,
                          help      = "Multiple sequence alignment (FASTA-format)")
        
        self.add_argument("-s", "--start",
                          metavar   = "",
                          type      = int,
                          default   = 1,
                          help      = "Start of the section you want to show (default: 1)")
        self.add_argument("-e", "--end",
                          metavar   = "",
                          type      = int,
                          default   = 999_999,
                          help      = "End of the section you want to show (default: alignment end)")
        self.add_argument("-r", "--row_length",
                          metavar   = "",
                          type      = int,
                          default   = 100,
                          help      = "Number of alignment columns to show per row (default: 100)")
        self.add_argument("-fs", "--font_style",
                          metavar   = "",
                          type      = str,
                          default   ="normal",
                          help      = "Font style for sequence names ( default: \'normal\", choices: [\'normal\', \'italic\', \'oblique\'])")
        self.add_argument("-cmap", "--color_map",
                          metavar   = "",
                          type      = str,
                          default   = "tab20",
                          help      = "Choose one of the many colormaps available to matplotlib (default: \'tab20\' )")
        self.add_argument("-cval", "--color_value",
                          metavar   = "",
                          type      = str,
                          help      = "A character to value map (JSON-format)")
        self.add_argument("-hl", "--highlight_last",
                          action="store_true",
                          help="Highlight last sequence by witing it's name in bold letters")
        self.add_argument("-of", "--output_format",
                          metavar   = "",
                          type      = str,
                          default   = "png",
                          help      = "Output format (default: \'png\' choices: ['eps', 'jpg', 'jpeg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'])")
        
def load_msa(file_path: str) -> list[SeqRecord]:
    
    return list(parse(file_path, "fasta"))

def load_color_values(file_path: str | None) -> dict[str, float]:
    
    if file_path is None:
        color_values = dict({c: float(i) for i, c in enumerate("ACDEFGHIKLMNPQRSTVWYX")}) 
    else:
        color_values = loads(file_path)
        
    color_values["-"] = nan
    
    return color_values
    
def column_conservation(column: list[str]) -> float:
    """Computes the conservation level at a position in the msa

    Args:
        column (list[str]): A list of single-character strings representing a column of a multiple sequence alignment

    Returns:
        float: conservation level (ranges from 0 to 1)
    """

    counts          = Counter(column)                                       # a dictionary that holds the frequency of each aminoacid
    
    if counts["-"] == len(column):
        
        return 0
    else:
        
        probabilities   = [count/len(column) for count in counts.values()]  # the frequency of each aminoacid 
        max_entropy     = entropy([1/len(column) for _ in column])          # the theoretical maximum shannon entropy (if each aminoacid was different)
        this_entropy    = entropy(probabilities)                            # the actual entropy in this msa column
        
        return float((max_entropy - this_entropy) / max_entropy)
    
def plot_conservation(msa_section: list[list[str]], axes: axes._axes.Axes) -> None:
    
    conservation = [column_conservation([row[j] for row in msa_section])
                    for j in range(len(msa_section[0]))]
    
    conservation = [[j/256 if j/256<c else nan for _, c in enumerate(conservation)]
                    for j in reversed(range(0,256))]

    heatmap(array(conservation), cmap="Greys", cbar=False, ax=axes)
    
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xticklabels([])
    axes.set_yticklabels([])
    axes.set_ylabel("Conservation", rotation=0, ha="right")
    
def plot_alignment(msa_section:     list[str],
                   axes:            axes._axes.Axes,
                   names:           list[str],
                   section_start:   int,
                   section_end:     int,
                   end:             int,
                   color_map:       str,
                   color_value:     dict[str, float],
                   font_style:      str):
    
    section_colors = [[color_value[c] for c in sequence] for sequence in msa_section]
    x              = range(min(section_end, end)-section_start)
    xticks         = [i+0.5 for i in x if i%10==0]
    xticklabels    = [str(i+section_start+1) for i in x if i%10==0]
    
    heatmap(section_colors, ax=axes, cmap=color_map, annot=msa_section, fmt='', cbar=False)
    
    axes.set_xticks(xticks)
    axes.set_xticklabels(xticklabels)
    axes.set_yticks([0.5+i for i in range(len(msa_section))])
    axes.set_yticklabels([name.replace("_"," ") for name in names], rotation=0, fontstyle=font_style)
    
def main():
    
    args        = MyArgumentParser().parse_args()
    
    print("1. Loading alignment ...")
    msa         = load_msa(args.msa)
    print("\tDone!")
    
    fig, axes = pyplot.subplots()
    
    color_value   = load_color_values(args.color_value)
    names         = [entry.name.replace("_", " ") for entry in msa]
    start         = max(args.start, 1) - 1
    end           = min(len(msa[0]), args.end) - 1
    row_length    = min(args.row_length, end-start)
    n_rows        = int((min(end,len(msa[0])) - start-1) / row_length) + 1
    section_start = start
    section_end   = min(section_start+args.row_length,args.end) - 1
    
    print("2. Calculating plot parameters")
    print(f"\tStart:      {start+1}")
    print(f"\tEnd:        {end+1}")
    print(f"\tRow length: {row_length}")
    print(f"\tRows:       {n_rows}")
    
    # Plotting
    # Longer multiple sequence alignments are divided into several sections and displayed
    # one below the other. There are three subplots for each section: 
    # a) Showing the conservation
    # b) Showing the alignment
    # c) An empty place-holder plot
    
    print("3. Plotting ...")
    
    fig, axes = pyplot.subplots(3*n_rows, 1,
                                figsize       = (15,n_rows*5),
                                height_ratios = [n_rows if i%3==1 else 1 for i in range(3*n_rows)])
    
    for i in range(n_rows):
        
        print(f"\tProcessing section {i+1} / {n_rows} (Positions {section_start+1:>4}-{min(section_end,args.end)+1:>4})")
        
        msa_section   = [[c for c in entry.seq[section_start:section_end+1]] for entry in msa]
        
        if len(msa_section[0]) < row_length:
            for sequence in msa_section:
                sequence += ["-" for _ in range(row_length-len(msa_section[-1])-1)]
                
        plot_conservation(msa_section, axes[3*i])
        plot_alignment(msa_section,
                       axes[3*i+1],
                       names,
                       section_start,
                       section_end,
                       args.end,
                       args.color_map,
                       color_value,
                       args.font_style)
        if args.highlight_last:
            axes[3*i+1].get_yticklabels()[-1].set_fontweight("bold")
        axes[3*i+2].axis("off")
        
        section_start = section_end + 1
        section_end   = min(section_end + row_length, args.end)

    print("\tDone!")
    print("4. Saving ...")
        
    pyplot.savefig(f"{args.msa}.ssfsbt_msa.{args.output_format}", bbox_inches='tight')
    
    print("\tDone!\n")
    print("##############################################")
    print("#    Simon says: Thanks for using SSfSBT!    #")
    print("##############################################")
    
if __name__ == "__main__":
    main()