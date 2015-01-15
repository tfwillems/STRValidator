import matplotlib as mpl
mpl.use('Agg')

import collections
import sys
import numpy
import matplotlib.pyplot as plt
import vcf
from matplotlib.backends.backend_pdf import PdfPages
from fractions import Fraction

class TRIO:
    def __init__(self, child, mother, father):
        self.child  = child
        self.mother = mother
        self.father = father
    def __str__(self):
        return "%s\t%s\t%s"%(self.child, self.mother, self.father)

class FATHER_SON_PAIR:
    def __init__(self, son, father):
        self.son    = son
        self.father = father
    def __str__(self):
        return "%s\t%s"%(self.son, self.father)

def read_1kg_pedigree_file(input_file, header=True):
    data = open(input_file, "r")
    if header:
        data.readline()

    trios, father_son_pairs = [], []
    for line in data:
        tokens = line.strip().split()
        if tokens[2] != "0" and tokens[3] != "0":
            child, dad, mom = tokens[1:4]
            trios.append(TRIO(child,  dad, mom))
        if tokens[2] != "0" and tokens[4] == "1":
            father_son_pairs.append(FATHER_SON_PAIR(tokens[1], tokens[2]))
    data.close()
    print("There are %d trios and %d father-son-pairs in the pedigree file"%(len(trios), len(father_son_pairs)))
    return trios, father_son_pairs

# Find the index for the highest bin which is less than
# or equal to the provided value
def find_index(bins, value):
    low  = 0
    high = len(bins)-1
    while high > low + 1:
        midval = bins[(low+high)/2]    
        if value > midval:
            low  = (low+high)/2 
        elif value < midval:
            high = (low+high)/2 - 1
        else:
            return (low+high)/2
    if value < bins[low]:
        exit("Unable to find index. Exiting...")
    if value >= bins[high]:
        return high
    else:
        return low

def is_discordant(a11, a12, a21, a22):
    if (a11 == a21 and a12 == a22) or (a11 == a22 and a12 == a21):
        return False
    else:
        return True
 
def is_mendelian(a11, a12, a21, a22, a31, a32):
    if (a31 == a11 or a31 == a12) and (a32 == a21 or a32 == a22):
        return True
    elif (a31 == a21 or a31 == a22) and (a32 == a11 or a32 == a12):
        return True
    else:
        return False

def draw_bp_histogram(discordant_counts, pdfpage):
    # Create histogram of father-son differences                                                                                                                                 
    bp_diff_counts     = [collections.defaultdict(int) for _ in xrange(6)]
    repeat_diff_counts = [collections.defaultdict(int) for _ in xrange(6)]
    out_frame_count    = 0
    in_frame_count     = 0
    for key,val in discordant_counts.items():
        bp_diff_counts[key[2]-1][key[1]-key[0]] += val
        repeat_diff_counts[key[2]-1][Fraction(key[1]-key[0], key[2])] += val
    for xlabel,diff_counts,in_frame in zip(["bps", "repeats"],
                                           [bp_diff_counts, repeat_diff_counts],
                                           [lambda bp,period: bp%period == 0, lambda rep,period: int(rep)==float(rep) ]):
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        diffs = sorted(list(set(reduce(lambda x,y:x+y, map(lambda z: z.keys(), diff_counts)))))
        colors = ['c', 'r', 'g', 'y', 'b', 'm']
        heights = numpy.zeros(len(diffs))
        for i in xrange(6):
            vals = [diff_counts[i][x] for x in diffs]
            if sum(vals) == 0:
                continue

            in_frame_trips  = filter(lambda x: in_frame(x[0], i+1), zip(diffs, vals, heights))
            out_frame_trips = filter(lambda x: not in_frame(x[0], i+1), zip(diffs, vals, heights))
            if len(in_frame_trips) != 0:
                x,y,h = zip(*in_frame_trips)
                in_frame_count += sum(y)
                ax.bar(x, y, bottom=h, align='center', color=colors[i], width=0.25, label=str(i+1))
            if len(out_frame_trips) != 0:
                x,y,h = zip(*out_frame_trips)
                out_frame_count += sum(y)
                ax.bar(x, y, bottom=h, align='center', color=colors[i], width=0.25, label=str(i+1), hatch='//')
            heights += vals
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xlabel(r"$father-son ("+xlabel+")$")
        ax.set_ylabel(r"$n_{calls}$")
        ax.legend()
        pdfpage.savefig(fig)
    print("IN FRAME=%d, OUT FRAME=%d"%(in_frame_count/2, out_frame_count/2))


class CHRY_STATS:
    def __init__(self, father_son_pairs, call_output):
        self.pairs        = father_son_pairs
        self.output_calls = open(call_output, "w")

    def initialize(self, vcf_reader):
        sample_indices    = dict(zip(vcf_reader.samples, range(len(vcf_reader.samples))))
        self.pair_indices = []
        for i in xrange(len(self.pairs)):
            if self.pairs[i].son not in sample_indices:
                exit("Unable to assess chrY inheritance because no data was found for " + self.pairs[i].son)
            if self.pairs[i].father not in sample_indices:
                exit("Unable to assess chrY inheritance because no data was found for " + self.pairs[i].father)
            self.pair_indices.append([sample_indices[self.pairs[i].father], sample_indices[self.pairs[i].son]])
        self.missing_data_skip_counts = numpy.zeros(len(self.pair_indices))
        self.het_gt_skip_counts       = numpy.zeros(len(self.pair_indices))

        self.num_concordant    = 0
        self.num_discordant    = 0
        self.pair_info         = {}
        self.discordant_counts = collections.defaultdict(int)
        self.call_count        = 0

    def process_record(self, record):
        motif_len = len(record.INFO['MOTIF'])
        for i in xrange(len(self.pair_indices)):
            if any(map(lambda x: record.samples[x]['GT'] is None, self.pair_indices[i])):
                self.missing_data_skip_counts[i] += 1
                continue

            self.call_count += 1
            father = record.samples[self.pair_indices[i][0]]
            son    = record.samples[self.pair_indices[i][1]]
            gb_1a, gb_1b = map(int, father['GB'].split("/"))
            gb_2a, gb_2b = map(int, son['GB'].split("/"))
            self.output_calls.write("%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n"%(self.call_count, record.CHROM, record.POS, record.INFO['END'], 
                                                                                        gb_1a + gb_1b, gb_2a + gb_2b,
                                                                                        gb_1a,  gb_1b, gb_2a,  gb_2b, father.sample, son.sample))
            
            if gb_1a != gb_1b or gb_2a != gb_2b:
                self.het_gt_skip_counts[i] += 1
                if gb_1a != gb_1b:
                    print("chrY\t%d\t%d\t%s\t%s\t%s"%(record.POS, record.INFO["END"], father.sample, str(gb_1a) + "|" + str(gb_1b), "HET"))
                if gb_2a != gb_2b:
                    print("chrY\t%d\t%d\t%s\t%s\t%s"%(record.POS, record.INFO["END"], father.sample, str(gb_2a) + "|" + str(gb_2b), "HET"))
                continue

            if gb_1a != gb_2a:
                self.num_discordant += 1
                self.discordant_counts[(gb_1a, gb_2a, motif_len)] +=1
                print("chrY\t%d\t%d\t%s\t%s\t%s"%(record.POS, record.INFO["END"], 
                                                  father.sample + "," + son.sample,
                                                  str(gb_1a) + "," + str(gb_2b), "DISCORDANT"))
            else:
                self.num_concordant += 1
            if (gb_1a, gb_2a) not in self.pair_info:
                self.pair_info[(gb_1a, gb_2a)] = []
            self.pair_info[(gb_1a, gb_2a)].append((record.CHROM, record.POS, record.INFO['END'], father.sample+"-"+son.sample))

    def finish(self, pdfpage, output_prefix):
        print("WARNING: Skipped " + str(self.missing_data_skip_counts) + " comparisons due to missing data for one or more individuals")
        print("WARNING: Skipped " + str(self.het_gt_skip_counts)       + " comparisons due to heterozygous genotypes for one or more individuals")

        if self.num_discordant + self.num_concordant != 0:
            print("%d vs. %d = %f Percent"%(self.num_discordant, self.num_concordant, 100.0*self.num_discordant/(self.num_discordant+self.num_concordant)))
        else:
            print("WARNING: No chrY calls were applicable for comparison")

        # Create bubble plot using all data
        fig  = plt.figure()
        ax   = fig.add_subplot(111)
        x, y = zip(*self.pair_info.keys())
        s    = numpy.array(map(len, self.pair_info.values()))*10
        ax.scatter(x, y, s=s, alpha=0.7)
        ax.set_xlabel("Father's genotype (bp)")
        ax.set_ylabel("Son's genotype (bp)")
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.plot(numpy.arange(min(x)-5, max(x)+5, 1.0), numpy.arange(min(y)-5, max(y)+5, 1.0), linestyle='--', color='k')
        pdfpage.savefig(fig)
        
        # Create histogram of father-son differences
        draw_bp_histogram(self.discordant_counts, pdfpage)

        viz_output = open(output_prefix+"_chrY.csv", "w")
        viz_output.write(",".join(["X","Y", "CHROMS", "STARTS", "STOPS", "SAMPLES"]) + "\n")
        for key,val in self.pair_info.items():
            chroms, positions, ends, samples = map(list, zip(*val))
            viz_output.write(",".join([str(key[0]), str(key[1]), "_".join(chroms), "_".join(map(str, positions)), "_".join(map(str, ends)), "_".join(map(str, samples))]) + "\n")
        viz_output.close()

        self.output_calls.close()


class MENDELIAN_STATS:
    def __init__(self, trios, coverage_bins, quality_bins, max_coverage, quality_thresholds):
        self.trios         = trios
        self.coverage_bins = coverage_bins
        self.quality_bins  = quality_bins
        self.max_coverage  = max_coverage
        self.qual_thresh   = quality_thresholds

    def initialize(self, vcf_reader):
        sample_indices    = dict(zip(vcf_reader.samples, range(len(vcf_reader.samples))))
        self.trio_indices = []
        for i in xrange(len(self.trios)):
            if self.trios[i].child not in sample_indices:
                exit("Unable to calculate Mendelian inheritance because no data was found for " + self.trios[i].child)
            if self.trios[i].father not in sample_indices:
                exit("Unable to calculate Mendelian inheritance because no data was found for " + self.trios[i].father)
            if self.trios[i].mother not in sample_indices:
                exit("Unable to calculate Mendelian inheritance because no data was found for " + self.trios[i].mother)
            
            # Father, Mother, Child
            self.trio_indices.append(map(lambda x: sample_indices[x], [self.trios[i].father, self.trios[i].mother, self.trios[i].child]))

        self.coverage_bins = numpy.concatenate(([-100000], self.coverage_bins))
        self.quality_bins  = numpy.concatenate(([-100000], self.quality_bins))
        
        # Quality/Coverage x Trios x Period x Thresholds
        self.all_loci_nstrs  = [numpy.zeros((len(self.trios), 5, len(self.coverage_bins))), numpy.zeros((len(self.trios), 5, len(self.quality_bins)))]
        self.all_loci_nmend  = [numpy.zeros((len(self.trios), 5, len(self.coverage_bins))), numpy.zeros((len(self.trios), 5, len(self.quality_bins)))]
        self.disc_loci_nstrs = [numpy.zeros((len(self.trios), 5, len(self.coverage_bins))), numpy.zeros((len(self.trios), 5, len(self.quality_bins)))]
        self.disc_loci_nmend = [numpy.zeros((len(self.trios), 5, len(self.coverage_bins))), numpy.zeros((len(self.trios), 5, len(self.quality_bins)))]

        self.missing_data_skip_counts = numpy.zeros(len(self.trios))
        self.coverage_skip_counts     = numpy.zeros(len(self.trios))

        # Trios x Period x Thresholds
        self.all_loci_nstrs_min_q  = numpy.zeros((len(self.trios), 5, len(self.coverage_bins)))
        self.all_loci_nmend_min_q  = numpy.zeros((len(self.trios), 5, len(self.coverage_bins)))
        self.disc_loci_nstrs_min_q = numpy.zeros((len(self.trios), 5, len(self.coverage_bins)))
        self.disc_loci_nmend_min_q = numpy.zeros((len(self.trios), 5, len(self.coverage_bins)))
        
    def process_record(self, record):
        for i in xrange(len(self.trios)):
            if any(map(lambda x: record.samples[x]['GT'] is None, self.trio_indices[i])):
                self.missing_data_skip_counts[i] += 1
                continue            
            if 'X' in record.CHROM or 'x' in record.CHROM or 'Y' in record.CHROM or 'y' in record.CHROM:
                continue
            
            q1, q2, q3 = map(lambda x: record.samples[x]["Q"],  self.trio_indices[i])
            c1, c2, c3 = map(lambda x: record.samples[x]["DP"], self.trio_indices[i])
            a11, a21   = record.samples[self.trio_indices[i][0]]["GT"].split("/")
            a21, a22   = record.samples[self.trio_indices[i][1]]["GT"].split("/")
            a31, a32   = record.samples[self.trio_indices[i][2]]["GT"].split("/")
            discordant = is_discordant(a11, a12, a21, a22)
            mendelian  = is_mendelian(a11, a12, a21, a22, a31, a32)

            # Filter out loci with too high of coverage
            if max(c1, c2, c3) > self.max_coverage:
                self.coverage_skip_counts[i] += 1
                continue

            coverage  = min(c1, c2, c3)
            bin_idx   = find_index(self.coverage_bins, coverage)
            motif_len = len(record.INFO["MOTIF"])-2
            self.all_loci_nstrs [0][i][motif_len][bin_idx] += 1
            self.all_loci_nmend [0][i][motif_len][bin_idx] += mendelian*1
            self.disc_loci_nstrs[0][i][motif_len][bin_idx] += discordant*1
            self.disc_loci_nmend[0][i][motif_len][bin_idx] += discordant*mendelian*1

            quality   = min(q1, q2, q3)
            bin_idx   = find_index(self.quality_bins, quality)
            self.all_loci_nstrs [1][i][motif_len][bin_idx] += 1
            self.all_loci_nmend [1][i][motif_len][bin_idx] += mendelian*1
            self.disc_loci_nstrs[1][i][motif_len][bin_idx] += discordant*1
            self.disc_loci_nmend[1][i][motif_len][bin_idx] += discordant*mendelian*1

            coverage  = min(c1, c2, c3)
            bin_idx   = find_index(self.coverage_bins, coverage)
            if quality > self.qual_thresh[motif_len]:
                self.all_loci_nstrs_min_q  [i][motif_len][bin_idx] += 1
                self.all_loci_nmend_min_q  [i][motif_len][bin_idx] += mendelian*1
                self.disc_loci_nstrs_min_q [i][motif_len][bin_idx] += discordant*1
                self.disc_loci_nmend_min_q [i][motif_len][bin_idx] += discordant*mendelian*1

    def finish(self, pdfpage):
        print("WARNING: Skipped " + str(self.missing_data_skip_counts) + " loci due to missing data for one or more individual")
        print("WARNING: Skipped " + str(self.coverage_skip_counts)     + " loci due to too high coverage")

        # Iterate over coverage and quality stats
        types = ['Coverage', 'Quality', 'Coverage']
        bins  = [self.coverage_bins, self.quality_bins, self.coverage_bins]
        for n in xrange(3):
            # Sum across all trios
            if n == 0 or n == 1:
                all_loci_nstrs  = numpy.sum(self.all_loci_nstrs [n],  axis=0)
                all_loci_nmend  = numpy.sum(self.all_loci_nmend [n],  axis=0)
                disc_loci_nstrs = numpy.sum(self.disc_loci_nstrs[n],  axis=0)
                disc_loci_nmend = numpy.sum(self.disc_loci_nmend[n],  axis=0)
            else:
                all_loci_nstrs  = numpy.sum(self.all_loci_nstrs_min_q,  axis=0)
                all_loci_nmend  = numpy.sum(self.all_loci_nmend_min_q,  axis=0)
                disc_loci_nstrs = numpy.sum(self.disc_loci_nstrs_min_q, axis=0)
                disc_loci_nmend = numpy.sum(self.disc_loci_nmend_min_q, axis=0)

            # Create plots for individual periods
            fig  = plt.figure()
            ax1 = fig.add_subplot(221)
            ax1.set_ylabel("Fraction Mendelian")
            ax1.set_title("All sites")

            ax2 = fig.add_subplot(222, sharey=ax1)
            ax2.set_title("Discordant parental sites")

            ax3 = fig.add_subplot(223, sharex=ax1)
            ax3.set_xlabel(types[n] + " threshold")
            ax3.set_ylabel("# genotypes")
            ax3.set_yscale('log')

            ax4 = fig.add_subplot(224, sharex=ax2, sharey=ax3)
            ax4.set_xlabel(types[n] + " threshold")
            ax4.set_yscale('log')

            box1 = ax1.get_position()
            ax1.set_position([box1.x0, box1.y0, box1.width*0.9, box1.height])
            ax2.set_position([box1.x0 + box1.width*1.15, box1.y0, box1.width*0.9, box1.height])
            box3 = ax3.get_position()
            ax3.set_position([box3.x0, box3.y0, box3.width*0.9, box3.height])
            ax4.set_position([box3.x0 + box3.width*1.15, box3.y0, box3.width*0.9, box3.height])
            

            font_size = 9
            for i in xrange(5):
                nstrs_all  = numpy.cumsum(all_loci_nstrs [i][::-1])[::-1]
                nmend_all  = numpy.cumsum(all_loci_nmend [i][::-1])[::-1]
                nstrs_disc = numpy.cumsum(disc_loci_nstrs[i][::-1])[::-1]
                nmend_disc = numpy.cumsum(disc_loci_nmend[i][::-1])[::-1]
                all_fracs  = (1.0*nmend_all/nstrs_all)[1:]
                disc_fracs = (1.0*nmend_disc/nstrs_disc)[1:]
                ax1.plot(bins[n][1:], all_fracs,      '-o', label=str(i+1))
                ax2.plot(bins[n][1:], disc_fracs,     '-o', label=str(i+1))
                ax3.plot(bins[n][1:], nstrs_all[1:],  '-o', label=str(i+1))
                ax4.plot(bins[n][1:], nstrs_disc[1:], '-o', label=str(i+1))
            ax4.legend(bbox_to_anchor=(1.05, 0.9, 0.25, 0.2), loc='center left')
            for ax in [ax1, ax2, ax3, ax4]:
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(font_size)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(font_size)
            pdfpage.savefig(fig)


            # Create plots using all periods

            # Sum across all periods
            all_loci_nstrs  = numpy.sum(all_loci_nstrs,  axis=0)
            all_loci_nmend  = numpy.sum(all_loci_nmend,  axis=0)
            disc_loci_nstrs = numpy.sum(disc_loci_nstrs, axis=0)
            disc_loci_nmend = numpy.sum(disc_loci_nmend, axis=0)
        
            # Transform into running sums
            all_loci_nstrs  = numpy.cumsum(all_loci_nstrs[::-1])[::-1]
            all_loci_nmend  = numpy.cumsum(all_loci_nmend[::-1])[::-1]
            disc_loci_nstrs = numpy.cumsum(disc_loci_nstrs[::-1])[::-1]
            disc_loci_nmend = numpy.cumsum(disc_loci_nmend[::-1])[::-1]
            
            # Calculate the fraction of Mendelian inheritance for all loci and discordant loci
            all_loci_fracs  = (1.0*all_loci_nmend/all_loci_nstrs)[1:]
            disc_loci_fracs = (1.0*disc_loci_nmend/disc_loci_nstrs)[1:]
        
            fig = plt.figure()
            ax1  = fig.add_subplot(221)
            ax1.set_ylabel("Fraction Mendelian")
            ax1.set_title("All sites")
            ax1.plot(bins[n][1:], all_loci_fracs, '-o')

            ax2  = fig.add_subplot(222, sharey=ax1)
            ax2.plot(bins[n][1:], disc_loci_fracs, '-o')
            ax2.set_title("Discordant parental sites")
           
            ax3  = fig.add_subplot(223, sharex=ax1)
            ax3.set_xlabel(types[n] + " threshold")
            ax3.set_ylabel("# genotypes")
            ax3.set_yscale('log')
            ax3.plot(bins[n][1:], all_loci_nstrs[1:], '-o')
           
            ax4  = fig.add_subplot(224, sharex=ax2, sharey=ax3)
            ax4.set_xlabel(types[n] + " threshold")
            ax4.set_yscale('log')
            ax4.plot(bins[n][1:], disc_loci_nstrs[1:], '-o')

            for ax in [ax1, ax2, ax3, ax4]:
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(font_size)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(font_size)
            pdfpage.savefig(fig)


            mpl.rcParams['xtick.labelsize'] = 10
            mpl.rcParams['ytick.labelsize'] = 10
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(bins[n][1:], all_loci_fracs, '-o', color='b')
            ax1.set_ylabel("Fraction Mendelian")
            ax1.set_xlabel(types[n] + " threshold")
            ax2 = ax1.twinx()
            ax2.set_yscale('log')
            ax2.plot(bins[n][1:], all_loci_nstrs[1:], '-o', color='g')
            pdfpage.savefig(fig)
            ax1.axis('equal')
            pdfpage.savefig(fig)
            
            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            ax1.plot(bins[n][1:], all_loci_fracs, '-o', color='b')
            ax1.set_ylabel("Fraction Mendelian")
            ax1.xaxis.set_ticks_position('bottom')
            ax1.yaxis.set_ticks_position('left')
            ax2.set_xlabel(types[n] + " threshold")
            ax2.plot(bins[n][1:], all_loci_nstrs[1:], '-o', color='g')
            ax2.set_yscale('log')
            ax2.set_ylabel("# Called loci across trios")
            ax2.xaxis.set_ticks_position('bottom')
            ax2.yaxis.set_ticks_position('left')
            pdfpage.savefig(fig)

def main():
    print("Invocation syntax: python pedigree_analysis.py 1kg_pedigree_file.txt vcf_file.vcf output_file.pdf") 
    trios, father_son_pairs = read_1kg_pedigree_file(sys.argv[1], header=True)
    vcf_reader = vcf.Reader(filename=sys.argv[2])
    call_stats = sys.argv[3]
    samples    = vcf_reader.samples
    trios_with_data = []
    pairs_with_data = []

    for trio in trios:
        if trio.child in samples and trio.mother in samples and trio.father in samples:
            trios_with_data.append(trio)
    print("There are %d trios with data"%len(trios_with_data))
    
    for pair in father_son_pairs:
        if pair.father in samples and pair.son in samples:
            pairs_with_data.append(pair)
    print("There are %d father-son pairs with data"%(len(pairs_with_data)))

    coverage_bins  = numpy.append(numpy.arange(1.001, 5.0011, 1.0), numpy.arange(6.001, 18.0011, 2.0))
    quality_bins   = numpy.arange(0.0, 1.0, 0.1)
    quality_thresh = [0.9, 0.5, 0.5, 0.5, 0.5, 0.5]
    max_coverage   = 100
    processors     = [CHRY_STATS(pairs_with_data, call_stats)]
    #mend_stats     = MENDELIAN_STATS(trios_with_data, coverage_bins, quality_bins, max_coverage, quality_thresh)
    for proc in processors:
        proc.initialize(vcf_reader)
    for record in vcf_reader:
        for proc in processors:
            proc.process_record(record)
    pp = PdfPages(sys.argv[3]+".pdf")
    for proc in processors:
        proc.finish(pp, sys.argv[3])
    pp.close()
    return 0


if __name__ == "__main__":
    main()
